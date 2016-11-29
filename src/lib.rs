#![feature(conservative_impl_trait)]
//! A fast parser for fastq.
//!
//! The records in a file can be accessed through three different interfaces:
//!
//! - `Parser::each`: This function takes a closure that is executed for each
//!   fastq record. It is the fastest way to iterate over the records, since
//!   no copying of the records is needed and we don't allocate anything during
//!   parsing.
//! - `Parser::record_sets`. This function returns an iterator over record sets.
//!   All records in a record set share the same data array, so we only need
//!   one allocation per record set.
//! - `Parser::parallel_each`. This is a conveniance function that wraps
//!   `Parser::record_sets` but passes the record sets to a number of
//!   background threads. For each of the record sets a closure is executed.
//!   Results from the threads are collected and returned to the caller.
//!
//! Since fastq file are usually compressed, this crate also includes a function
//! `thread_reader` to offload the decompression to a different core.
//!
//! # Examples
//!
//! A minimal program that reads uncompressed fastq from stdin and counts the
//! number of fastq records. Since we do not need ownership of the records we
//! can use the fastest `Parser::each`.
//!
//! ```rust,no_run
//! use std::io::stdin;
//! use fastq::Parser;
//!
//! let mut parser = Parser::new(stdin());
//! let mut total: usize = 0;
//! parser.each(|_| {
//!     total += 1
//! }).expect("Invalid fastq file");
//! println!("{}", total);
//! ```
//!
//! If the file is already in the system cache, this runs on my laptop about
//! 10% slower than `wc -l` (at about 2GB/s), and at disk speed (~500MB/s) otherwise.
//!
//! If we add lz4 decompression in a background thread:
//!
//! ```rust,no_run
//! extern crate lz4;
//! extern crate fastq;
//!
//! use std::io::stdin;
//! use fastq::{Parser, thread_reader};
//!
//! # fn main() {
//! // lz4 decompression is faster with a large 4MB buffer
//! let BUFSIZE: usize = 1<<22;
//! let QUEUELEN: usize = 2;
//!
//! let reader = lz4::Decoder::new(stdin()).expect("Input is not lz4 compressed");
//! let total = thread_reader(BUFSIZE, QUEUELEN, reader, |reader| {
//!     let mut total: u64 = 0;
//!     let mut parser = Parser::new(reader);
//!     parser.each(|_| {
//!         total += 1;
//!      }).expect("Invalid fastq file");
//!     total
//! }).expect("Reader thread paniced");
//! println!("{}", total);
//! # }
//! ```
//!
//! This gets us up to ~1.5GB/s on a fresh page cache.
//!
//! If we want to do more than just count the number of records (in this
//! example count how many sequences align to ATTAATCCAT with a score
//! better than 7), we probably wantto use more cores:
//!
//! ```rust,no_run
//! # extern crate lz4;
//! # extern crate fastq;
//! extern crate parasailors;
//!
//! # use std::io::stdin;
//! use fastq::{Parser, thread_reader, Record};
//! use parasailors as align;
//!
//! # fn main() {
//! # let BUFSIZE: usize = 1<<22;
//! # let QUEUELEN: usize = 2;
//! # let reader = lz4::Decoder::new(stdin()).expect("Input is not lz4 compressed");
//! const N_THREADS: usize = 2;
//!
//! let results = thread_reader(BUFSIZE, 3, reader, |reader| {
//!     let parser = Parser::new(reader);
//!     let results: Vec<u64> = parser.parallel_each(N_THREADS, |record_sets| {
//!         let matrix = align::Matrix::new(align::MatrixType::Identity);
//!         let profile = align::Profile::new(b"ATTAATCCAT", &matrix);
//!         let mut thread_total: u64 = 0;
//!         for record_set in record_sets {
//!             for record in record_set.iter() {
//!                 let score = align::local_alignment_score(&profile, record.seq(), 11, 1);
//!                 if score > 7 {
//!                     thread_total += 1;
//!                 }
//!             }
//!         }
//!         thread_total
//!     }).expect("Invalid fastq file");
//!     results
//! }).expect("Reader thread paniced");
//! println!("{}", results.iter().sum::<u64>());
//! # }
//! ```
//!
//! On my feeeble 2 core laptop this ends up being bound by the alignment at ~300MB/s,
//! but it should scale well to a larger number of cores.

use std::io::{Write, Result, Read, Error, ErrorKind};
use std::thread;
use std::sync::mpsc::{SyncSender, sync_channel};
use std::sync::Arc;
use std::iter::FromIterator;
use memchr::memchr;

extern crate memchr;

mod thread_reader;
mod buffer;

pub use thread_reader::thread_reader;

//const BUFSIZE: usize = 17 * 1024;
const BUFSIZE: usize = 68 * 1024;


/// Trait to be implemented by types that represent fastq records.
pub trait Record {
    /// Return the fastq sequence as byte slice
    fn seq(&self) -> &[u8];
    /// Return the id-line of the record as byte slice
    fn head(&self) -> &[u8];
    /// Return the quality of the bases as byte slice
    fn qual(&self) -> &[u8];
    /// Write the record to a writer
    fn write<W: Write>(&self, writer: &mut W) -> Result<usize>;

    /// Return true if the sequence contains only A, C, T and G.
    fn validate_dna(&self) -> bool {
        self.seq().iter().all(|&x| x == b'A' || x == b'C' || x == b'T' || x == b'G')
    }

    /// Return true if the sequence contains only A, C, T, G and N.
    fn validate_dnan(&self) -> bool {
        self.seq().iter().all(|&x| x == b'A' || x == b'C' || x == b'T' || x == b'G' || x == b'N')
    }
}


/// A fastq record that borrows data from an array.
#[derive(Debug)]
pub struct RefRecord<'a> {
    // (start, stop), but might include \r at the end
    head: usize,
    seq: usize,
    sep: usize,
    qual: usize,
    data: &'a [u8],
}

/// A fastq record that ownes its data arrays.
#[derive(Debug)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub sep: Vec<u8>,
    pub qual: Vec<u8>,
}

#[derive(Debug)]
struct IdxRecord {
    head: usize,
    seq: usize,
    sep: usize,
    qual: usize,
    data: (usize, usize),
}


/// Remove a final '\r' from a byte slice
#[inline]
fn trim_winline(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', remaining)) = line.split_last() {
        remaining
    } else {
        line
    }
}


impl<'a> Record for RefRecord<'a> {
    #[inline]
    fn head(&self) -> &[u8] {
        // skip the '@' at the beginning
        trim_winline(&self.data[1 .. self.head])
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        trim_winline(&self.data[self.head + 1 .. self.seq])
    }

    #[inline]
    fn qual(&self) -> &[u8] {
        trim_winline(&self.data[self.sep + 1 .. self.qual])
    }

    #[inline]
    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        writer.write_all(&self.data).map(|_| self.data.len())
    }
}


impl Record for OwnedRecord {
    fn head(&self) -> &[u8] {
        // skip the '@' at the beginning
        &self.head
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn qual(&self) -> &[u8] {
        &self.qual
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        let mut written = 0;
        written += writer.write(b"@")?;
        written += writer.write(self.head())?;
        written += writer.write(b"\n")?;
        written += writer.write(self.seq())?;
        written += writer.write(b"\n")?;
        written += writer.write(&self.sep)?;
        written += writer.write(b"\n")?;
        written += writer.write(self.qual())?;
        written += writer.write(b"\n")?;
        Ok(written)
    }
}


enum IdxRecordResult {
    Incomplete,
    EmptyBuffer,
    Record(IdxRecord),
}


#[inline]
fn read_header(buffer: &[u8]) -> Result<Option<usize>> {
    match buffer.first() {
        None => { Ok(None) },
        Some(&b'@') => {
            Ok(memchr(b'\n', buffer))
        },
        Some(_) => {
            return Err(Error::new(ErrorKind::InvalidData,
                                  "Fastq headers must start with '@'"))
        }
    }
}


#[inline]
fn read_seq(buffer: &[u8]) -> Result<Option<usize>> {
    Ok(memchr(b'\n', buffer))
}


#[inline]
fn read_sep(buffer: &[u8]) -> Result<Option<usize>> {
    match buffer.first() {
        None => { return Ok(None) },
        Some(&b'+') => { Ok(memchr(b'\n', buffer)) },
        Some(_) => {
            return Err(Error::new(ErrorKind::InvalidData,
                                  "Sequence and quality not separated by +"));
        }
    }
}


#[inline]
fn read_qual(buffer: &[u8]) -> Result<Option<usize>> {
    Ok(memchr(b'\n', buffer))
}


impl<'a> RefRecord<'a> {
    /// Copy the borrowed data array and return an owned record.
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            seq: self.seq().to_vec(),
            qual: self.qual().to_vec(),
            head: self.head().to_vec(),
            sep: trim_winline(&self.data[self.seq + 1..self.sep]).to_vec()
        }
    }
}

impl IdxRecord {
    #[inline]
    fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecord<'a> {
        let data = &buffer[self.data.0..self.data.1];
        let datalen = data.len();
        debug_assert!(datalen == self.data.1 - self.data.0);

        debug_assert!(self.head < datalen);
        debug_assert!(self.qual < datalen);
        debug_assert!(self.seq < datalen);
        debug_assert!(self.sep < datalen);
        debug_assert!(self.head < self.seq);
        debug_assert!(self.seq < self.sep);
        debug_assert!(self.sep < self.qual);

        RefRecord {
            data: data,
            head: self.head,
            seq: self.seq,
            sep: self.sep,
            qual: self.qual,
        }
    }

    #[inline]
    fn from_buffer(buffer: &[u8]) -> Result<IdxRecordResult> {
        if buffer.len() == 0 {
            return Ok(IdxRecordResult::EmptyBuffer);
        }


        let head_end = match read_header(buffer)? {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some(val) => val
        };
        let pos = head_end + 1;

        let buffer_ = &buffer[pos..];
        let seq_end = match read_seq(buffer_)? {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some(end) => end + pos
        };
        let pos = seq_end + 1;

        let buffer_ = &buffer[pos..];
        let sep_end = match read_sep(buffer_)? {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some(end) => end + pos,
        };
        let pos = sep_end + 1;

        let buffer_ = &buffer[pos..];
        let qual_end = match read_qual(buffer_)? {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some(end) => end + pos,
        };

        if qual_end - sep_end != seq_end - head_end {
            return Err(Error::new(ErrorKind::InvalidData,
                                  "Sequence and quality length mismatch"));
        }

        Ok(IdxRecordResult::Record(
            IdxRecord {
                data: (0, qual_end + 1),
                head: head_end,
                seq: seq_end,
                sep: sep_end,
                qual: qual_end,
            }
        ))
    }
}


/// Parser for fastq files.
pub struct Parser<R: Read> {
    reader: R,
    buffer: buffer::Buffer,
}


impl<R: Read> Parser<R> {
    /// Create a new fastq parser.
    pub fn new(reader: R) -> Parser<R> {
        Parser {
            reader: reader,
            buffer: buffer::Buffer::new(BUFSIZE),
        }
    }

    /// Apply a function to each fastq record in an file.
    #[inline]
    pub fn each<F>(&mut self, mut func: F) -> Result<()> where F: FnMut(RefRecord) {
        loop {
            match IdxRecord::from_buffer(self.buffer.data())? {
                IdxRecordResult::EmptyBuffer => {
                    self.buffer.clean();
                    if self.buffer.read_into(&mut self.reader)? == 0 {
                        return Ok(())
                    };
                }
                IdxRecordResult::Incomplete => {
                    let freed = self.buffer.clean();
                    if freed == 0 {
                        return Err(Error::new(ErrorKind::InvalidData,
                                              "Fastq record is too long"));
                    }
                    if self.buffer.read_into(&mut self.reader)? == 0 {
                        return Err(Error::new(ErrorKind::InvalidData,
                                              "Possibly truncated input file."))

                    }
                }
                IdxRecordResult::Record(record) => {
                    func(record.to_ref_record(self.buffer.data()));
                    self.buffer.consume(record.data.1 - record.data.0);
                }
            }
        }
    }
}


/// A collection of fastq records used to iterate over records in chunks.
#[derive(Debug)]
pub struct RecordSet {
    buffer: Box<[u8]>,
    records: Vec<IdxRecord>,
}


impl RecordSet {
    fn from_records(buffer: Box<[u8]>, records: Vec<IdxRecord>) -> RecordSet {
        RecordSet {
            buffer: buffer,
            records: records,
        }
    }

    /// Return an iterator over all fastq records in this record set.
    #[inline]
    pub fn iter<'a>(&'a self) -> impl Iterator<Item=RefRecord<'a>> + 'a {
        self.records.iter().map(move |r| r.to_ref_record(&self.buffer))
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}


pub struct RecordSetIter<R: Read> {
    parser: Parser<R>,
    num_records_guess: usize,
    reader_at_end: bool,
}


impl<R: Read> Iterator for RecordSetIter<R> {
    type Item = Result<RecordSet>;

    fn next(&mut self) -> Option<Result<RecordSet>> {
        if self.reader_at_end {
            return None
        }

        let mut records: Vec<IdxRecord> = Vec::with_capacity(self.num_records_guess);

        loop {
            let parse_result = {
                match IdxRecord::from_buffer(self.parser.buffer.data()) {
                    Ok(val) => { val },
                    Err(e) => { return Some(Err(e)) }
                }
            };
            use IdxRecordResult::*;
            match parse_result {
                EmptyBuffer => {
                    self.num_records_guess = records.len() + 1;

                    let buffer = vec![0u8; BUFSIZE].into_boxed_slice();
                    let buffer = self.parser.buffer.replace_buffer(buffer);
                    match self.parser.buffer.read_into(&mut self.parser.reader) {
                        Err(e) => { return Some(Err(e)) },
                        Ok(0) => { self.reader_at_end = true },
                        _ => { }
                    }
                    return Some(Ok(RecordSet::from_records(buffer, records)))
                }
                Incomplete => {
                    if self.reader_at_end {
                        return Some(Err(Error::new(ErrorKind::InvalidData,
                                                   "Truncated input file.")))
                    }
                    self.num_records_guess = records.len() + 1;

                    let buffer = vec![0u8; BUFSIZE].into_boxed_slice();
                    let buffer = self.parser.buffer.replace_buffer(buffer);
                    match self.parser.buffer.read_into(&mut self.parser.reader) {
                        Err(e) => { return Some(Err(e)) }
                        Ok(0) => { return Some(Err(Error::new(ErrorKind::InvalidData,
                                                              "Fastq record is too long."))) }
                        _ => { }
                    }
                    return Some(Ok(RecordSet::from_records(buffer, records)))
                }
                Record(record) => {
                    let (start, end) = (record.data.0, record.data.1);
                    records.push(record);
                    self.parser.buffer.consume(end - start);
                }
            }
        }
    }
}


impl<R: Read> Parser<R> {
    /// Return the fastq records in chunks.
    pub fn record_sets(self) -> RecordSetIter<R> {
        RecordSetIter {
            parser: self,
            reader_at_end: false,
            num_records_guess: 20,
        }
    }

    pub fn parallel_each<O, I, F>(self, n_threads: usize, func: F) -> Result<I>
        where
            I: FromIterator<O>,
            O: Send + 'static,
            F: Send + Sync + 'static,
            F: Fn(Box<Iterator<Item=RecordSet>>) -> O,
    {
        let mut senders: Vec<SyncSender<_>> = vec![];
        let mut threads: Vec<thread::JoinHandle<_>> = vec![];

        let func = Arc::new(func);

        for i in 0..n_threads {
            let (tx, rx): (SyncSender<RecordSet>, _) = sync_channel(2);
            let func = func.clone();

            let thread = thread::Builder::new()
                .name(format!("worker-{}", i))
                .spawn(move || {
                    func(Box::new(rx.into_iter()))
                })
                .expect("Could not start worker threads");

            senders.push(tx);
            threads.push(thread);
        }

        let mut io_error = None;
        let mut send_error = None;
        for (record_set, sender) in self.record_sets().zip(senders.iter().cycle()) {
            match record_set {
                Ok(records) => {
                    if let Err(e) = sender.send(records) {
                        send_error = Some(e);
                        break;
                    }
                },
                Err(e) => {
                    io_error = Some(e);
                    break;
                }
            }
        }

        if send_error.is_some() {
            panic!("Worker thread died.");
        }

        ::std::mem::drop(senders);
        let results = threads.into_iter()
            .map(|thread| thread.join())
            .collect::<Vec<_>>()
            .into_iter()
            .map(|res| res.expect("Panic in worker thread"))
            .collect();

        match io_error {
            Some(e) => { Err(e) },
            None => { Ok(results) }
        }
    }
}


#[cfg(test)]
mod tests {
    use std::io::{Cursor, Write, Seek, SeekFrom, ErrorKind};
    use super::{Parser, Record};

    #[test]
    fn correct() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hallo\nTCC\n+\nabc\n");
        let mut parser = Parser::new(data);
        let mut i: u64 = 0;
        let ok = parser.each(move |record| {
            if i == 0 {
                assert_eq!(record.head(), b"hi");
                assert_eq!(record.seq(), b"NN");
                assert_eq!(record.qual(), b"++");

                let record = record.to_owned_record();
                assert_eq!(record.head(), b"hi");
                assert_eq!(record.seq(), b"NN");
                assert_eq!(record.qual(), b"++");
                let mut out = Cursor::new(Vec::new());
                assert_eq!(record.write(&mut out).unwrap(), 12);
                assert_eq!(&out.into_inner()[..], b"@hi\nNN\n+\n++\n");
            } else {
                assert_eq!(record.head(), b"hallo");
                assert_eq!(record.seq(), b"TCC");
                assert_eq!(record.qual(), b"abc");

                let record = record.to_owned_record();
                assert_eq!(record.head(), b"hallo");
                assert_eq!(record.seq(), b"TCC");
                assert_eq!(record.qual(), b"abc");
                let mut out = Cursor::new(Vec::new());
                assert_eq!(record.write(&mut out).unwrap(), 17);
                assert_eq!(&out.into_inner()[..], b"@hallo\nTCC\n+\nabc\n");
            }
            assert!(i < 2);
            i += 1;
        });
        ok.unwrap();
    }

    #[test]
    fn missing_lines() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN");
        let mut parser = Parser::new(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");
        });
        match ok {
            Err(e) => { assert!(e.kind() == ErrorKind::InvalidData) },
            Ok(()) => { panic!("should fail") },
        }
    }

    #[test]
    fn truncated() {
        let data = Cursor::new(b"@hi\nNN\n+\n++");
        let mut parser = Parser::new(data);
        let ok = parser.each(|_| { assert!(false) });
        assert!(ok.is_err());
    }

    #[test]
    fn second_idline() {
        let data = Cursor::new(b"@hi\nNN\n+hi\n++\n@hi\nNN\n+hi\n++\n");
        let mut parser = Parser::new(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");

            let mut out = Cursor::new(Vec::new());
            assert_eq!(record.write(&mut out).unwrap(), 14);
            assert_eq!(&out.into_inner()[..], b"@hi\nNN\n+hi\n++\n");
        });
        ok.unwrap();
    }

    #[test]
    fn windows_lineend() {
        let data = Cursor::new(b"@hi\r\nNN\r\n+\r\n++\r\n@hi\r\nNN\r\n+\r\n++\r\n");
        let mut parser = Parser::new(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");
        });
        ok.unwrap();
    }

    #[test]
    fn length_mismatch() {
        let data = Cursor::new(b"@hi\nNN\n+\n+\n");
        let mut parser = Parser::new(data);
        let ok = parser.each(|_| { assert!(false) });
        assert!(ok.is_err());
    }

    #[test]
    fn huge_incomplete() {
        let mut data: Cursor<Vec<u8>> = Cursor::new(vec![]);
        data.write(b"@").unwrap();
        for _ in 0..super::BUFSIZE {
            data.write(b"longid").unwrap();
        };
        data.seek(SeekFrom::Start(0)).unwrap();
        let mut parser = Parser::new(data);
        assert!(parser.each(|_| ()).is_err());
    }


    #[test]
    fn refset() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN\n+\n++\n");
        let parser = Parser::new(data);
        let mut count: usize = 0;
        for set in parser.record_sets() {
            let set = set.unwrap();
            for record in set.iter() {
                count += 1;
                assert_eq!(record.head(), b"hi");
                assert_eq!(record.seq(), b"NN");
                assert_eq!(record.qual(), b"++");
            }
        }
        assert_eq!(count, 2)
    }

    #[test]
    fn refset_incomplete() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN\n+\n++");
        let parser = Parser::new(data);
        assert!(parser.record_sets().any(|x| x.is_err()));
    }

    #[test]
    fn refset_huge_incomplete() {
        let mut data: Cursor<Vec<u8>> = Cursor::new(vec![]);
        data.write(b"@").unwrap();
        for _ in 0..super::BUFSIZE {
            data.write(b"longid").unwrap();
        };
        data.seek(SeekFrom::Start(0)).unwrap();
        let parser = Parser::new(data);
        assert!(parser.record_sets().any(|x| x.is_err()));
    }
}
