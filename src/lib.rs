#![feature(conservative_impl_trait)]

//! # Examples
//! Count the number of fastq records that contain an `N`
//!
//! ```rust
//! use fastq::{Parser, Record};
//! let reader = ::std::io::Cursor::new(b"@hi\nNN\n+\n++\n");
//! let mut parser = Parser::new(reader);
//! let mut total: usize = 0;
//! parser.each(|record| {
//!     if record.seq().contains(&b'N') {
//!         total += 1
//!     }
//! }).unwrap();
//! assert_eq!(total, 1);
//! ```


use std::io::{Write, Result, Read, Error, ErrorKind};
use std::borrow::{ToOwned, Borrow};
use std::thread;
use std::sync::mpsc::{SyncSender, sync_channel};
use std::sync::Arc;
use std::iter::FromIterator;
use memchr::memchr;

extern crate memchr;

pub mod thread_reader;
mod buffer;

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
    fn write<W: Write>(&self, writer: &mut W) -> Result<()>;
}


/// A fastq record that borrows data from an array.
#[derive(Debug)]
pub struct RefRecord<'a> {
    seq: (usize, usize),
    qual: (usize, usize),
    head: (usize, usize),
    data: &'a [u8],
}

/// A fastq record that ownes its data array.
#[derive(Debug)]
pub struct OwnedRecord {
    seq: (usize, usize),
    qual: (usize, usize),
    head: (usize, usize),
    data: Vec<u8>,
}

#[derive(Debug)]
struct IdxRecord {
    seq: (usize, usize),
    qual: (usize, usize),
    head: (usize, usize),
    data: (usize, usize),
}

impl<'a> Record for RefRecord<'a> {
    fn seq(&self) -> &[u8] {
        &self.data[self.seq.0 .. self.seq.1]
    }

    fn head(&self) -> &[u8] {
        &self.data[self.head.0 .. self.head.1]
    }

    fn qual(&self) -> &[u8] {
        &self.data[self.qual.0 .. self.qual.1]
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(self.data.borrow())
    }
}


impl Record for OwnedRecord {
    fn seq(&self) -> &[u8] {
        &self.data[self.seq.0 .. self.seq.1]
    }

    fn head(&self) -> &[u8] {
        &self.data[self.head.0 .. self.head.1]
    }

    fn qual(&self) -> &[u8] {
        &self.data[self.qual.0 .. self.qual.1]
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(self.data.borrow())
    }
}


enum IdxRecordResult {
    Incomplete,
    EmptyBuffer,
    Record(IdxRecord),
}


fn end_of_line(buffer: &[u8]) -> Option<(usize, usize)> {
    match memchr(b'\n', buffer) {
        None => { None },
        Some(pos) => {
            if unsafe { pos > 0 && *buffer.get_unchecked(pos - 1) == b'\r' } {
                Some((pos - 1, pos + 1))
            } else {
                Some((pos, pos + 1))
            }
        }
    }
}


fn parse_header(buffer: &[u8]) -> Result<Option<(usize, usize, usize)>> {
    match buffer.first() {
        None => { return Ok(None) },
        Some(&b'@') => { },
        Some(_) => {
            return Err(Error::new(ErrorKind::InvalidData,
                                  "Fastq headers must start with '@'"))
        }
    }
    Ok(end_of_line(buffer).map(|(stop, end)| (1, stop, end)))
}


fn parse_seq(buffer: &[u8]) -> Result<Option<(usize, usize, usize)>> {
    let (stop_seq, start_sep) = match end_of_line(buffer) {
        None => { return Ok(None) },
        Some(val) => val
    };

    match buffer[start_sep..].first() {
        None => { return Ok(None) },
        Some(&b'+') => { },
        Some(_) => {
            return Err(Error::new(ErrorKind::InvalidData,
                                  "Sequence and quality not separated by +"));
        }
    }

    let (_, end) = match end_of_line(&buffer[start_sep..]) {
        None => { return Ok(None) },
        Some((a, b)) => (a + start_sep, b + start_sep),
    };
    Ok(Some((0, stop_seq, end)))
}


fn parse_qual(buffer: &[u8]) -> Result<Option<(usize, usize, usize)>> {
    Ok(end_of_line(buffer).map(|(stop, end)| (0, stop, end)))
}


impl<'a> RefRecord<'a> {
    /// Copy the borrowed data array and return an owned record.
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            data: self.data.to_owned(),
            seq: self.seq,
            qual: self.qual,
            head: self.head,
        }
    }
}

impl IdxRecord {
    fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecord<'a> {
        let data = &buffer[self.data.0..self.data.1];
        let datalen = data.len();
        debug_assert!(datalen == self.data.1 - self.data.0);

        debug_assert!(self.head.0 < datalen);
        debug_assert!(self.head.1 < datalen);
        debug_assert!(self.qual.0 < datalen);
        debug_assert!(self.qual.1 < datalen);
        debug_assert!(self.seq.0 < datalen);
        debug_assert!(self.seq.1 < datalen);
        debug_assert!(self.head.0 < self.head.1);
        debug_assert!(self.qual.0 < self.qual.1);
        debug_assert!(self.seq.0 < self.seq.1);
        debug_assert!(self.seq.1 - self.seq.0 == self.qual.1 - self.qual.0);

        RefRecord {
            data: data,
            head: self.head,
            seq: self.seq,
            qual: self.qual,
        }
    }

    #[inline]
    fn from_buffer(buffer: &[u8]) -> Result<IdxRecordResult> {
        if buffer.len() == 0 {
            return Ok(IdxRecordResult::EmptyBuffer);
        }

        let (head0, head1, stop) = match try!(parse_header(buffer)) {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some(val) => val
        };


        let buffer_ = &buffer[stop..];
        let (seq0, seq1, stop) = match try!(parse_seq(buffer_)) {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some((seq0, seq1, stop_)) => (seq0 + stop, seq1 + stop, stop_ + stop)
        };

        let buffer_ = &buffer[stop..];
        let (qual0, qual1, stop) = match try!(parse_qual(buffer_)) {
            None => { return Ok(IdxRecordResult::Incomplete) },
            Some((qual0, qual1, stop_)) => (qual0 + stop, qual1 + stop, stop_ + stop)
        };

        if seq1 - seq0 != qual1 - qual0 {
            return Err(Error::new(ErrorKind::InvalidData,
                                  "Sequence and quality length mismatch"));
        }

        Ok(IdxRecordResult::Record(
            IdxRecord {
                data: (0, stop),
                head: (head0, head1),
                seq: (seq0, seq1),
                qual: (qual0, qual1),
            }
        ))
    }
}


impl OwnedRecord {
    #[inline]
    pub fn to_ref_record(&self) -> RefRecord {
        RefRecord {
            data: &self.data,
            seq: self.seq,
            qual: self.qual,
            head: self.head,
        }
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

    pub fn apply_threaded<O, I, F>(self, n_threads: usize, func: F) -> Result<I>
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
    use std::io::{Cursor, Write, Seek, SeekFrom};
    use super::{Parser, Record};

    #[test]
    fn correct() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN\n+\n++\n");
        let mut parser = Parser::new(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");
        });
        ok.unwrap();
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
        let data = Cursor::new(b"@hi\nNN\n+hi\n++\n@hi\nNN\n+blubb\n++\n");
        let mut parser = Parser::new(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");
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
