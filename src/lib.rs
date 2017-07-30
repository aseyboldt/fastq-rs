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
//! - `Parser::parallel_each`. This is a convenience function that wraps
//!   `Parser::record_sets` but passes the record sets to a number of
//!   background threads. A closure is executed on each thread with an iterator
//!   over record sets. Results from the threads are collected and returned
//!   to the caller.
//!
//! Since fastq file are usually compressed, this crate also includes a function
//! `thread_reader` to offload the decompression to a different core, and
//! `from_path` to guess the compression.
//!
//! # The FastQ standard
//!
//! This library supports Windows and Unix line endings, it does not support
//! the old MAC line ending `\r`. It allows arbitrary data on the separator
//! line between sequence and quality as long as it starts with a `+` (some
//! fastq files repeat the id on this line). It does not validate that the
//! sequence or the quality contain only allowed characters. Sequence and
//! quality must have the same length. They are not allowed to contain
//! newline characters.
//!
//! At the moment it does not make any effort to pair reads. This means that
//! pairs that belong together might end up on different cores in a
//! multithreaded setup. (TODO This should change it the future!).
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
//!     total += 1;
//!     // return `true` if you want to continue iterating
//!     true
//! }).expect("Invalid fastq file");
//! println!("{}", total);
//! ```
//!
//! If we want to do more than just count the number of records (in this
//! example, count how many sequences align to an illumina adapter with a score
//! better than 10), we probably want to use more cores:
//!
//! ```rust,no_run
//! use fastq::{parse_path, Record};
//! use std::env::args;
//! use parasailors as align;
//!
//! extern crate fastq;
//! extern crate parasailors;
//!
//! fn main() {
//!     let filename = args().nth(1);
//!     // Treat "-" as stdin
//!     let path = match filename.as_ref().map(String::as_ref) {
//!         None | Some("-") => { None },
//!         Some(name) => Some(name)
//!     };
//!
//!     parse_path(path, |parser| {
//!         let nthreads = 4;
//!         let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {
//!             // we can initialize thread local variables here.
//!             let adapter = b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
//!             let matrix = align::Matrix::new(align::MatrixType::Identity);
//!             let profile = align::Profile::new(adapter, &matrix);
//!             let mut thread_total = 0;
//!
//!             for record_set in record_sets {
//!                 for record in record_set.iter() {
//!                     let score = align::local_alignment_score(&profile, record.seq(), 5, 1);
//!                     if score > 8 {
//!                         thread_total += 1;
//!                     }
//!                 }
//!             }
//!
//!             // The values we return (it can be any type implementing `Send`)
//!             // are collected from the different threads by
//!             // `parser.parallel_each` and returned. See doc for a description of
//!             // error handling.
//!             thread_total
//!         }).expect("Invalid fastq file");
//!         println!("{}", results.iter().sum::<usize>());
//!     }).expect("Invalid compression");
//! }
//! ```
//!
//! On my feeble 2 core laptop this ends up being bound by the alignment at ~300MB/s,
//! but it should scale well to a larger number of cores.

use std::io::{Write, Result, Read, Error, ErrorKind, Cursor};
use std::thread;
use std::sync::mpsc::{SyncSender, sync_channel};
use std::sync::Arc;
use std::iter::FromIterator;
use std::path::Path;
use memchr::memchr;
use lz4::Decoder;
use flate2::read::GzDecoder;

extern crate memchr;
extern crate lz4;
extern crate flate2;

mod thread_reader;
mod buffer;

pub use thread_reader::thread_reader;

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
    ///
    /// FIXME This might be much faster with a [bool; 256] array
    /// or using some simd instructions (eg with the jetscii crate).
    fn validate_dna(&self) -> bool {
        self.seq().iter().all(|&x| x == b'A' || x == b'C' || x == b'T' || x == b'G')
    }

    /// Return true if the sequence contains only A, C, T, G and N.
    ///
    /// FIXME This might be much faster with a [bool; 256] array
    /// or using some simd instructions (eg with the jetscii crate).
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
    pub sep: Option<Vec<u8>>,
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
        writer.write_all(&self.data)?;
        Ok(self.data.len())
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
        match self.sep {
            Some(ref s) => { written += writer.write(s)? }
            None => { written += writer.write(b"+")? }
        }
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


impl<'a> RefRecord<'a> {
    /// Copy the borrowed data array and return an owned record.
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            seq: self.seq().to_vec(),
            qual: self.qual().to_vec(),
            head: self.head().to_vec(),
            sep: Some(trim_winline(&self.data[self.seq + 1..self.sep]).to_vec())
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
        let seq_end = match memchr(b'\n', buffer_) {
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
        let qual_end = match memchr(b'\n', buffer_) {
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


/// Create a parser and guess the compression.
///
/// At the moment this supports gzip, lz4 and plain fastq files.
/// If the path is None, read from stdin.
///
/// # Examples
///
/// ```rust,no_run
/// use fastq::{parse_path, Record};
/// use std::env::args;
///
/// let filename = args().nth(1);
/// // Accept "-" as stdin
/// let path = match filename.as_ref().map(String::as_ref) {
///     None | Some("-") => { None },
///     Some(name) => Some(name)
/// };
///
/// parse_path(path, |mut parser| {
///     let stopped = parser.each(|record| {
///         // stop parsing if we find a sequnce containing 'N'
///         record.validate_dnan()
///     }).expect("Invalid fastq file");
///     if stopped {
///         println!("The file contains only sequences with ACTGN");
///     } else {
///         println!("The file contains invalid sequences");
///     }
/// }).expect("Invalid compression");
/// ```
pub fn parse_path<P, F, O>(path: Option<P>, func: F) -> Result<O>
        where P: AsRef<Path>,
              F: FnOnce(Parser<Box<&mut Read>>) -> O
{
    let mut reader: Box<Read + Send> = match path {
        None => {
            Box::new(std::io::stdin())
        },
        Some(path) => {
            Box::new(std::fs::File::open(path)?)
        }
    };
    let mut magic_bytes = [0u8; 4];
    reader.read_exact(&mut magic_bytes)?;
    let mut reader = Cursor::new(magic_bytes.to_vec()).chain(reader);
    if unsafe { std::mem::transmute::<_, u32>(magic_bytes.clone()) }.to_le() ==  0x184D2204 {
        let bufsize = 1<<22;
        let queuelen = 2;
        return Ok(thread_reader(bufsize, queuelen, Decoder::new(reader)?, |reader| {
            func(Parser::new(Box::new(reader)))
        }).expect("lz4 reader thread paniced"))
    } else if &magic_bytes[..2] == b"\x1f\x8b" {
        let bufsize = 1<<22;
        let queuelen = 2;
        let reader = GzDecoder::new(reader)?;
        return Ok(thread_reader(bufsize, queuelen, reader, |reader| {
            func(Parser::new(Box::new(reader)))
        }).expect("gzip reader thread paniced"))
    } else if magic_bytes[0] == b'@' {
        Ok(func(Parser::new(Box::new(&mut reader))))
    } else {
        return Err(Error::new(ErrorKind::InvalidData, "Not a gzip, lz4 or plain fastq file"))
    }
}


impl<'a, R: 'a + Read> Parser<R> {
    /// Create a new fastq parser.
    pub fn new(reader: R) -> Parser<R> {
        Parser {
            reader: reader,
            buffer: buffer::Buffer::new(BUFSIZE),
        }
    }

    fn ref_iter(self) -> RecordRefIter<R> {
        RecordRefIter {
            parser: self,
            current_length: None,
            current: None,
        }
    }

    /// Apply a function to each fastq record in an file.
    ///
    /// Stop the parser if the closure returns `false`.
    /// Return `true`, if the parser reached the end of the file.
    #[inline]
    pub fn each<F>(self, mut func: F) -> Result<bool> where F: FnMut(RefRecord) -> bool {
        let mut iter = self.ref_iter();
        loop {
            iter.advance()?;
            match iter.get() {
                None => { return Ok(true) },
                Some(record) => {
                    let go_on = func(record);
                    if !go_on {
                        return Ok(false)
                    }
                }
            }
        }
    }
}


pub struct RecordRefIter<R: Read> {
    parser: Parser<R>,
    current: Option<IdxRecord>,
    current_length: Option<usize>,
}


impl<R: Read> RecordRefIter<R> {
    pub fn get(&self) -> Option<RefRecord> {
        match self.current {
            None => None,
            Some(ref rec) => Some(rec.to_ref_record(self.parser.buffer.data()))
        }
    }

    pub fn advance(&mut self) -> Result<()> {
        let mut buffer = &mut self.parser.buffer;
        let mut reader = &mut self.parser.reader;
        if let Some(len) = self.current_length.take() {
            buffer.consume(len);
        }
        loop {
            match IdxRecord::from_buffer(buffer.data()) {
                Err(e) => { return Err(e) },
                Ok(IdxRecordResult::EmptyBuffer) => {
                    buffer.clean();
                    match buffer.read_into(reader) {
                        Err(e) => { return Err(e) },
                        Ok(0) => {
                            self.current = None;
                            self.current_length = None;
                            return Ok(())
                        },
                        _ => { continue }
                    }
                },
                Ok(IdxRecordResult::Incomplete) => {
                    buffer.clean();
                    if buffer.n_free() == 0 {
                        return Err(Error::new(ErrorKind::InvalidData,
                                              "Fastq record is too long"));
                    }
                    match buffer.read_into(reader) {
                        Err(e) => { return Err(e) }
                        Ok(0) => {
                            return Err(Error::new(ErrorKind::InvalidData,
                                                  "Possibly truncated input file"));
                        },
                        _ => { continue }
                    }
                },
                Ok(IdxRecordResult::Record(record)) => {
                    let length = record.data.1.checked_sub(record.data.0).unwrap();
                    self.current = Some(record);
                    self.current_length = Some(length);
                    return Ok(());
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
    pub fn iter<'a>(&'a self) -> RecordSetItems<'a> {
        RecordSetItems { idx_records: self.records.iter(), buffer: &self.buffer }
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}


pub struct RecordSetItems<'a> {
    idx_records: ::std::slice::Iter<'a, IdxRecord>,
    buffer: &'a [u8],
}


impl<'a> Iterator for RecordSetItems<'a> {
    type Item = RefRecord<'a>;

    #[inline]
    fn next(&mut self) -> Option<RefRecord<'a>> {
        match self.idx_records.next() {
            Some(idx_record) => {
                Some(idx_record.to_ref_record(self.buffer))
            },
            None => None
        }
    }
}


struct RecordSetIter<R: Read> {
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
            let buffer_pos = self.parser.buffer.pos();
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
                    self.num_records_guess = records.len() + 1;

                    let buffer = vec![0u8; BUFSIZE].into_boxed_slice();
                    let buffer = self.parser.buffer.replace_buffer(buffer);
                    if self.parser.buffer.n_free() == 0 {
                        return Some(Err(Error::new(ErrorKind::InvalidData,
                                                   "Fastq record is too long.")))
                    }
                    match self.parser.buffer.read_into(&mut self.parser.reader) {
                        Err(e) => { return Some(Err(e)) }
                        Ok(0) => { return Some(Err(Error::new(ErrorKind::InvalidData,
                                                              "Truncated input file."))) }
                        _ => { }
                    }
                    return Some(Ok(RecordSet::from_records(buffer, records)))
                }
                Record(mut record) => {
                    record.data.0 += buffer_pos;
                    record.data.1 += buffer_pos;
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
    fn record_sets(self) -> RecordSetIter<R> {
        RecordSetIter {
            parser: self,
            reader_at_end: false,
            num_records_guess: 100,
        }
    }

    /// Apply a function to each record, but distribute the work on a number of threads.
    ///
    /// # Parameters
    ///
    /// - `n_threads`: The number of worker threads to start.
    /// - `func`: A closure that is executed on each new thread and takes an
    ///    Iterator of RecordSets as argument.
    ///
    /// This function parses the fastq file and passes record sets to the worker threads.
    ///
    /// It terminates if one of the following happes:
    ///
    /// - The parser exhauts the fastq file. This function waits for all worker
    ///   threads to terminate (their iterator will not yield new values after
    ///   they finish buffered ones). It collects their return values and returns them.
    /// - The parser finds a syntax error. The iterators in the worker threads stop
    ///   yielding new record sets. (The worker threads are not notified of the error).
    ///   The function waits for them to terminate, discards their return values and
    ///   returns an IO error with error kind `std::error::ErrorKind::InvalidData`.
    /// - The underlying Reader yields an io error other an `Interrupted`. The behaviour
    ///   is the same as the previous, but it returns the error value of the inner reader.
    /// - A worker thread terminates before its iterator is exhausted. The parser stops,
    ///   waits for all workers to exit and returns the collected return values. If the
    ///   caller wants to know whether we parsed the whole file, this information must
    ///   be encoded in the return values of the worker threads.
    /// - A worker panics. The parser stops and this function panics.
    ///
    /// # Panics
    ///
    /// Panics, if one of the worker threads panics.
    ///
    /// # Examples
    ///
    /// Parse a fastq file and print a sequence starting with ATTAATTA if the file
    /// contains one (it might not be the first one):
    ///
    /// ```rust
    /// use std::io::{Result, ErrorKind, Cursor};
    /// use fastq::{Parser, Record};
    ///
    /// let reader = Cursor::new(b"@hi\nATTAATTAATTA\n+\n++++++++++++\n");
    /// let parser = Parser::new(reader);
    /// let result: Result<Vec<_>> = parser.parallel_each(4, |record_sets| {
    ///     for record_set in record_sets {
    ///         for record in record_set.iter() {
    ///             if record.seq().starts_with(b"ATTAATTA") {
    ///                 // Early return stops the parser
    ///                 return Some(record.seq().to_vec());
    ///             }
    ///         }
    ///     };
    ///     None
    /// });
    /// match result {
    ///     Ok(res) => {
    ///         match res.iter().filter(|x| x.is_some()).next() {
    ///             None => { assert!(false) } // nothing found
    ///             Some(seq) => {
    ///                 // Yay! we found it.
    ///                 assert_eq!(seq.as_ref().unwrap(), b"ATTAATTAATTA")
    ///             }
    ///         }
    ///     },
    ///     Err(e) => {
    ///         if e.kind() == ErrorKind::InvalidData {
    ///             assert!(false);  // this is not a valid fastq file.
    ///         } else {
    ///             assert!(false);  // some other io error.
    ///         }
    ///     }
    /// }
    pub fn parallel_each<O, S, F>(self, n_threads: usize, func: F) -> Result<S>
        where
            S: FromIterator<O>,
            O: Send + 'static,
            F: Send + Sync + 'static,
            F: Fn(Box<Iterator<Item=RecordSet>>) -> O,
    {
        let mut senders: Vec<SyncSender<_>> = vec![];
        let mut threads: Vec<thread::JoinHandle<_>> = vec![];

        let func = Arc::new(func);

        for i in 0..n_threads {
            let (tx, rx): (SyncSender<RecordSet>, _) = sync_channel(10);
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
        for (record_set, sender) in self.record_sets().zip(senders.iter().cycle()) {
            match record_set {
                Ok(records) => {
                    // We ignore send errors. If the caller wants to know if we exhausted
                    // the file they have to return that information in the worker thread.
                    if let Err(_) = sender.send(records) {
                        break;
                    }
                },
                Err(e) => {
                    io_error = Some(e);
                    break;
                }
            }
        }
        // Make iterators in the workers should stop yielding values
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


/// Step through two fastq files and call a callback for pairs of Records.
///
/// The callback returns a tuple of bools that indicate for each parser
/// if it should advance the iterator. This can be used to deal with
/// single reads that are interleaved in paired sequences and would otherwise
/// lead to misaligned pairs.
pub fn each_zipped<R1, R2, F>(parser1: Parser<R1>, parser2: Parser<R2>, mut callback: F)
    -> Result<(bool, bool)>
    where
        R1: Read,
        R2: Read,
        F: FnMut(Option<RefRecord>, Option<RefRecord>) -> (bool, bool)
{
    let mut iter1 = parser1.ref_iter();
    let mut iter2 = parser2.ref_iter();
    let mut finished = (false, false);
    iter1.advance()?;
    iter2.advance()?;
    loop {
        let advance_flags =  {
            let val1 = if finished.0 { None } else { iter1.get() };
            let val2 = if finished.1 { None } else { iter2.get() };
            finished = (val1.is_none(), val2.is_none());
            callback(val1, val2)
        };
        if advance_flags == (false, false) || finished == (true, true) {
            return Ok(finished);
        }
        if advance_flags.0 && !finished.0 {
            iter1.advance()?;
        }
        if advance_flags.1 && !finished.1 {
            iter2.advance()?;
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
            true
        });
        ok.unwrap();
    }

    #[test]
    fn empty_id() {
        let data = Cursor::new(b"@\nNN\n+\n++\n");
        let mut parser = Parser::new(data);
        parser.each(|record| {
            assert_eq!(record.head(), b"");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");
            true
        }).unwrap();
    }

    #[test]
    fn missing_lines() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN");
        let mut parser = Parser::new(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual(), b"++");
            true
        });
        match ok {
            Err(e) => { assert!(e.kind() == ErrorKind::InvalidData) },
            Ok(_) => { panic!("should fail") },
        }
    }

    #[test]
    fn truncated() {
        let data = Cursor::new(b"@hi\nNN\n+\n++");
        let mut parser = Parser::new(data);
        let ok = parser.each(|_| { assert!(false); true });
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
            true
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
            true
        });
        ok.unwrap();
    }

    #[test]
    fn length_mismatch() {
        let data = Cursor::new(b"@hi\nNN\n+\n+\n");
        let mut parser = Parser::new(data);
        let ok = parser.each(|_| { assert!(false); true });
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
        assert!(parser.each(|_| true).is_err());
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
