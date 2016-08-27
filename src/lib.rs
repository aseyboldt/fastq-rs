#![feature(conservative_impl_trait)]

use std::io::{Write, Result, Read, Error, ErrorKind};
use std::fs::File;
use std::borrow::{ToOwned, Borrow};
use std::path::Path;
use memchr::memchr;

extern crate lz4;
extern crate memchr;

pub mod decode;

const BUFSIZE: usize = 128 * 1024;

pub struct Record<T> {
    record: T,
    seq: (usize, usize),
    qual: (usize, usize),
    head: (usize, usize),
}


pub type RefRecord<'a> = Record<&'a [u8]>;
pub type OwnedRecord = Record<Vec<u8>>;
type IdxRecord = Record<(usize, usize)>;


impl<T: Borrow<[u8]>> Record<T> {
    pub fn seq(&self) -> &[u8] {
        &self.record.borrow()[self.seq.0 .. self.seq.1]
    }

    pub fn head(&self) -> &[u8] {
        &self.record.borrow()[self.head.0 .. self.head.1]
    }

    pub fn qual_str(&self) -> &[u8] {
        &self.record.borrow()[self.qual.0 .. self.qual.1]
    }

    pub fn qual(&self) -> Vec<u16> {
        self.qual_str().iter().map(|x| *x as u16).collect()
    }

    pub fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        try!(writer.write(self.record.borrow()));
        Ok(())
    }

    pub fn fastq(&self) -> &[u8] {
        self.record.borrow()
    }
}


enum IdxRecordResult {
    Incomplete,
    EmptyBuffer,
    Record(IdxRecord),
}


fn parse_header(buffer: &[u8]) -> Result<Option<(usize, usize, usize)>> {
    if buffer.len() < 2 {
        return Ok(None)
    }

    if buffer[0] != b'@' {
        return Err(Error::new(ErrorKind::InvalidData, "Fastq headers must start with '@'"));
    }

    match memchr(b'\n', buffer) {
        None => { Ok(None) },
        Some(pos) => { Ok(Some((1, pos, pos + 1))) }
    }
}


fn parse_seq(buffer: &[u8]) -> Result<Option<(usize, usize, usize)>> {
    let linestop = match memchr(b'\n', buffer) {
        None => { return Ok(None) },
        Some(i) => { i }
    };

    if buffer.len() < linestop + 3 {
        return Ok(None)
    }
    if &buffer[linestop .. linestop + 3] != b"\n+\n" {
        Err(Error::new(ErrorKind::InvalidData,
                       "Fastq sequence and quality must be separated by a '+'"))
    } else {
        Ok(Some((0, linestop, linestop + 3)))
    }
}


fn parse_qual(buffer: &[u8]) -> Result<Option<(usize, usize, usize)>> {
    let linestop = match memchr(b'\n', buffer) {
        None => { return Ok(None) },
        Some(i) => { i }
    };

    Ok(Some((0, linestop, linestop + 1)))
}


impl<'a> RefRecord<'a> {
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            record: self.record.to_owned(),
            seq: self.seq,
            qual: self.qual,
            head: self.head,
        }
    }
}

impl IdxRecord {
    fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecord<'a> {
        let data = &buffer[self.record.0..self.record.1];
        let datalen = data.len();
        assert!(datalen == self.record.1 - self.record.0);

        /*
        assert!(self.head.0 < datalen);
        assert!(self.head.1 < datalen);
        assert!(self.qual.0 < datalen);
        assert!(self.qual.1 < datalen);
        assert!(self.seq.0 < datalen);
        assert!(self.seq.1 < datalen);
        assert!(self.head.0 < self.head.1);
        assert!(self.qual.0 < self.qual.1);
        assert!(self.seq.0 < self.seq.1);
        assert!(self.seq.1 - self.seq.0 == self.qual.1 - self.qual.0);
        assert!(self.seq.1 - self.seq.0 == 75);
        */

        RefRecord {
            record: data,
            head: self.head,
            seq: self.seq,
            qual: self.qual,
        }
    }

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
            return Err(Error::new(ErrorKind::InvalidData, "Sequence and quality length mismatch"));
        }

        Ok(IdxRecordResult::Record(
            IdxRecord {
                record: (0, stop),
                head: (head0, head1),
                seq: (seq0, seq1),
                qual: (qual0, qual1),
            }
        ))
    }
}


impl OwnedRecord {
    pub fn to_ref_record(&self) -> RefRecord {
        RefRecord {
            record: &self.record,
            seq: self.seq,
            qual: self.qual,
            head: self.head,
        }
    }
}


pub struct Parser<R: Read> {
    reader: R,
    buffer: Box<[u8]>,
    read_start: usize,
    read_end: usize,
}


impl<R: Read> Parser<R> {
    pub fn from_reader(reader: R) -> Parser<R> {
        Parser {
            reader: reader,
            buffer: vec![0; BUFSIZE].into_boxed_slice(),
            read_start: 0,
            read_end: 0,
        }
    }

    pub fn from_lz4<P: AsRef<Path>>(path: P) -> Result<Parser<lz4::Decoder<File>>> {
        let file = try!(File::open(path));
        let file = try!(lz4::Decoder::new(file));
        Ok(Parser::from_reader(file))
    }

    fn copy_to_start(&mut self) {
        if self.read_start == 0 {
            return
        }

        let n_in_buffer = self.read_end - self.read_start;
        let new_end = (n_in_buffer + 15) & !0x0f;  // make sure next read is aligned
        let new_start = new_end - n_in_buffer;

        let dest = self.buffer[new_start..].as_mut_ptr();
        let src = self.buffer[self.read_start..].as_ptr();

        unsafe { std::ptr::copy(src, dest, n_in_buffer); }
        self.read_start = new_start;
        self.read_end = new_end;
    }

    fn fill_buffer(&mut self) -> Result<bool> {
        self.copy_to_start();

        let n_free = self.buffer.len() - self.read_end;
        let aligned_max = self.read_end + n_free - n_free % 16;

        let n_read = try!(self.reader.read(&mut self.buffer[self.read_end..aligned_max]));
        self.read_end += n_read;
        Ok(n_read > 0)
    }

    pub fn each<F>(mut self, mut func: F) -> Result<()> where F: FnMut(RefRecord) {
        let mut incomplete_record = false;

        'outer: loop {
            let has_new = try!(self.fill_buffer());
            if !has_new {
                if incomplete_record {
                    return Err(Error::new(ErrorKind::InvalidData,
                                          "Possibly truncated input file."))
                } else {
                    return Ok(())
                }
            }

            loop {
                let record_length;
                let buffer = &self.buffer[self.read_start..self.read_end];
                let result = try!(IdxRecord::from_buffer(buffer));
                match result {
                    IdxRecordResult::EmptyBuffer => {
                        incomplete_record = false;
                        continue 'outer;
                    }
                    IdxRecordResult::Incomplete => {
                        if self.read_end - self.read_start > self.buffer.len() {
                            return Err(Error::new(ErrorKind::InvalidData,
                                                  "Fastq record is too long"));
                        }
                        incomplete_record = true;
                        continue 'outer;
                    },
                    IdxRecordResult::Record(record) => {
                        record_length = record.record.1 - record.record.0;
                        func(record.to_ref_record(buffer));
                    }
                }
                self.read_start += record_length;
                assert!(self.read_end >= self.read_start);
            }
        }
    }
}


#[allow(dead_code)]
pub struct RecordSet {
    pub buffer: Box<[u8]>,
    records: Vec<IdxRecord>,
}


impl RecordSet {
    #[allow(needless_lifetimes)]
    #[inline]
    pub fn records<'a>(&'a self) -> impl Iterator<Item=RefRecord<'a>> + 'a {
        self.records.iter().map(move |r| r.to_ref_record(&self.buffer))
    }

    fn from_records(buffer: &[u8], records: Vec<IdxRecord>) -> RecordSet {
        RecordSet {
            buffer: Vec::from(buffer).into_boxed_slice(),
            records: records,
        }
    }
}


pub struct RecordSetIter<R: Read> {
    parser: Parser<R>,
    maxlen: usize,
    incomplete_record: bool,
}


impl<R: Read> Iterator for RecordSetIter<R> {
    type Item = Result<RecordSet>;

    fn next(&mut self) -> Option<Result<RecordSet>> {
        let mut records: Vec<IdxRecord> = Vec::with_capacity(self.maxlen);
        let mut totalbytes: usize = 0;
        let is_empty = self.parser.read_end == self.parser.read_start;

        if is_empty || self.incomplete_record {
            match self.parser.fill_buffer() {
                Err(e) => { return Some(Err(e)) },
                Ok(true) => { },  // new data in buffer
                Ok(false) => {  // end of stream
                    if self.incomplete_record {
                        return Some(Err(Error::new(ErrorKind::InvalidData,
                                                   "Incomplete fastq record".to_string())))
                    } else {
                        return None
                    }
                },
            }
        }
        let buffer = &self.parser.buffer;
        let start = self.parser.read_start;

        for _ in 0..self.maxlen {
            let record_length;
            let buffer = &buffer[self.parser.read_start..self.parser.read_end];
            match IdxRecord::from_buffer(buffer) {
                Err(e) => { return Some(Err(e)) },
                Ok(IdxRecordResult::EmptyBuffer) => {
                    self.incomplete_record = false;
                    break
                },
                Ok(IdxRecordResult::Incomplete) => {
                    self.incomplete_record = true;
                    break
                },
                Ok(IdxRecordResult::Record(mut record)) => {
                    self.incomplete_record = false;
                    record_length = record.record.1 - record.record.0;
                    record.record.0 += totalbytes;
                    record.record.1 += totalbytes;
                    totalbytes += record_length;
                    records.push(record);
                }
            }
            self.parser.read_start += record_length;
            assert!(self.parser.read_end >= self.parser.read_start);
        }

        Some(Ok(RecordSet::from_records(&buffer[start..start+totalbytes], records)))
    }
}


impl<R: Read> Parser<R> {
    pub fn record_sets(self, max_size: usize) -> RecordSetIter<R> {
        RecordSetIter {
            parser: self,
            maxlen: max_size,
            incomplete_record: false,
        }
    }
}


#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use super::Parser;

    #[test]
    fn correct() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN\n+\n++\n");
        let parser = Parser::from_reader(data);
        let ok = parser.each(|record| {
            assert_eq!(record.head(), b"hi");
            assert_eq!(record.seq(), b"NN");
            assert_eq!(record.qual_str(), b"++");
        });
        ok.unwrap();
    }

    #[test]
    fn truncated() {
        let data = Cursor::new(b"@hi\nNN\n+\n++");
        let parser = Parser::from_reader(data);
        let ok = parser.each(|_| { assert!(false) });
        assert!(ok.is_err());
    }

    #[test]
    fn length_mismatch() {
        let data = Cursor::new(b"@hi\nNN\n+\n+\n");
        let parser = Parser::from_reader(data);
        let ok = parser.each(|_| { assert!(false) });
        assert!(ok.is_err());
    }

    #[test]
    fn refset() {
        let data = Cursor::new(b"@hi\nNN\n+\n++\n@hi\nNN\n+\n++\n");
        let parser = Parser::from_reader(data);
        let mut count: usize = 0;
        for set in parser.record_sets(1) {
            let set = set.unwrap();
            for record in set.records() {
                count += 1;
                assert_eq!(record.head(), b"hi");
                assert_eq!(record.seq(), b"NN");
                assert_eq!(record.qual_str(), b"++");
            }
        }
        assert_eq!(count, 2)
    }
}
