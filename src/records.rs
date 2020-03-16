use memchr::memchr;
use std::io::{Error, ErrorKind, Result, Write};

pub enum RecordTypes {
    Fastq,
    Fasta,
}

/// Trait to be implemented by types that represent fastq records.
pub trait Record {
    /// Return the fastq sequence as byte slice
    fn seq(&self) -> &[u8];

    /// Return id line of the record as a byte slice
    fn head(&self) -> &[u8];

    /// Write the record to a writer
    fn write<W: Write>(&self, writer: &mut W) -> Result<usize>;

    /// Return true if the sequence contains only A, C, T and G.
    ///
    /// FIXME This might be much faster with a [bool; 256] array
    /// or using some simd instructions (eg with the jetscii crate).
    fn validate_dna(&self) -> bool {
        self.seq()
            .iter()
            .all(|&x| x == b'A' || x == b'C' || x == b'T' || x == b'G')
    }

    /// Return true if the sequence contains only A, C, T, G and N.
    ///
    /// FIXME This might be much faster with a [bool; 256] array
    /// or using some simd instructions (eg with the jetscii crate).
    fn validate_dnan(&self) -> bool {
        self.seq()
            .iter()
            .all(|&x| x == b'A' || x == b'C' || x == b'T' || x == b'G' || x == b'N')
    }
}

/// A fasta record that borrows data from an array.
#[derive(Debug)]
pub struct RefRecordFa<'a> {
    // (start, stop), but might include \r at the end
    head: usize,
    seq: usize,
    data: &'a [u8],
}

/// A fasta record that ownes its data arrays.
#[derive(Debug)]
pub struct OwnedRecordFa {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
}

#[derive(Debug)]
pub struct IdxRecordFa {
    head: usize,
    seq: usize,
    pub data: (usize, usize),
}

/// A fastq record that borrows data from an array.
#[derive(Debug)]
pub struct RefRecordFq<'a> {
    // (start, stop), but might include \r at the end
    head: usize,
    seq: usize,
    sep: usize,
    qual: usize,
    data: &'a [u8],
}

/// A fastq record that ownes its data arrays.
#[derive(Debug)]
pub struct OwnedRecordFq {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub sep: Option<Vec<u8>>,
    pub qual: Vec<u8>,
}

#[derive(Debug)]
pub struct IdxRecordFq {
    head: usize,
    seq: usize,
    sep: usize,
    qual: usize,
    pub data: (usize, usize),
}

// TODO: Why is this even needed?
/// Remove a final '\r' from a byte slice
#[inline]
fn trim_winline(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', remaining)) = line.split_last() {
        remaining
    } else {
        line
    }
}

impl<'a> RefRecordFq<'a> {
    #[inline]
    fn qual(&self) -> &[u8] {
        trim_winline(&self.data[self.sep + 1..self.qual])
    }
}

impl<'a> Record for RefRecordFq<'a> {
    #[inline]
    fn head(&self) -> &[u8] {
        // skip the '@' at the beginning
        trim_winline(&self.data[1..self.head])
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        trim_winline(&self.data[self.head + 1..self.seq])
    }

    #[inline]
    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        writer.write_all(&self.data)?;
        Ok(self.data.len())
    }
}

impl OwnedRecordFq {
    fn qual(&self) -> &[u8] {
        &self.qual
    }
}

impl Record for OwnedRecordFq {
    fn head(&self) -> &[u8] {
        // skip the '@' at the beginning
        &self.head
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        let mut written = 0;
        written += writer.write(b"@")?;
        written += writer.write(self.head())?;
        written += writer.write(b"\n")?;
        written += writer.write(self.seq())?;
        written += writer.write(b"\n")?;
        match self.sep {
            Some(ref s) => written += writer.write(s)?,
            None => written += writer.write(b"+")?,
        }
        written += writer.write(b"\n")?;
        written += writer.write(self.qual())?;
        written += writer.write(b"\n")?;
        Ok(written)
    }
}

impl<'a> Record for RefRecordFa<'a> {
    #[inline]
    fn head(&self) -> &[u8] {
        // skip the '@' at the beginning
        trim_winline(&self.data[1..self.head])
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        trim_winline(&self.data[self.head + 1..self.seq])
    }

    #[inline]
    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        writer.write_all(&self.data)?;
        Ok(self.data.len())
    }
}

impl Record for OwnedRecordFa {
    fn head(&self) -> &[u8] {
        // skip the '@' at the beginning
        &self.head
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        let mut written = 0;
        written += writer.write(b">")?;
        written += writer.write(self.head())?;
        written += writer.write(b"\n")?;
        written += writer.write(self.seq())?;
        written += writer.write(b"\n")?;
        Ok(written)
    }
}

pub enum IdxRecordFaResult {
    Incomplete,
    EmptyBuffer,
    Record(IdxRecordFa),
}

pub enum IdxRecordFqResult {
    Incomplete,
    EmptyBuffer,
    Record(IdxRecordFq),
}

trait RecordReader {
    fn read_header(buffer: &[u8]) -> Result<Option<usize>>;
}
struct FastqReader {}
struct FastaReader {}

impl RecordReader for FastaReader {
    #[inline]
    fn read_header(buffer: &[u8]) -> Result<Option<usize>> {
        match buffer.first() {
            None => Ok(None),
            Some(&b'>') => Ok(memchr(b'\n', buffer)),
            Some(_) => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    "Fasta headers must start with '>'",
                ))
            }
        }
    }
}

impl RecordReader for FastqReader {
    #[inline]
    fn read_header(buffer: &[u8]) -> Result<Option<usize>> {
        match buffer.first() {
            None => Ok(None),
            Some(&b'@') => Ok(memchr(b'\n', buffer)),
            Some(_) => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    "Fastq headers must start with '@'",
                ))
            }
        }
    }
}

impl FastqReader {
    #[inline]
    fn read_sep(buffer: &[u8]) -> Result<Option<usize>> {
        match buffer.first() {
            None => return Ok(None),
            Some(&b'+') => Ok(memchr(b'\n', buffer)),
            Some(_) => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    "Sequence and quality not separated by +",
                ));
            }
        }
    }
}

impl<'a> RefRecordFq<'a> {
    /// Copy the borrowed data array and return an owned record.
    pub fn to_owned_record(&self) -> OwnedRecordFq {
        OwnedRecordFq {
            seq: self.seq().to_vec(),
            qual: self.qual().to_vec(),
            head: self.head().to_vec(),
            sep: Some(trim_winline(&self.data[self.seq + 1..self.sep]).to_vec()),
        }
    }
}

impl<'a> RefRecordFa<'a> {
    /// Copy the borrowed data array and return an owned record.
    pub fn to_owned_record(&self) -> OwnedRecordFa {
        OwnedRecordFa {
            seq: self.seq().to_vec(),
            head: self.head().to_vec(),
        }
    }
}

impl IdxRecordFq {
    #[inline]
    pub fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecordFq<'a> {
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

        RefRecordFq {
            data: data,
            head: self.head,
            seq: self.seq,
            sep: self.sep,
            qual: self.qual,
        }
    }

    #[inline]
    pub fn from_buffer(buffer: &[u8]) -> Result<IdxRecordFqResult> {
        if buffer.len() == 0 {
            return Ok(IdxRecordFqResult::EmptyBuffer);
        }

        let head_end = match FastqReader::read_header(buffer)? {
            None => return Ok(IdxRecordFqResult::Incomplete),
            Some(val) => val,
        };
        let pos = head_end + 1;

        let buffer_ = &buffer[pos..];
        let seq_end = match memchr(b'\n', buffer_) {
            None => return Ok(IdxRecordFqResult::Incomplete),
            Some(end) => end + pos,
        };
        let pos = seq_end + 1;

        let buffer_ = &buffer[pos..];
        let sep_end = match FastqReader::read_sep(buffer_)? {
            None => return Ok(IdxRecordFqResult::Incomplete),
            Some(end) => end + pos,
        };
        let pos = sep_end + 1;

        let buffer_ = &buffer[pos..];
        let qual_end = match memchr(b'\n', buffer_) {
            None => return Ok(IdxRecordFqResult::Incomplete),
            Some(end) => end + pos,
        };

        if qual_end - sep_end != seq_end - head_end {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "Sequence and quality length mismatch",
            ));
        }

        Ok(IdxRecordFqResult::Record(IdxRecordFq {
            data: (0, qual_end + 1),
            head: head_end,
            seq: seq_end,
            sep: sep_end,
            qual: qual_end,
        }))
    }
}

impl IdxRecordFa {
    #[inline]
    pub fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecordFa<'a> {
        let data = &buffer[self.data.0..self.data.1];
        let datalen = data.len();
        debug_assert!(datalen == self.data.1 - self.data.0);

        debug_assert!(self.head < datalen);
        debug_assert!(self.seq < datalen);
        debug_assert!(self.head < self.seq);

        RefRecordFa {
            data: data,
            head: self.head,
            seq: self.seq,
        }
    }

    #[inline]
    pub fn from_buffer(buffer: &[u8]) -> Result<IdxRecordFaResult> {
        if buffer.len() == 0 {
            return Ok(IdxRecordFaResult::EmptyBuffer);
        }

        let head_end = match FastaReader::read_header(buffer)? {
            None => return Ok(IdxRecordFaResult::Incomplete),
            Some(val) => val,
        };
        let pos = head_end + 1;

        // TODO: handle multiline seqs
        let buffer_ = &buffer[pos..];
        let seq_end = match memchr(b'\n', buffer_) {
            None => return Ok(IdxRecordFaResult::Incomplete),
            Some(end) => end + pos,
        };
        Ok(IdxRecordFaResult::Record(IdxRecordFa {
            data: (0, seq_end + 1),
            head: head_end,
            seq: seq_end,
        }))
    }
}
