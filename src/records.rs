use memchr::memchr;
use std::io::{Error, ErrorKind, Result, Write};

#[derive(Debug)]
pub struct RefRecord<'a> {
    head: usize,
    seq: usize,
    data: &'a [u8],
    kind: RefRecordKind,
}

#[derive(Debug)]
pub enum RefRecordKind {
    Fastq { sep: usize, qual: usize },
    Fasta,
}

#[derive(Debug)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub kind: OwnedRecordKind,
}

#[derive(Debug)]
pub enum OwnedRecordKind {
    Fastq { sep: Option<Vec<u8>>, qual: Vec<u8> },
    Fasta,
}

#[derive(Debug)]
pub struct IdxRecord {
    head: usize,
    seq: usize,
    pub data: (usize, usize),
    kind: IdxRecordKind,
}

#[derive(Debug)]
pub enum IdxRecordKind {
    Fastq { sep: usize, qual: usize },
    Fasta,
}

// TODO: implementing everything on the enum types above

/// Trait to be implemented by types that represent fastq records.
pub trait Record {
    /// Return the fastq sequence as byte slice
    fn seq(&self) -> &[u8];

    /// Return id line of the record as a byte slice
    fn head(&self) -> &[u8];

    /// Return quality line of a record, if it exists, as a byte slice
    fn qual(&self) -> Option<&[u8]>;

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
    // skip the '@' at the beginning
    fn head(&self) -> &[u8] {
        trim_winline(&self.data[1..self.head])
    }
    #[inline]
    fn qual(&self) -> Option<&[u8]> {
        match self.kind {
            RefRecordKind::Fastq { sep, qual, .. } => Some(trim_winline(&self.data[sep + 1..qual])),
            RefRecordKind::Fasta { .. } => None,
        }
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

impl Record for OwnedRecord {
    fn head(&self) -> &[u8] {
        &self.head
    }

    fn qual(&self) -> Option<&[u8]> {
        match &self.kind {
            OwnedRecordKind::Fastq { qual, .. } => Some(&qual),
            OwnedRecordKind::Fasta => None,
        }
    }
    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        match &self.kind {
            OwnedRecordKind::Fastq { qual, sep } => {
                let mut written = 0;
                written += writer.write(b"@")?;
                written += writer.write(self.head())?;
                written += writer.write(b"\n")?;
                written += writer.write(self.seq())?;
                written += writer.write(b"\n")?;
                match sep {
                    Some(ref s) => written += writer.write(s)?,
                    None => written += writer.write(b"+")?,
                }
                written += writer.write(b"\n")?;
                written += writer.write(&qual)?;
                written += writer.write(b"\n")?;
                Ok(written)
            }
            OwnedRecordKind::Fasta => {
                let mut written = 0;
                written += writer.write(b">")?;
                written += writer.write(self.head())?;
                written += writer.write(b"\n")?;
                written += writer.write(self.seq())?;
                written += writer.write(b"\n")?;
                Ok(written)
            }
        }
    }
}

impl<'a> RefRecord<'a> {
    /// Copy the borrowed data array and return an owned record.
    pub fn to_owned_record(&self) -> OwnedRecord {
        match self.kind {
            RefRecordKind::Fastq { sep, .. } => OwnedRecord {
                seq: self.seq().to_vec(),
                head: self.head().to_vec(),
                kind: OwnedRecordKind::Fastq {
                    qual: self.qual().unwrap().to_vec(),
                    sep: Some(trim_winline(&self.data[self.seq + 1..sep]).to_vec()),
                },
            },
            RefRecordKind::Fasta => OwnedRecord {
                seq: self.seq().to_vec(),
                head: self.head().to_vec(),
                kind: OwnedRecordKind::Fasta,
            },
        }
    }
}

pub enum IdxRecordResult {
    Incomplete,
    EmptyBuffer,
    Record(IdxRecord),
}

impl IdxRecord {
    #[inline]
    pub fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecord<'a> {
        match self.kind {
            IdxRecordKind::Fastq { qual, sep } => {
                let inner_data = &buffer[self.data.0..self.data.1];
                let datalen = inner_data.len();
                debug_assert!(datalen == self.data.1 - self.data.0);

                debug_assert!(self.head < datalen);
                debug_assert!(qual < datalen);
                debug_assert!(self.seq < datalen);
                debug_assert!(sep < datalen);
                debug_assert!(self.head < self.seq);
                debug_assert!(self.seq < sep);
                debug_assert!(sep < qual);

                RefRecord {
                    data: inner_data,
                    head: self.head,
                    seq: self.seq,
                    kind: RefRecordKind::Fastq {
                        sep: sep,
                        qual: qual,
                    },
                }
            }
            IdxRecordKind::Fasta => {
                let inner_data = &buffer[self.data.0..self.data.1];
                let datalen = inner_data.len();
                debug_assert!(datalen == self.data.1 - self.data.0);

                debug_assert!(self.head < datalen);
                debug_assert!(self.seq < datalen);
                debug_assert!(self.head < self.seq);

                RefRecord {
                    data: inner_data,
                    head: self.head,
                    seq: self.seq,
                    kind: RefRecordKind::Fasta,
                }
            }
        }
    }

    #[inline]
    pub fn from_buffer_fq(buffer: &[u8]) -> Result<IdxRecordResult> {
        if buffer.len() == 0 {
            return Ok(IdxRecordResult::EmptyBuffer);
        }

        let head_end = match FastqReader::read_header(buffer)? {
            None => return Ok(IdxRecordResult::Incomplete),
            Some(val) => val,
        };
        let pos = head_end + 1;

        let buffer_ = &buffer[pos..];
        let seq_end = match memchr(b'\n', buffer_) {
            None => return Ok(IdxRecordResult::Incomplete),
            Some(end) => end + pos,
        };
        let pos = seq_end + 1;

        let buffer_ = &buffer[pos..];
        let sep_end = match FastqReader::read_sep(buffer_)? {
            None => return Ok(IdxRecordResult::Incomplete),
            Some(end) => end + pos,
        };
        let pos = sep_end + 1;

        let buffer_ = &buffer[pos..];
        let qual_end = match memchr(b'\n', buffer_) {
            None => return Ok(IdxRecordResult::Incomplete),
            Some(end) => end + pos,
        };

        if qual_end - sep_end != seq_end - head_end {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "Sequence and quality length mismatch",
            ));
        }

        Ok(IdxRecordResult::Record(IdxRecord {
            data: (0, qual_end + 1),
            head: head_end,
            seq: seq_end,
            kind: IdxRecordKind::Fastq {
                sep: sep_end,
                qual: qual_end,
            },
        }))
    }

    #[inline]
    pub fn from_buffer_fa(buffer: &[u8]) -> Result<IdxRecordResult> {
        if buffer.len() == 0 {
            return Ok(IdxRecordResult::EmptyBuffer);
        }

        let head_end = match FastaReader::read_header(buffer)? {
            None => return Ok(IdxRecordResult::Incomplete),
            Some(val) => val,
        };
        let pos = head_end + 1;

        // TODO: handle multiline seqs
        let buffer_ = &buffer[pos..];
        let seq_end = match memchr(b'\n', buffer_) {
            None => return Ok(IdxRecordResult::Incomplete),
            Some(end) => end + pos,
        };
        Ok(IdxRecordResult::Record(IdxRecord {
            data: (0, seq_end + 1),
            head: head_end,
            seq: seq_end,
            kind: IdxRecordKind::Fasta,
        }))
    }
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
