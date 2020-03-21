#![no_main]

use fastq::Record;
use std::io::Cursor;

#[macro_use]
extern crate libfuzzer_sys;
extern crate criterion;
extern crate fastq;

fuzz_target!(|data: &[u8]| {
    let reader = Cursor::new(data);
    let mut parser = fastq::Parser::new(reader, fastq::ParserKind::Fasta);

    let mut sum: usize = 0;
    let _ = parser.each(|rec| {
        sum += rec.seq().len();
        true
    });
    criterion::black_box(sum);
});
