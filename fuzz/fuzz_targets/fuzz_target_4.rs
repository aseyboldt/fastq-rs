#![no_main]

use fastq::Record;
use std::io::Cursor;

#[macro_use]
extern crate libfuzzer_sys;
extern crate fastq;

fuzz_target!(|data: &[u8]| {
    let reader = Cursor::new(data);
    let mut parser = fastq::Parser::new(reader, fastq::ParserKind::Fasta);

    let _: Result<Vec<_>, _> = parser.parallel_each(3, |sets| {
        for set in sets {
            for rec in set.iter() {
                rec.seq();
            }
        }
        true
    });
});
