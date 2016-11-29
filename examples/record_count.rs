use fastq::{Parser, thread_reader};
use std::io::{stdin, Read};
use std::fs::File;
use std::env::args;

extern crate fastq;
extern crate lz4;
extern crate parasailors;

const BUFSIZE: usize = 1 << 22;


fn main() {
    let file: Box<Read + Send> = match args().nth(1).as_ref().map(String::as_ref) {
        None | Some("-") => { Box::new(stdin()) },
        Some(path) => { Box::new(File::open(path).unwrap()) }
    };

    let file = lz4::Decoder::new(file).unwrap();
    let total = thread_reader(BUFSIZE, 3, file, |reader| {
        let mut parser = Parser::new(reader);

        let mut total: u64 = 0;
        parser.each(|_| {
            total += 1
        }).expect("Invalid fastq file");
        total
    }).expect("Reader thread paniced");
    println!("{}", total);
}
