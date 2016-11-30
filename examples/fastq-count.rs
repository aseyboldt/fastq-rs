use fastq::Parser;
use std::io::stdin;

extern crate fastq;


fn main() {
    let mut parser = Parser::new(stdin());

    let mut total: u64 = 0;
    parser.each(|_| {
        total += 1
    }).expect("Invalid fastq file");
    println!("{}", total);
}
