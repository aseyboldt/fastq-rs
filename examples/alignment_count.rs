use fastq::{parse_path, Record};
use std::env::args;
use bio::alignment::pairwise::*;

extern crate fastq;
extern crate bio;

const ADAPTER: &'static [u8] = b"AATGATACGGCGACCACCGAGA\
                                 TCTACACTCTTTCCCTACACGA\
                                 CGCTCTTCCGATCT";

fn main() {
    let filename = args().nth(1);
    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => { None },
        Some(name) => Some(name)
    };

    let results = parse_path(path, |parser| {
        let nthreads = 2;
        let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {

            let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
            let mut aligner = Aligner::new(-5, -1, &score);
            
            let mut thread_total = 0;

            for record_set in record_sets {
                for record in record_set.iter() {
                    let score = aligner.semiglobal(ADAPTER, record.seq()).score;
                    if score > 10 {
                        thread_total += 1;
                    }
                }
            }
            thread_total
        }).expect("Invalid fastq file");
        results
    }).expect("Invalid compression");
    println!("{}", results.iter().sum::<usize>());
}
