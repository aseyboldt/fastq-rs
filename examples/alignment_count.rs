use fastq::{parse_path, Record};
use std::env::args;
use parasailors as align;

extern crate fastq;
extern crate parasailors;

fn main() {
    let filename = args().nth(1);
    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => { None },
        Some(name) => Some(name)
    };

    let results = parse_path(path, |parser| {
        let nthreads = 2;
        let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {
            let adapter = b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
            let matrix = align::Matrix::new(align::MatrixType::Identity);
            let profile = align::Profile::new(adapter, &matrix);
            let mut thread_total = 0;

            for record_set in record_sets {
                for record in record_set.iter() {
                    let score = align::local_alignment_score(&profile, record.seq(), 5, 1);
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
