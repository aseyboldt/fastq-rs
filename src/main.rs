use fastq::Parser;
use std::io::Read;
use fastq::RecordSet;
use fastq::decode::ThreadedReader;
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread::{Builder, JoinHandle};
use parasailors as align;

extern crate fastq;
extern crate libc;
extern crate twoway;
extern crate lz4;
extern crate bio;
extern crate parasailors;

// const BUFFSIZE: usize = 4194304;
const BUFFSIZE: usize = 8 * 1024;


fn counter_thread(rx: Receiver<Option<RecordSet>>, thread_idx: usize) -> JoinHandle<usize> {
    let name = format!("counter-{}", thread_idx);
    Builder::new()
        .name(name)
        .spawn(move || {
            let matrix = align::Matrix::new(align::MatrixType::Identity);
            let profile = align::Profile::new(b"ATTAATCCAT", &matrix);

            let mut total: usize = 0;
            while let Some(record_set) = rx.recv().unwrap() {
                for record in record_set.records() {
                    let score = align::local_alignment_score(&profile, record.seq(), 11, 1);
                    if score >= 8 {
                        total += 1;
                    }
                }
            }
            total
        })
        .unwrap()
}


fn count_seq<R: Read + Send + 'static>(reader: R) -> usize {
    let parser = Parser::from_reader(reader);

    let n_threads = 2;
    let (txs, rxs): (Vec<_>, Vec<_>) = (0..n_threads)
        .map(|_| sync_channel::<Option<RecordSet>>(2))
        .unzip();

    let threads: Vec<_> = rxs.into_iter()
        .enumerate()
        .map(|(i, rx)| counter_thread(rx, i))
        .collect();

    for (set, tx) in parser.record_sets(50).zip(txs.iter().cycle()) {
        let set = set.unwrap();
        tx.send(Some(set)).unwrap();
    }

    for tx in txs {
        tx.send(None).unwrap();
    }

    let vals = threads.into_iter().map(|t| t.join().unwrap());
    vals.sum()
}


fn main() {
    let reader = std::io::stdin();
    let lz4 = lz4::Decoder::new(reader).unwrap();
    let threaded = ThreadedReader::new(lz4, BUFFSIZE);
    let total = count_seq(threaded);
    println!("total = {}", total);
}
