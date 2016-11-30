use std::io::stdin;
use std::sync::mpsc::sync_channel;
use std::thread::Builder;
use bio::io::fastq::{Reader, Record};

extern crate bio;

fn main() {
    let mut reader = Reader::new(stdin());
    let mut record = Record::new();
    let (tx, rx) = sync_channel(100);

    let handle = Builder::new().name("worker".to_string()).spawn(move || {
        let mut count: usize = 0;
        while let Ok(_) = rx.recv() {
            count += 1;
        }
        count
    }).unwrap();

    while let Ok(_) = reader.read(&mut record) {
        if record.is_empty() { break };
        tx.send(record.clone()).unwrap();
    }
    //for record in reader.records() {
    //    tx.send(record).unwrap();
    //    count += 1;
    //}
    ::std::mem::drop(tx);
    println!("{}", handle.join().unwrap());
}
