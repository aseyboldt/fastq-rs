# A fast parser for fastq.

This library can process fastq files at about the speed of the
coreutils `wc -l` (about 2GB/s on my laptop). It also makes it
easy to distribute the processing of fastq records to many
cores.

See the [documentation](https://docs.rs/fastq/0.2.0/fastq/) for details.

# Examples

Count the number of fastq records that contain an `N`

```rust
use fastq::{Parser, Record};
let reader = ::std::io::stdin();
let mut parser = Parser::new(reader);
let mut total: usize = 0;

parser.each(|record| {
   if record.seq().contains(&b'N') {
       total += 1
   }
}).unwrap();
println!("{}", total);
```

And an (unnecessarily) parallel version of this

```rust
const n_threads: usize = 2;

use fastq::{Parser, Record};
let reader = ::std::io::stdin();
let parser = Parser::new(reader);

let results: Vec<u64> = parser.parallel_each(n_threads, |record_sets| {
    let mut thread_total = 0;
    for record_set in record_sets {
        for record in record_set.iter() {
            if record.seq().contains(&b'N') {
                thread_total += 1;
            }
        }
    }
    thread_total
}).expect("Invalid fastq file");

let total: u64 = results.iter().sum();
println!("{}", total);
```
