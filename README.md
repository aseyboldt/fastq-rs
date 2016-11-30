# A fast parser for fastq.

This library can process fastq files at about the speed of the
coreutils `wc -l` (about 2GB/s on my laptop, `seqan` runs at
about 150MB/s). It also makes it easy to distribute the
processing of fastq records to many cores, without losing much
of the performance.

See the [documentation](https://docs.rs/fastq/) for details.

# Benchmarks

We compare this library with the fastq parser in `rust-bio`,
the C++ library `seqan` 2.0.0, with `kseq.h` and with `wc -l`.

We test 4 scenarios:
- A 2GB test file is uncompressed on a ramdisk. The program
  counts the number of records in the file.
- The test file lz4 compressed on disk, with an empty page
  cache. Again, the program should just count the number
  of records.
- The test file is lz4 compressed on disk with empty page
  cache, but the program sends records to a different
  thread. This thread counts the number of records.
- The same as scenario 3, but with gzip compression.

All measurements are taken with a 2GB test file (TODO describe!)
on a Haskwell i7-4510U @ 2GH. Each program is executed three
times (clearing the os page cache where appropriate) and the best
time is used. Libraries without native support for a compression
algorithm get the input via a pipe from `zcat` or `lz4 -d`.
The C and C++ programs are compiled with gcc 6.2.1 with the
fags `-O3 -march=native`. All programs can be found in the
`examples` directory of this repository.


|           |   ramdisk  |   lz4     | lz4 + thread | gzip    | gzip + thread |
| ----------| -----------| --------- | ------------ | ------- | ------------- |
| `wc -l`   | **2.3GB/s**|  1.2GB/s  |  NA          | 300MB/s |  NA           |
| `fastq`   |   1.9GB/s  |**1.9GB/s**| **1.6GB/s**  | 300MB/s |  300MB/s      |
| `rust-bio`|   730MB/s  |     NA    |  150MB/s     |   NA    |    NA         |
| `seqan`   |   150MB/s  |     NA    |    NA        |   NA    |    NA         |
| `kseq.h`  |   980MB/s  |  680MB/s  |    NA        |   NA    |    NA         |

Some notes from checking `perf record`:

- `wc -l` and `fastq` spend most of the time in `memchr()`, but in contrast
  to `wc`, `fastq` has to check that headers begin with `@` and separator
  lines with `+`. This seems to explain most of the difference in scenario 1.
  `lz4 -d` uses a large buffer size (default 4MB), which seems to prevent
  the operating system from running `lz4` and `wc`  concurrently.
  `fastq` avoids this problem with an internal queue.
- `rust-bio` looses some time copying data and validating utf8.
  The large slowdown in the threaded version stems from the fact, that it
  sends each record to the other thread individually. Each send (I use a
  `sync_channel` from the rust stdlib) requires the use of synchronisation
  primitives, and three allocations for header, sequence and quality.
- `seqan` is busy allocating stuff, and uses (I think) a naive
  implementation of `memchr()` to find line breaks.
- gzip decompression runs at about 320MB/s, so there is not much we can do
  on that front.

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
    // we can initialize thread local variables here.
    let mut thread_total = 0;

    // We iterate over sets of records
    for record_set in record_sets {
        for record in record_set.iter() {
            if record.seq().contains(&b'N') {
                thread_total += 1;
            }
        }
    }

    // The values we return (it can be any type implementing `Send`)
    // are collected from the different threads by
    // `parser.parallel_each` and returned. See doc for a description of
    // the error handling.
    thread_total
}).expect("Invalid fastq file");

// Add up the results from the individual worker threads
let total: u64 = results.iter().sum();
println!("{}", total);
```
