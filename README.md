[![Build Status](https://travis-ci.org/aseyboldt/fastq-rs.svg?branch=master)](https://travis-ci.org/aseyboldt/fastq-rs)

# A fast parser for fastq.

This library can process fastq files at about the speed of the
coreutils `wc -l` (about 2GB/s on my laptop, `seqan` runs at
about 150MB/s). It also makes it easy to distribute the
processing of fastq records to many cores, without losing much
of the performance.

See the [documentation](https://docs.rs/fastq/) for details and examples.

# Benchmarks

We compare this library with the fastq parser in `rust-bio`,
the C++ library `seqan` 2.2.0, with `kseq.h` and with `wc -l`.

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
| `fastq`   |   1.9GB/s  |**1.9GB/s**| **1.6GB/s**  |**650MB/s**|**620MB/s** |
| `rust-bio`|   730MB/s  |     NA    |  250MB/s     |   NA    |    NA         |
| `seqan`   |   150MB/s  |     NA    |    NA        |   NA    |    NA         |
| `kseq.h`  |   980MB/s  |  680MB/s  |    NA        |   NA    |    NA         |

Some notes from checking `perf record`:

- `wc -l` and `fastq` spend most of the time in `memchr()`, but in contrast
  to `wc`, `fastq` has to check that headers begin with `@` and separator
  lines with `+` and do some more bookeeping.
  `lz4 -d` uses a large buffer size (default 4MB), which seems to prevent
  the operating system from running `lz4` and `wc` concurrently when connected
  by a pipe. `fastq` avoids this problem with an internal queue.
- `rust-bio` looses some time copying data and validating utf8.
  The large slowdown in the threaded version stems from the fact, that it
  sends each record to the other thread individually. Each send (I use a
  `sync_channel` from the rust stdlib) requires the use of synchronisation
  primitives, and three allocations for header, sequence and quality.
  Collecting records in a `Vec` and sending only after a large number of
  them is available speeds this up from 150MB/s to 250MB/s.
- `seqan` is busy allocating stuff, and uses (I think) a naive
  implementation of `memchr()` to find line breaks.
