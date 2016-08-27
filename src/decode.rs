use std::io::{Result, Read};
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread;
use std::cmp::min;

const QUEUELEN: usize = 1;


pub struct ThreadedReader {
    pub thread_handle: thread::JoinHandle<()>,
    receiver: Receiver<Result<(Box<[u8]>, usize)>>,
    last_buffer: Option<Box<[u8]>>,
    last_buffer_len: usize,
    n_read: usize,
    finished: bool,
}


impl Read for ThreadedReader {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        if self.finished {
            return Ok(0);
        }
        if self.last_buffer.is_none() || (self.n_read == self.last_buffer_len) {
            let (new_buf, new_buf_len) = try!(self.receiver.recv().unwrap());
            if new_buf_len == 0 {
                self.finished = true;
                return Ok(0);
            }

            assert!(new_buf_len <= new_buf.len());

            self.last_buffer_len = new_buf_len;
            self.last_buffer = Some(new_buf);
            self.n_read = 0;
        }

        let n_copy = min(buf.len(), self.last_buffer_len - self.n_read);
        match self.last_buffer {
            None => { assert!(false) },
            Some(ref buffer) => {
                buf[..n_copy].copy_from_slice(&buffer[self.n_read..self.n_read + n_copy]);
            }
        };
        self.n_read += n_copy;
        Ok(n_copy)
    }
}


impl ThreadedReader {
    pub fn new<R: Read + Send + 'static>(mut reader: R, buffsize: usize) -> ThreadedReader {
        let (tx, rx) = sync_channel(QUEUELEN);

        let handle = thread::Builder::new().name("reader-thread".into()).spawn(move || {
            loop {
                let mut buf = vec![0; buffsize].into_boxed_slice();
                match reader.read(&mut buf) {
                    Ok(n_read) => {
                        tx.send(Ok((buf, n_read))).unwrap();
                        if n_read == 0 {
                            break;
                        }
                    }
                    Err(e) => {
                        tx.send(Err(e)).unwrap();
                        break;
                    }
                }
            }
        }).unwrap();
        ThreadedReader {
            thread_handle: handle,
            receiver: rx,
            last_buffer: None,
            last_buffer_len: 0,
            n_read: 0,
            finished: false,
        }
    }
}
