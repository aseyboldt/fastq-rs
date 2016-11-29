use std::io::{Read, Result, ErrorKind};


pub struct Buffer {
    data: Box<[u8]>,
    start: usize,
    end: usize
}


impl Buffer {
    pub fn new(size: usize) -> Buffer {
        Buffer {
            data: vec![0u8; size].into_boxed_slice(),
            start: 0,
            end: 0
        }
    }

    pub fn len(&self) -> usize {
        self.end.checked_sub(self.start).unwrap()
    }

    pub fn n_free(&self) -> usize {
        self.data.len().checked_sub(self.end).unwrap()
    }

    pub fn replace_buffer(&mut self, mut buffer: Box<[u8]>) -> Box<[u8]> {
        let n_in_buffer = self.len();
        let new_end = (n_in_buffer + 15) & !0x0f;  // make sure next read is aligned
        let num_move = self.end.checked_sub(new_end).unwrap();
        let new_start = self.start.checked_sub(num_move).unwrap();

        assert!(buffer.len() >= new_end);

        {
            let dest = &mut buffer[new_start..new_end];
            let src = &self.data[self.start..self.end];
            dest.copy_from_slice(src);
        }

        ::std::mem::swap(&mut self.data, &mut buffer);
        self.start = new_start;
        self.end = new_end;

        buffer
    }

    /// Move data to the start of the buffer, freeing space at the end.
    ///
    /// Return the amount of *aditional* free space.
    pub fn clean(&mut self) -> usize {
        if self.start == 0 {
            return 0
        }

        let n_in_buffer = self.len();
        let new_end = (n_in_buffer + 15) & !0x0f;  // make sure next read is aligned
        let num_move = self.end.checked_sub(new_end).unwrap();
        let new_start = self.start.checked_sub(num_move).unwrap();

        let dest = self.data[new_start..].as_mut_ptr();
        let src = self.data[self.start..].as_ptr();

        unsafe { ::std::ptr::copy(src, dest, n_in_buffer); }
        self.start = new_start;
        self.end = new_end;

        num_move
    }

    pub fn read_into<R: Read>(&mut self, reader: &mut R) -> Result<usize> {
        let n_free = self.n_free();
        let num_read = if n_free < 4096 { n_free } else { n_free - n_free % 4096 };

        let dest = &mut self.data[self.end..self.end + num_read];

        let n_read;
        loop {
            match reader.read(dest) {
                Err(e) => {
                    if e.kind() != ErrorKind::Interrupted {
                        return Err(e)
                    }
                },
                Ok(val) => { n_read = val; break }
            }
        };
        self.end += n_read;
        Ok(n_read)
    }

    #[inline]
    pub fn data(&self) -> &[u8] {
        &self.data[self.start..self.end]
    }

    #[inline]
    pub fn consume(&mut self, count: usize) {
        self.start = self.start.checked_add(count).unwrap();
        debug_assert!(self.start <= self.end);
    }
}
