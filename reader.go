// Package goseq provides a DNA/Protein sequence reader interface for fasta/fastq
// files. It can efficiently read from files with minimal allocations.
//
// Example usage to convert a fastq file to fasta on stdout:
//
//    rdr, err := goseq.Open("test.fastq")
//    if err != nil {
//      log.Fatalf(err)
//    }
//    for rdr.Next() {
//      fmt.Println(">"+ rdr.Identifier())
//      fmt.Println(rdr.Sequence())
//    }
//    if err := rdr.Err(); err != nil && err != io.EOF {
//      log.Fatalf(err)
//    }
//
package goseq

import (
	"bufio"
	"compress/bzip2"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// Reader is a sequence file reader interface for Fasta/Fastq files.
type Reader interface {
	// Next advances the reader to the next sequence. It returns false if no more
	// sequences are available, or an error occurs.
	Next() bool

	// Identifier returns the identifier text for the current sequence record.
	Identifier() string

	// Sequence returns the sequence content for the current sequence record.
	Sequence() string

	// Err returns the last error that occured during reading.
	Err() error

	// Progress returns the percentage progress through the input file (0.0-100.0)
	// NB especially for compressed files this may not update due to buffering.
	// If an error occurs, returns -1.0
	Progress() float64
}

// ByteReader provides a high-performance byte-level reader for each sequence
// in the input, in addition to the standard Reader methods.
type ByteReader interface {
	Reader

	// SequenceBytes calls the callback function for every byte of sequence data
	// in the input, and does zero allocations.
	SequenceBytes(func(byte)) error
}

//////////////////////

type fastFastaReader struct {
	f       *os.File
	r       io.Reader
	br      *bufio.Reader
	lastErr error

	identifierBytes []byte
	lastBytes       []byte
}

func (f *fastFastaReader) Progress() float64 {
	info, err := f.f.Stat()
	if err != nil || info.Size() <= 0 {
		return -1.0
	}
	pos, err := f.f.Seek(0, os.SEEK_CUR)
	if err != nil {
		return -1.0
	}

	return float64(pos*100.0) / float64(info.Size())
}

var errNotStarted = fmt.Errorf("goseq: Next() not yet called")

func Open(filename string) (Reader, error) {
	ff, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	r := io.Reader(ff)
	if strings.HasSuffix(filename, ".gz") {
		r, err = gzip.NewReader(ff)
		if err != nil {
			return nil, err
		}
	}
	if strings.HasSuffix(filename, ".bz2") {
		r = bzip2.NewReader(ff)
	}
	f := &fastFastaReader{f: ff, r: r}
	f.br = bufio.NewReader(r)
	f.lastErr = errNotStarted

	bb, err := f.br.Peek(1)
	if len(bb) != 1 {
		return f, err
	}
	if bb[0] == '@' {
		return fastqOpenFrom(f)
	}
	return f, err
}

func (f *fastFastaReader) Next() bool {
	if f.lastErr == errNotStarted {
		f.identifierBytes, f.lastErr = f.br.ReadBytes('\n')
	}
	if f.lastErr == io.EOF {
		f.f.Close()
		return false
	}
	if f.lastErr != nil {
		return false
	}
	if len(f.identifierBytes) == 0 {
		return false
	}
	return f.identifierBytes[0] == '>'
}

func (f *fastFastaReader) Err() error {
	return f.lastErr
}

func (f *fastFastaReader) Identifier() string {
	// slice off '>' and '\n'
	return string(f.identifierBytes[1 : len(f.identifierBytes)-1])
}

func (f *fastFastaReader) Sequence() string {
	seq := make([]byte, 0, 65535)
	for {
		f.lastBytes, f.lastErr = f.br.ReadSlice('\n')
		if f.lastErr == bufio.ErrBufferFull {
			f.lastErr = nil
		}

		if f.lastErr != nil {
			return string(seq)
		}

		if len(f.lastBytes) == 0 {
			return string(seq)
		}

		if f.lastBytes[0] == '>' {
			f.identifierBytes = f.lastBytes
			return string(seq)
		}

		seq = append(seq, f.lastBytes[:len(f.lastBytes)-1]...)
	}
}

func (f *fastFastaReader) SequenceBytes(eachbyte func(byte)) error {
	for {
		f.lastBytes, f.lastErr = f.br.ReadSlice('\n')
		if f.lastErr == bufio.ErrBufferFull {
			f.lastErr = nil
		}

		if len(f.lastBytes) > 1 {
			if f.lastBytes[0] == '>' {
				f.identifierBytes = f.lastBytes
				return nil
			}

			for _, b := range f.lastBytes[:len(f.lastBytes)-1] {
				eachbyte(b)
			}
		}
		if f.lastErr == io.EOF {
			return nil
		}
		if f.lastErr != nil {
			return f.lastErr
		}
	}
}
