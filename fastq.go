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

type fastFastqReader struct {
	f  *os.File
	r  io.Reader
	br *bufio.Reader

	identifierBytes []byte
	lastSeqLen      int
	lastBytes       []byte
	lastErr         error
}

func fastqOpenFrom(f *fastFastaReader) (Reader, error) {
	fq := &fastFastqReader{
		f: f.f, r: f.r, br: f.br,
		lastErr: f.lastErr,
	}
	return fq, nil
}

func OpenFastq(filename string) (Reader, error) {
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
	f := &fastFastqReader{f: ff, r: r}
	f.br = bufio.NewReader(r)
	f.lastErr = errNotStarted

	bb, err := f.br.Peek(1)
	if len(bb) != 1 {
		return f, err
	}
	if bb[0] != '@' {
		return nil, fmt.Errorf("goseq: not a fastq file")
	}
	return f, err
}

func (f *fastFastqReader) Progress() float64 {
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

func (f *fastFastqReader) Next() bool {
	var skipped int
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

	if f.lastSeqLen > 0 { // need to skip over the quality scores

		// NB lastSeqLen doesn't include any \n so this discards
		// strictly less than or equal to next sequence position
		for f.lastSeqLen > 0 && f.lastErr == nil {
			skipped, f.lastErr = f.br.Discard(f.lastSeqLen)
			f.lastSeqLen -= skipped
		}

		if f.lastErr != nil {
			return false
		}

		// if quality scores are line-wrapped, then this might be a short read, so
		// read until an '@' starts the line (or an error occurs)
		f.identifierBytes, f.lastErr = f.br.ReadBytes('\n')
		for len(f.identifierBytes) > 0 && f.identifierBytes[0] != '@' && f.lastErr == nil {
			f.identifierBytes, f.lastErr = f.br.ReadBytes('\n')
		}
	}

	if len(f.identifierBytes) == 0 {
		return false
	}
	return f.identifierBytes[0] == '@'
}

func (f *fastFastqReader) Err() error {
	return f.lastErr
}

func (f *fastFastqReader) Identifier() string {
	// slice off '@' and '\n'
	return string(f.identifierBytes[1 : len(f.identifierBytes)-1])
}

func (f *fastFastqReader) Sequence() string {
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

		if f.lastBytes[0] == '+' {
			f.lastSeqLen = len(seq)
			return string(seq)
		}

		seq = append(seq, f.lastBytes[:len(f.lastBytes)-1]...)
	}
}

func (f *fastFastqReader) SequenceBytes(eachbyte func(byte)) error {
	for {
		f.lastBytes, f.lastErr = f.br.ReadSlice('\n')
		if f.lastErr == bufio.ErrBufferFull {
			f.lastErr = nil
		}

		if len(f.lastBytes) > 1 {
			if f.lastBytes[0] == '+' {
				return nil
			}

			for _, b := range f.lastBytes[:len(f.lastBytes)-1] {
				eachbyte(b)
				f.lastSeqLen++
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
