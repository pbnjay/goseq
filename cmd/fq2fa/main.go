// Command fq2fa is a tool to split and/or convert fasta/fastq files into multiple
// fasta files. It will open many input files at once (up to the number of CPU
// cores), and distribute sequences evenly to all output files.
//
//    USAGE: fq2fa [options] file1.fastq file2.fastq ...
//
// Options:
//
//    -n int
//          number of output files to split into (default 1)
//    -o filename
//          output filename pattern (use %d for multiple output sequence number)
//    -z    gzip output files
//
//
// Examples:
//
//    # convert a fastq file into fasta (automatically named myseqs.fasta)
//    fq2fa /path/to/myseqs.fastq
//
//    # convert all fastq files in a folder to 64 gzipped fasta files
//    fq2fa -n 64 -z -o newfiles%03d.fasta.gz /path/to/*.fastq
//
//    # convert all fastq files in a folder to a single fasta file
//    fq2fa -o simple.fasta /path/to/*.fastq
//
package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/pbnjay/goseq"
)

type Seq struct {
	Identifier string
	Sequence   string
}

type st struct {
	cpu int
	pct float64
}

var (
	wg = &sync.WaitGroup{}
)

func statusline(ncpu int, reporter chan st) {
	statuses := make([]float64, ncpu)

	var fmtPct func(pct float64) string

	switch {
	case ncpu <= 10:
		fmtPct = func(pct float64) string {
			if pct > 100 {
				return "   -    "
			}
			// 100.00%
			return fmt.Sprintf("%6.2f%% ", pct)
		}
	case ncpu <= 16:
		fmtPct = func(pct float64) string {
			if pct > 100 {
				return "  -  "
			}
			// 100%
			return fmt.Sprintf("%3d%% ", int(pct))
		}

	case ncpu <= 40:
		fmtPct = func(pct float64) string {
			if pct > 100.0 {
				return "--"
			}
			if pct > 99.0 {
				return "99"
			}
			// 00-99
			return fmt.Sprintf("%02d", pct)
		}

	case ncpu <= 100:
		fmtPct = func(pct float64) string {
			if pct > 100 {
				return "-"
			}
			if pct > 90.0 {
				return "z"
			}
			// a-z
			return fmt.Sprintf("%c", 'a'+int(26.0*pct/100.0))
		}
	}

	tic := time.NewTicker(time.Second)
	for {
		select {
		case <-tic.C:
			line := make([]string, ncpu)
			for i, sp := range statuses {
				line[i] = fmtPct(sp)
			}
			fmt.Fprint(os.Stderr, "\r"+strings.Join(line, ""))
			os.Stderr.Sync()
		case s := <-reporter:
			statuses[s.cpu] = s.pct
		}
	}
}

func fanFiles(ncpu int, files []string, seqChan chan Seq, reporter chan st) {
	wg2 := &sync.WaitGroup{}
	fnChan := make(chan string)
	for i := 0; i < ncpu; i++ {
		wg2.Add(1)
		go convertFile(i, wg2, fnChan, seqChan, reporter)
	}
	for _, fn := range files {
		fnChan <- fn
	}
	close(fnChan)
	wg2.Wait()
}

func convertFile(cpu int, wg2 *sync.WaitGroup, fnChan chan string, seqChan chan Seq, reporter chan st) {
	for fn := range fnChan {
		rdr, err := goseq.Open(fn)
		if err != nil {
			panic(err)
		}
		for rdr.Next() {
			reporter <- st{cpu: cpu, pct: rdr.Progress()}
			seqChan <- Seq{Identifier: rdr.Identifier(), Sequence: rdr.Sequence()}
		}
		if err := rdr.Err(); err != nil && err != io.EOF {
			panic(err)
		}
	}
	wg2.Done()
}

func stripExtenstions(fn string) string {
	ext := filepath.Ext(fn)
	for ext != "" {
		fn = strings.TrimSuffix(fn, ext)
		ext = filepath.Ext(fn)
	}
	return fn
}

func seqWriter(outfilename string, seqChan chan Seq) {
	f, err := os.Create(outfilename)
	if err != nil {
		panic(err)
	}
	w := io.WriteCloser(f)
	if strings.HasSuffix(outfilename, "gz") {
		w = gzip.NewWriter(f)
	}
	for s := range seqChan {
		fmt.Fprintln(w, ">"+s.Identifier)
		for len(s.Sequence) > 80 {
			fmt.Fprintln(w, s.Sequence[:80])
			s.Sequence = s.Sequence[80:]
		}
		fmt.Fprintln(w, s.Sequence)
	}
	if zw, ok := w.(*gzip.Writer); ok {
		zw.Flush()
		w.Close()
	}
	f.Close()
	wg.Done()
}

func fanWriters(n int, pattern string, seqChan chan Seq) {
	subs := make([]chan Seq, n)
	for i := 0; i < n; i++ {
		wg.Add(1)
		subs[i] = make(chan Seq, 5)
		go seqWriter(fmt.Sprintf(pattern, i), subs[i])
	}

	idx := 0
	for s := range seqChan {
		subs[idx%n] <- s
		idx++
	}
	for i := 0; i < n; i++ {
		close(subs[i])
	}
	wg.Done()
	log.Printf("Split %d sequences.", idx)
}

func main() {
	outputpattern := flag.String("o", "", "output `filename` pattern (use %d for multiple output sequence number)")
	noutputs := flag.Int("n", 1, "number of output files to split into")
	compress := flag.Bool("z", false, "gzip output files")
	flag.Parse()

	files := flag.Args()
	if len(files) == 0 {
		fmt.Fprintln(os.Stderr, "no input files provided")
		os.Exit(1)
	}

	if *outputpattern == "" {
		if len(files) == 1 {
			*outputpattern = stripExtenstions(files[0])
		} else {
			*outputpattern = filepath.Base(files[0])
		}

		switch {
		case *noutputs == 1:
			// no leading zeroes
		case *noutputs < 10:
			*outputpattern += "%01d"
		case *noutputs < 100:
			*outputpattern += "%02d"
		case *noutputs < 1000:
			*outputpattern += "%03d"
		case *noutputs < 10000:
			*outputpattern += "%04d"
		default:
			*outputpattern += "%08d"
		}
		*outputpattern += ".fasta"
	}

	if *compress && !strings.HasSuffix(*outputpattern, "gz") {
		*outputpattern += ".gz"
	}

	nc := runtime.NumCPU()
	if len(files) < nc {
		nc = len(files)
	}
	fmt.Fprintf(os.Stderr, "Using %d CPUs for splitting\n", nc)
	fmt.Fprintf(os.Stderr, "Splitting %d inputs into %d outputs.\n", len(files), *noutputs)
	fmt.Fprintf(os.Stderr, "Output file pattern is: '%s'\n\n", *outputpattern)

	seqChan := make(chan Seq)
	wg.Add(1)

	if *noutputs == 1 {
		go seqWriter(*outputpattern, seqChan)
	} else {
		go fanWriters(*noutputs, *outputpattern, seqChan)
	}

	reporter := make(chan st)
	go statusline(nc, reporter)

	fanFiles(nc, files, seqChan, reporter)

	close(seqChan)
	wg.Wait()
}
