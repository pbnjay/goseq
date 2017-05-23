package goseq

import (
	"fmt"
	"strings"
)

// Kmer is a packed k-mer DNA representation.
type Kmer uint64

// InvalidKmer is a sentinal value for an invalid Kmer.
//
// WARNING: If k=32, then even InvalidKmer represents a valid 32-mer and
// should not be used.
const InvalidKmer = ^Kmer(0)

// MaxKmerSize is the maximum K (bases) that can be represented in a Kmer type.
const MaxKmerSize = 32

// BaseKmer tracks some basic values to make Kmer manipulation faster while
// still allowing code to use multiple values of K
type BaseKmer struct {
	k    int
	mask Kmer
}

// NewKmerBase returns a manipulator type that can be used to work with a
// specific kmer size of k. It will return an error if the specific k cannot be
// represented using this implementation (e.g. k > MaxKmerSize).
func NewKmerBase(k uint) (*BaseKmer, error) {
	if k > MaxKmerSize {
		return nil, fmt.Errorf("kmer size k=%d is too large (max %d)", k, MaxKmerSize)
	}
	b := &BaseKmer{
		k:    int(k),
		mask: Kmer(1<<(k*2)) - 1,
	}
	return b, nil
}

// Length returns the base length of the kmers.
func (b *BaseKmer) Length() int {
	return b.k
}

// Count returns the total possible kmers than can be represented.
func (b *BaseKmer) Count() uint64 {
	return 1 + uint64(b.mask)
}

// AppendBase shifts the current kmer left by one base (dropping the left-most),
// then appends the nucleotide to the right.
func (b *BaseKmer) AppendBase(k *Kmer, nuc byte) {
	*k = b.mask & (*k << 2)
	switch nuc {
	case 'A', 'a':
		// k |= 0
	case 'C', 'c':
		*k |= 1
	case 'G', 'g':
		*k |= 2
	case 'T', 't':
		*k |= 3
	}
}

// String returns a string representation of the Kmer.
func (b *BaseKmer) String(k Kmer) string {
	if k == InvalidKmer {
		return strings.Repeat("-", b.k)
	}
	r := ""
	for i := b.k - 1; i >= 0; i-- {
		switch (k >> uint(i*2)) & 3 {
		case 0:
			r += "a"
		case 1:
			r += "c"
		case 2:
			r += "g"
		case 3:
			r += "t"
		}
	}
	return r
}
