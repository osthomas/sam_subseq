# sam_subseq - Extract GFF Features From Aligned Reads

`sam_subseq` takes two inputs:

1. SAM file with reads (or sequences in general) aligned to one or more references
2. GFF file defining features for the reference(s)

`sam_subseq` will project the GFF coordinates (which refer to the reference)
onto the reads, extract subsequences corresponding to the GFF features (with
deletion, insertions, truncations, etc.), and output these subsequences in
FASTA format.


## Installation

```bash
pip install sam_subseq
```


## Usage

```
$ sam_subseq -h

usage: sam_subseq [-h] --gff GFF [infile] [outfile]

Extract features (subsequences) from aligned reads in a SAM file,
using annotations for the reference sequence.

sam_subseq parses the CIGAR string to determine which part of the
read sequence (the query) to output.
The SAM file must be sorted by coordinate (default for samtools sort)!

Example:
             80        180        290
                       |---CDS----|
             |----exon------------|
REF:  -------------------------------
QRY:         xxxxxxxxxxyyyyyyy--z

The reference has an exon annotation from position 80-290.
Extracting this feature from the query will yield: xxxxxxxxxxyyyyyyyz
The CDS in the query shows a deletion and is incompletely represented.
Extracting the CDS from 180-290 will yield yyyyyyyz.

Some information from the gff file is written into the header of each
output sequence. Coordinates conform to Python conventions, ie.
zero-based and end-exclusive.

These fields are of the form 'label=value;'. Currently, the following
information is output:
- the original sequence header
- qry_start: The start coordinate of the extracted feature in the
  query (ie. aligned, non-reference sequence)
- qry_stop: The end coordinate of the extracted feature in the query
- qry_len: The length of the extracted feature in the query
  The length can be zero, for example if a feature spans positions
  50-100, but the alignment of the query spans only positions 10-40
- gff_id: The ID of the gff record
- gff_type: The type of the gff record
- gff_start: The start coordinate as defined in the GFF (ie. for the
  reference)
- gff_end: The end coordinate as defined in the GFF
- gff_phase: The phase as defined in the GFF
- gff_name: If a 'Name' annotation is present in the GFF attribute
  field, it is output. If it is not available, this is set to NA.

The output is a FASTA file with one extracted feature per record.

positional arguments:
  infile      Input file (.sam). Default: stdin
  outfile     Output file (.fasta) Default: stdout

options:
  -h, --help  show this help message and exit
  --gff GFF   GFF files with features to extract. GFF SEQIDs (field 1) must correspond to SAM RNAMEs (field 1), or they will not be found.

```


## Example

**SAM input**
```
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref1	LN:67
@SQ	SN:ref2	LN:67
@PG	ID:minimap2	PN:minimap2	VN:2.26-r1175	CL:minimap2 -a -s1 -m1 -w1 -E1,0 refs.fa queries.fa
@PG	ID:samtools	PN:samtools	PP:minimap2	VN:1.17	CL:samtools sort -O sam
qry1	0	ref1	1	60	67M	*	0	0	ATCGAGTCGTAGCAGGCTGAGCGATGCGAGGCAGCGACGGACGAGTAGCAGCTAAAGCTAAGGAGCA	*	NM:i:0	ms:i:134	AS:i:134	nn:i:0	tp:A:P	cm:i:53	s1:i:67	s2:i:0	de:f:0	rl:i:0
qry3	0	ref1	1	46	25M19D23M	*	0	0	ATCGAGTCGTAGCAGGCTGAGCGATGTAGCAGCTAAAGCTAAGGAGCA	*	NM:i:19	ms:i:88	AS:i:73	nn:i:0	tp:A:P	cm:i:21	s1:i:44	s2:i:0	de:f:0.0204	rl:i:0
qry2	0	ref1	33	48	35M	*	0	0	AGCGACGGACGAGTAGCAGCTAAAGCTAAGGAGCA	*	NM:i:0	ms:i:70	AS:i:70	nn:i:0	tp:A:P	cm:i:21	s1:i:35	s2:i:0	de:f:0	rl:i:0
qry4	0	ref2	29	55	39M	*	0	0	GAGCTGATGCACGACACGACGATCGATCGACTGTATGTA	*	NM:i:0	ms:i:78	AS:i:78	nn:i:0	tp:A:P	cm:i:25	s1:i:39	s2:i:0	de:f:0	rl:i:0
```


**Alignment**
```
<< Alignments to ref1 >>
1         11        21        31        41        51        61     
ATCGAGTCGTAGCAGGCTGAGCGATGCGAGGCAGCGACGGACGAGTAGCAGCTAAAGCTAAGGAGCA
...................................................................
.........................*******************.......................
                                ...................................

<< Alignments to ref2 >>
1         11        21        31        41        51        61     
ACGACGTACGTAGCGAACGACGATCGACGAGCTGATGCACGACACGACGATCGATCGACTGTATGTA
                            .......................................

```


**GFF input**
```
##gff-version 3
##sequence-region ref1 1 67
ref1	.	gene	1	67	.	+	.	ID=ref1
ref1	.	exon	10	62	.	+	.	ID=ref1:exon;=ref1-exon;Parent=ref1
ref1	.	CDS	20	62	.	+	0	ID=ref1:CDS;Name=ref1-cds;Parent=ref1
##sequence-region ref2 1 67
ref2	.	gene	1	67	.	+	.	ID=ref2
ref2	.	exon	10	62	.	+	.	ID=ref2:exon;=ref2-exon;Parent=ref2
ref2	.	CDS	20	62	.	+	0	ID=ref2:CDS;Name=ref2-cds;Parent=ref2
```


**FASTA output**
```
>qry1;qry_start=0;qry_stop=67;qry_len=67;gff_id=ref1;gff_type=gene;gff_start=0;gff_end=67;gff_phase=.;gff_name=NA
ATCGAGTCGTAGCAGGCTGAGCGATGCGAGGCAGCGACGGACGAGTAGCAGCTAAAGCTAAGGAGCA
>qry1;qry_start=9;qry_stop=62;qry_len=53;gff_id=ref1;gff_type=exon;gff_start=9;gff_end=62;gff_phase=.;gff_name=NA
TAGCAGGCTGAGCGATGCGAGGCAGCGACGGACGAGTAGCAGCTAAAGCTAAG
>qry1;qry_start=19;qry_stop=62;qry_len=43;gff_id=ref1;gff_type=CDS;gff_start=19;gff_end=62;gff_phase=0;gff_name=ref1-cds
AGCGATGCGAGGCAGCGACGGACGAGTAGCAGCTAAAGCTAAG
>qry3;qry_start=0;qry_stop=48;qry_len=48;gff_id=ref1;gff_type=gene;gff_start=0;gff_end=67;gff_phase=.;gff_name=NA
ATCGAGTCGTAGCAGGCTGAGCGATGTAGCAGCTAAAGCTAAGGAGCA
>qry3;qry_start=9;qry_stop=43;qry_len=34;gff_id=ref1;gff_type=exon;gff_start=9;gff_end=62;gff_phase=.;gff_name=NA
TAGCAGGCTGAGCGATGTAGCAGCTAAAGCTAAG
>qry3;qry_start=19;qry_stop=43;qry_len=24;gff_id=ref1;gff_type=CDS;gff_start=19;gff_end=62;gff_phase=0;gff_name=ref1-cds
AGCGATGTAGCAGCTAAAGCTAAG
>qry2;qry_start=0;qry_stop=35;qry_len=35;gff_id=ref1;gff_type=gene;gff_start=0;gff_end=67;gff_phase=.;gff_name=NA
AGCGACGGACGAGTAGCAGCTAAAGCTAAGGAGCA
>qry2;qry_start=0;qry_stop=30;qry_len=30;gff_id=ref1;gff_type=exon;gff_start=9;gff_end=62;gff_phase=.;gff_name=NA
AGCGACGGACGAGTAGCAGCTAAAGCTAAG
>qry2;qry_start=0;qry_stop=30;qry_len=30;gff_id=ref1;gff_type=CDS;gff_start=19;gff_end=62;gff_phase=0;gff_name=ref1-cds
AGCGACGGACGAGTAGCAGCTAAAGCTAAG
>qry4;qry_start=0;qry_stop=39;qry_len=39;gff_id=ref2;gff_type=gene;gff_start=0;gff_end=67;gff_phase=.;gff_name=NA
GAGCTGATGCACGACACGACGATCGATCGACTGTATGTA
>qry4;qry_start=0;qry_stop=34;qry_len=34;gff_id=ref2;gff_type=exon;gff_start=9;gff_end=62;gff_phase=.;gff_name=NA
GAGCTGATGCACGACACGACGATCGATCGACTGT
>qry4;qry_start=0;qry_stop=34;qry_len=34;gff_id=ref2;gff_type=CDS;gff_start=19;gff_end=62;gff_phase=0;gff_name=ref2-cds
GAGCTGATGCACGACACGACGATCGATCGACTGT
```
