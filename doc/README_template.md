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

{{README_USAGE}}

```


## Example

**SAM input**
```
{{EXAMPLE_SAM}}
```


**Alignment**
```
{{EXAMPLE_TVIEW}}
```


**GFF input**
```
{{EXAMPLE_GFF}}
```


**FASTA output**
```
{{EXAMPLE_OUTPUT}}
```
