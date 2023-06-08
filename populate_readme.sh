#!/usr/bin/env bash

# Run example and populate README

set -euo pipefail


# Replace a placeholder with file contents
function replace() {
    sed -e "/$1/{r $2
        d
    }"
}

rootdir="$(realpath "$(dirname "$0")")"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cp -r "$rootdir/doc/." "$tmp"
cd "$tmp"

minimap2 -a -s1 -m1 -w1 -E1,0 refs.fa queries.fa \
    | samtools sort -O sam \
    > example.sam

# .bam and .bai for tview
samtools view -O bam example.sam > example.bam
samtools index example.bam


for ref in $(samtools view example.sam | cut -f3 | uniq); do
    echo "<< Alignments to $ref >>" >> "$tmp/tview"
    samtools tview -d text -p "$ref" example.bam refs.fa \
        | cut -c1-67 \
        | awk 'NR != 3' \
        >> "$tmp/tview"
    echo >> "$tmp/tview"
done

sam_subseq example.sam --gff example.gff > example_out.fasta

# Make sure generated files match test case
for f in "example.sam" "example_out.fasta"; do
    diff "$f" "$rootdir/tests/resources/$f" &> /dev/null
done

cat README_template.md \
    | replace "{{README_USAGE}}" <(sam_subseq -h) \
    | replace "{{EXAMPLE_SAM}}" "example.sam"  \
    | replace "{{EXAMPLE_TVIEW}}" "$tmp/tview" \
    | replace "{{EXAMPLE_GFF}}" "example.gff" \
    | replace "{{EXAMPLE_OUTPUT}}" "example_out.fasta"
