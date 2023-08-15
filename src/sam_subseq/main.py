#!/usr/bin/env pythonmain

"""
Extract reference features from aligned reads in a SAM file.
"""

import sys
import argparse
import textwrap

from sam_subseq import io
from sam_subseq.SamRefAlignment import SamRefAlignment


def parse_args():
    argparser = argparse.ArgumentParser(
        formatter_class = argparse.RawTextHelpFormatter,
        description = textwrap.dedent("""
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
        """)
    )
    argparser.add_argument(
        "infile", nargs = "?",
        type = argparse.FileType("r"), default = sys.stdin,
        help = "Input file (.sam). Default: stdin")
    argparser.add_argument(
        "outfile", nargs = "?",
        type = argparse.FileType("w"), default = sys.stdout,
        help = "Output file (.fasta) Default: stdout")
    argparser.add_argument(
        "--gff", required = True,
        help = "GFF files with features to extract. GFF SEQIDs (field 1) must "
        "correspond to SAM RNAMEs (field 1), or they will not be found."
    )
    args = argparser.parse_args()
    return args


def extract_features(record, record_features, header_delim = ";"):
    """
    Extract features from SAM sequence records.

    Yields:
        tuple : (header, subsequence)
            a header with information from the feature and sequence record and
            the extracted subsequence.
    """
    offset = record["pos"]
    cigar = record["cigar"]
    seq = record["seq"]
    alignment = SamRefAlignment(
        seq = seq,
        cigar = cigar,
        offset = offset,
        allow_oob = True)
    for feature in record_features:
        # Extract subsequence
        start, stop = feature["start"], feature["end"]
        qry_start, qry_stop = alignment._index_map.qry_range(start, stop)
        subseq = alignment[start:stop]
        # Prepare header for output
        header = record["qname"]
        if not header.endswith(header_delim):
            header += header_delim
        header += header_delim.join([
            # Feature information for record, ie. actual start/end positions
            f"qry_start={qry_start}",
            f"qry_stop={qry_stop}",
            f"qry_len={qry_stop - qry_start}",
            # Feature information from the gff
            io.gff_to_header(feature, delim = header_delim)
        ])
        yield (header, subseq)


def main(samfile, gfffile, outfile):
    features = {}
    # Store features for each reference sequence
    for feature in io.parse_gff(gfffile):
        if not feature["seqid"] in features.keys():
            features[feature["seqid"]] = []
        features[feature["seqid"]].append(feature)
    # Look for features and extract subsequences from alignments
    for record in io.parse_sam(samfile):
        try:
            record_features = features[record["rname"]]
            for header, subseq in extract_features(record, record_features):
                io.write_fasta(outfile, header, subseq)
        except KeyError:
            # No features to extract were found
            continue


def main_cmdline():
    args = parse_args()
    try:
        main(args.infile, args.gff, args.outfile)
    except BrokenPipeError:
        sys.exit(0)


if __name__ == "__main__":
    main()
