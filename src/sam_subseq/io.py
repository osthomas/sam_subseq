import sys
import io


def stdin_or_fh(f):
    """
    Read a line from stdin or a file on disk.
    """
    if f is sys.stdin:
        for line in f:
            yield line
    elif isinstance(f, str):
        with open(f, "r") as fh:
            for line in fh:
                yield line
    elif isinstance(f, io.TextIOWrapper):
        for line in f:
            yield line
    else:
        raise TypeError(f"Unrecognized input file of type {f.__class__.__name__}")


def parse_sam(src):
    """
    Read a SAM file.

    The one-based offset position is converted to zero-based.

    Parameters:
        src : Either sys.stdin or a fpath to a file.

    Yields:
        record : dict
            An alignment record from the SAM file, with the keys:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
            'pos' is converted to zero-based!

    Raises:
        ValueError if the SAM file is not sorted by coordinate.
    """
    is_sorted = False
    for line in stdin_or_fh(src):
        line = line.strip()
        fields = line.split("\t")
        if line.startswith("@"):
            # Still in the header block, no sort order defined
            if not is_sorted and line.startswith("@HD"):
                for field in fields:
                    is_sorted = "SO:coordinate" in fields
                    continue
        else:
            # Not in the header block
            if not is_sorted:
                raise ValueError("SAM file must be sorted by coordinate!")
            else:
                record = {
                    "qname": fields[0],
                    "flag": fields[1],
                    "rname": fields[2],
                    # Convert one-based position to zero-based
                    "pos": int(fields[3]) - 1,
                    "mapq": fields[4],
                    "cigar": fields[5].upper(),
                    "rnext": fields[6],
                    "pnext": fields[7],
                    "tlen": int(fields[8]),
                    "seq": fields[9],
                    "qual": fields[10]
                }
                yield record


def parse_gff(src):
    is_gff = False
    for line in stdin_or_fh(src):
        line = line.strip()
        if line.startswith("#"):
            # In the comment block
            if not is_gff:
                is_gff = "gff-version 3" in line
                continue
        else:
            # Not in the comment block
            fields = line.split("\t")
            attributes = {}
            for attribute in fields[8].split(";"):
                k, v = [x.strip() for x in attribute.split("=")]
                attributes[k] = v
            record = {
                "seqid": fields[0],
                "source": fields[1],
                "type": fields[2],
                # One-based to zero-based
                "start": int(fields[3]) - 1,
                # GFF stop is one-based and inclusive. Python indexing is
                # zero-based and end-exclusive -> no adjustment necesary
                "end": int(fields[4]),
                "score": fields[5],
                "strand": fields[6],
                "phase": fields[7],
                "attributes": attributes
            }
            yield record


def gff_to_header(gff_record, delim):
    """
    Format fields of a gff record for inclusion as a FASTX header.
    """
    header = [
        ("gff_id", gff_record["seqid"]),
        ("gff_type", gff_record["type"]),
        ("gff_start", gff_record["start"]),
        ("gff_end", gff_record["end"]),
        ("gff_phase", gff_record["phase"])
    ]
    if "Name" in gff_record["attributes"]:
        header.append(("gff_name", gff_record["attributes"]["Name"]))
    else:
        header.append(("gff_name", "NA"))
    header = ([f"{k}={v}" for k, v in header])
    header = ";".join([x for x in header])
    return header


def write_fasta(dst, header, seq):
    dst.write(f">{header}\n{seq}\n")
