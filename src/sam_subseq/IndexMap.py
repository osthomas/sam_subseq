import re


class IndexMap:
    """
    Build an index map, mapping reference coordinates to query coordinates.
    This allows retrieval of mapped read segments using reference positions.
    The index map is a list of tuples. Each list element corresponds to a
    reference position. The tuple at this list element contains the interval of
    the query positions corresponding to this reference position.

    Example:
    Reference sequence: AATTA
    Query sequence: ATTCA
    CIGAR: 1M1D2M1I1M

    The resulting alignment:
    REF: AATT-A
    QRY: A-TTCA
    CIG: 1M1D2M1I1M

    The resulting index map:
    [(0, 1), (1, 1), (1, 2), (2, 4), (4, 5)].

    ref[0] : "A" -> qry[0:1] : "A"   (match)
    ref[1] : "A" -> qry[1:1] : ""    (deletion)
    ref[2] : "T" -> qry[1:2] : "T"   (match)
    ref[3] : "T" -> qry[2:4] : "TC"  (insertion)
    ref[4] : "A" -> qry[4:5] : "A"   (match)
    """

    # Which CIGAR operations consume reference/query bases?
    C_CONSUME_REF = "MDN=X"
    C_CONSUME_QRY = "MIS=X"
    C_CONSUME_NONE = "HP"
    C_IS_CLIP = "HS"

    def __init__(self, cigar, offset, allow_oob = False):
        """
        Initiliaze an index map.

        Parameters:
            cigar : str
                The cigar string for the alignment.
            offset : int
                The start position of the alignment in the reference.
                This corresponds to the SAM field 'POS' minus 1.
            allow_oob : bool
                Allow indexing with out-of-bounds limits? If this is False,
                requesting a reference range that is not covered by the query
                is an error. If it is True, the out-of-bounds parts will be
                empty, and no error is raised.
        """
        self._cigar = cigar
        self._offset = offset
        self.allow_oob = allow_oob
        self._parse_cigar(cigar)
        # The reference range covered
        self.limits = (self._offset, len(self) + self._offset)
        self._qry_cur = 0
        self._index_map = []
        self.make_index_map()

    def qry_range(self, ref_start, ref_stop):
        """
        Return the range of the query sequence corresponding to specified
        reference range.

        Parameters:
            ref_start : int
                First reference position to include for the query
                (0-based, inclusive)
            ref_stop : int
                Last reference position for the query (0-based, exclusive)

        Returns:
            range : tuple
                Range (qry_start, qry_end) that corresponds to the requested
                reference range.
                Suitable to be used as qry_seq[qry_start:qry_end]
        """
        if ref_start is None:
            ref_start = self._offset
        if ref_stop is None:
            ref_stop = len(self) + self._offset
        if ref_start < 0 or ref_stop < 0:
            raise IndexError("Indices cannot be negative")
        if ref_start > ref_stop:
            raise IndexError("Start cannot be larger than stop")
        width = ref_stop - ref_start
        start = ref_start - self._offset
        if width > 0:
            stop = start + width - 1
        else:
            stop = start
        if not self.allow_oob:
            # Bounds checking
            for want in (start, stop):
                if want < 0 or want > len(self):
                    msg = "Requested reference range "
                    msg += f"{ref_start} - {ref_stop}, "
                    msg += f"but only {self.limits[0]} - {self.limits[1]} "
                    msg += "is covered by the query"
                    raise IndexError(msg)
        else:
            # Fix bounds
            # Both ends out of bounds to one side:
            # REF                     -----------------
            # QRY        xxxxxxxxxxx                     xxxxxxxxxxxx
            if start < 0 and stop < 0:
                # fully oob to the left
                # Return empty interval with earliest position of query
                interval = (self[0][0], self[0][0])
                print(interval)
                return (self[0][0], self[0][0])
            elif start >= len(self) and stop >= len(self):
                # fully oob to the right
                # Return empty interval with latest position of query
                return(self[-1][0], self[-1][0])
            # Clamp bounds to maxima
            # REF       ------------------------
            # QRY               xxxxxxxxxxx
            else:
                if start < 0:
                    # half oob left
                    start = 0
                if stop >= len(self):
                    # half oob right
                    stop = len(self) - 1

        if width == 0:
            return (self[start][0], self[start][0])
        elif width == 1:
            return (self[start][0], self[start][1])
        else:
            return (self[start][0], self[stop][1])

    def __getitem__(self, index):
        return self._index_map[index]

    def __len__(self):
        """
        The length of the index map is the length of the underlying reference
        sequence, ie. the sum of all cigar operations that consume the reference.
        """
        length = 0
        for count, op in self._cigar_ops:
            if op in self.C_CONSUME_REF:
                length += count
        return length

    def __repr__(self):
        return repr(self._index_map)

    def __str__(self):
        return str(self._index_map)

    def _parse_cigar(self, cigar):
        """
        Decompose a CIGAR string into its components, ie. count and operation
        """
        split_on = re.compile("([^0-9])")
        cigar = re.sub(split_on, r"\1 ", cigar).strip().split()
        out = []
        for c in cigar:
            count, op = re.split(split_on, c)[:2]
            if not count:
                count = 1
            count = int(count)
            out.append((count, op))
        # Sanity checks
        ops = [x[1] for x in out]
        for i, op in enumerate(ops):
            if op == "H" and i != 0 and i != len(ops) - 1:
                raise ValueError("Invalid CIGAR: H is only allowed first or last")
            if op == "S":
                if i != 0 and i != len(out)  - 1 and (ops[i - 1] != "H" and ops[i + 1] != "H"):
                    raise ValueError("Invalid CIGAR: S can only be preceded/followed by H")
            if op not in set(self.C_CONSUME_QRY + self.C_CONSUME_REF + self.C_CONSUME_NONE):
                raise ValueError(f"Invalid CIGAR: Unknown operation '{op}'")

        self._cigar_ops = out

    def advance_ref(self, by):
        # Place zero-length intervals at advanced ref indices
        addon = [(self._qry_cur, self._qry_cur)] * by
        self._index_map += addon

    def advance_qry(self, by, is_clip = False):
        # Extend interval of qry at current ref index
        # If there is no previous corresponding match position (ie, we're still
        # consuming clipping operations), just advance the query position.
        if is_clip:
            self._qry_cur += by
            return
        before = self._index_map[-1]
        after = (before[0], before[1] + by)
        self._index_map[-1] = after
        self._qry_cur += by

    def advance_both(self, by):
        addon = [(i, i + 1) for i in range(self._qry_cur, self._qry_cur + by)]
        self._index_map += addon
        self._qry_cur += by

    def advance(self, cigar_op):
        count, op = cigar_op
        is_clip = op in self.C_IS_CLIP
        if op in self.C_CONSUME_REF and op in self.C_CONSUME_QRY:
            self.advance_both(by = count)
        elif op in self.C_CONSUME_QRY:
            self.advance_qry(by = count, is_clip = is_clip)
        elif op in self.C_CONSUME_REF:
            self.advance_ref(by = count)
        elif op in self.C_CONSUME_NONE:
            pass
        else:
            raise ValueError(f"Unexpected CIGAR operation: '{op}'")

    def make_index_map(self, until = None):
        for cigar_op in self._cigar_ops:
            self.advance(cigar_op)
