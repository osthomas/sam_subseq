from sam_subseq.IndexMap import IndexMap


class SamRefAlignment:
    def __init__(self, seq, cigar, offset = 0, allow_oob = False):
        self._seq = seq
        # Initialize cursor position for stepping through the CIGAR string
        # The index map maps an index position on the reference to the
        # corresponding position(s) on the seq.
        self._index_map = IndexMap(cigar, offset, allow_oob)
        self.limits = self._index_map.limits

    def __len__(self):
        return len(self._index_map)

    def __getitem__(self, index):
        """
        Index with coordinates of the reference sequence and return the
        corresponding aligned sequence.

        Incomplete slices (:3, 3:5, :) are extended to the range of the *query*
        sequences. Explicitly requesting reference coordinates that are not
        covered by the alignment raises an IndexError.
        """
        if isinstance(index, int):
            start = index
            stop = start + 1
        elif isinstance(index, slice):
            start, stop, step = index.start, index.stop, index.step
            if step is not None:
                raise NotImplementedError("Stepped indexing is not implemented")
        else:
            raise TypeError(f"Indices must be integers or slices, not {type(index).__name__}")
        if start and start < 0 or stop and stop < 0:
            raise IndexError("Indices cannot be negative")

        qry_start, qry_stop = self._index_map.qry_range(start, stop)
        return self._seq[qry_start:qry_stop]
