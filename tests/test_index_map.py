from sam_subseq.IndexMap import IndexMap
import pytest

def test_parse_cigar():
    assert IndexMap("10M", 0)._cigar_ops == [(10, "M")]
    assert IndexMap("5H10M1I1D5S", 5)._cigar_ops == [
        (5, "H"),
        (10, "M"),
        (1, "I"),
        (1, "D"),
        (5, "S")
    ]
    # S only terminally or with H until end of CIGAR
    with pytest.raises(ValueError):
        IndexMap("5M10S5M", 0)
    # Unknown op
    with pytest.raises(ValueError):
        IndexMap("5M1K", 0)


def test_parse_cigar_no_count():
    assert IndexMap("MIDXSH", 0)._cigar_ops == [
        (1, "M"),
        (1, "I"),
        (1, "D"),
        (1, "X"),
        (1, "S"),
        (1, "H")
    ]


def test_length():
    assert len(IndexMap("10M", 0)) == 10
    assert len(IndexMap("10M10I10D", 0)) == 20
    assert len(IndexMap("10S10M", 10)) == 10


def test_limits():
    assert IndexMap("1M", 0).limits == (0, 1)
    assert IndexMap("10M", 0).limits == (0, 10)
    assert IndexMap("10S10M", 10).limits == (10, 20)
    assert IndexMap("10S10M20S", 10).limits == (10, 20)


def test_index_map_no_clips():
    assert IndexMap("3M", 0)._index_map == [(0, 1), (1, 2), (2, 3)]
    assert IndexMap("3M1D2M", 0)._index_map == [(0, 1), (1, 2), (2, 3), (3, 3), (3, 4), (4, 5)]

    # TC-----T
    # T-AAAAAT
    assert IndexMap("1M1D5I1M", 0)._index_map == [(0, 1), (1, 6), (6, 7)]

    # TCC--T
    # T--AAT
    assert IndexMap("1M2D2I1M", 0)._index_map == [(0, 1), (1, 1), (1, 3), (3, 4)]

    assert IndexMap("1M1D2M1I1M", 0)._index_map == [(0, 1), (1, 1), (1, 2), (2, 4), (4, 5)]


def test_index_map_softclips():
    assert IndexMap("2S2M", 2)._index_map == [(2, 3), (3, 4)]
    assert IndexMap("2S1M2S", 2)._index_map == [(2, 3)]


def test_index_map_hardclips():
    assert IndexMap("10H3M", 10)._index_map == [(0, 1), (1, 2), (2, 3)]
    assert IndexMap("10H3M10H", 10)._index_map == [(0, 1), (1, 2), (2, 3)]


def test_index_map_mixed_clips():
    assert IndexMap("10H4S1M", 10)._index_map == [(4, 5)]
    assert IndexMap("10H4S1M10S", 10)._index_map == [(4, 5)]


def test_index_map_oob_ranges():
    a = IndexMap("10S10M10S", offset = 10, allow_oob = True)
    # Full oob left
    assert a.qry_range(0, 5) == (10, 10)
    # Full oob right
    assert a.qry_range(30, 35) == (19, 19)
    # Half oob left
    assert a.qry_range(5, 10) == (10, 10)
    assert a.qry_range(9, 10) == (10, 10)
    assert a.qry_range(10, 10) == (10, 10)
    # in range
    assert a.qry_range(10, 11) == (10, 11)
    assert a.qry_range(10, 20) == (10, 20)
    # Half oob right
    assert a.qry_range(10, 25) == (10, 20)
    assert a.qry_range(11, 21) == (11, 20)
    assert a.qry_range(11, 25) == (11, 20)
    # oob left and right
    assert a.qry_range(5, 25) == (10, 20)
    assert a.qry_range(5, 20) == (10, 20)
    assert a.qry_range(10, 25) == (10, 20)
