import pytest
from sam_subseq.SamRefAlignment import SamRefAlignment

def test_sam_ref_alignment_basic_indexing():
    a = SamRefAlignment("ATG", "3M")
    assert a[0] == "A"
    assert a[1] == "T"
    assert a[2] == "G"
    assert a[1:1] == ""
    assert a[0:1] == "A"
    assert a[0:3] == "ATG"
    assert a[:2] == "AT"
    assert a[:3] == "ATG"
    assert a[:] == "ATG"


def test_sam_ref_alignment_raises():
    a = SamRefAlignment("ATG", "3M")
    with pytest.raises(IndexError):
        a[4]
    with pytest.raises(IndexError):
        a[3:4]
    with pytest.raises(IndexError):
        a[4:5]
    with pytest.raises(TypeError):
        a["a"]
    with pytest.raises(NotImplementedError):
        a[::1]
    with pytest.raises(IndexError):
        a[-1]


def test_sam_ref_alignment_allow_oob():
    a = SamRefAlignment("atgccATGC", "5S4M", offset=5, allow_oob = True)
    assert a[5] == "A"
    assert a[:] == "ATGC"
    assert a[:7] == "AT"
    assert a[2:7] == "AT"
    assert a[5:] == "ATGC"
    assert a[5:6] == "A"
    assert a[5:7] == "AT"
    assert a[5:20] == "ATGC"
    assert a[0:20] == "ATGC"
    assert a[8:9] == "C"
    assert a[9:9] == ""
    # a spans 5-8
    # OOB left
    assert a[0:2] == ""
    # Half OOB left, just in
    assert a[4:6] == "A"
    assert a[5:6] == "A"
    assert a[6:7] == "T"
    assert a[7:8] == "G"
    assert a[8:9] == "C"
    # Half OOB right, just out
    assert a[8:10] == "C"
    assert a[8:11] == "C"
    # OOB right
    assert a[9:10] == ""
    assert a[10:10] == ""
    assert a[10:11] == ""
    with pytest.raises(IndexError):
        # Invalid indexing should not be rescued by allow_oob
        a[10:7]
    b = SamRefAlignment("ATGC", "4M", allow_oob = True)
    assert b[0] == "A"
    assert b[:] == "ATGC"
    assert b[:7] == "ATGC"
    assert b[0:3] == "ATG"
    assert b[5:7] == ""

def test_sam_ref_alignment_complex_indexing_no_offset():
    # AAATTTGGG CCC
    # AAATT GGGACCC
    a = SamRefAlignment("AAATTGGGACCC", cigar = "5MD3MI3M")
    assert a[:] == "AAATTGGGACCC"
    assert a[5] == ""
    assert a[8] == "GA"
    assert a[9] == "C"
    assert a[0:5] == "AAATT"
    assert a[0:6] == "AAATT"
    assert a[0:7] == "AAATTG"
    assert a[0:8] == "AAATTGG"
    assert a[0:9] == "AAATTGGGA"
    assert a[0:10] == "AAATTGGGAC"
    assert a[0:11] == "AAATTGGGACC"
    assert a[0:12] == "AAATTGGGACCC"


def test_sam_ref_alignment_complex_indexing_with_clipped():
    # xxxxAAATTTGGG CCC
    # atgcAAATT GGGACCC
    a = SamRefAlignment("atgcAAATTGGGACCC", cigar = "4S5MD3MI3M", offset = 4)
    assert a[4] == "A"
    assert a[:] == "AAATTGGGACCC"
    assert a[:7] == "AAA"
    assert a[:9] == "AAATT"
    assert a[6:13] == "ATTGGGA"
    #   xxATGC
    # aattATGC
    b = SamRefAlignment("aattATGC", cigar = "4S4M", offset = 0)
    assert b[0] == "A"


def test_sam_ref_alignment_complex_indexing_with_offset():
    #       6
    # xxxxxxxxxxAAATTTGGG CCC
    #       atgcAAATT GGGACCC
    #           0
    a = SamRefAlignment("atgcAAATTGGGACCC", cigar = "4S5MD3MI3M", offset = 10)
    assert a[10] == "A"
    assert a[:] == "AAATTGGGACCC"
    assert a[:11] == "A"
    assert a[:19] == "AAATTGGGA"
    assert a[12:19] == "ATTGGGA"
