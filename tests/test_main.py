import io
import pytest
from importlib import resources

from sam_subseq.main import main


@pytest.fixture
def sam():
    return str(resources.files("tests.resources").joinpath("example.sam"))

    # return pkg_resources.resource_filename("tests.resources", "example.sam")


@pytest.fixture
def gff():
    return str(resources.files("tests.resources").joinpath("example.gff"))


@pytest.fixture
def expect():
    fpath = str(resources.files("tests.resources").joinpath("example_out.fasta"))
    with open(fpath, "r") as f:
        expect = f.read()
    return expect


def test_main(sam, gff, expect):
    with io.StringIO() as outf:
        main(sam, gff, outf)
        assert outf.getvalue() == expect
