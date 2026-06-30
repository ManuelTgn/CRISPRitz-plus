"""tests/core/test_genome_io.py — unit tests for genome_io.py"""

from __future__ import annotations

import pytest

from crispritz_plus.genome_io import GenomeReader, GenomeWriter, INDELOFFSET
from crispritz_plus.crispritz_errors import GenomeReaderError, GenomeWriterError


# ===========================================================================
# Constants
# ===========================================================================


class TestConstants:
    def test_indeloffset_is_50(self):
        assert INDELOFFSET == 50


# ===========================================================================
# GenomeReader
# ===========================================================================


class TestGenomeReaderRead:
    def test_reads_header_and_sequence(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        r.read()
        assert r.header == "chr1"
        assert len(r.sequence) == 100

    def test_sequence_uppercased(self, tmp_path):
        p = tmp_path / "lc.fa"
        p.write_text(">chrX\nacgt\n")
        r = GenomeReader(str(p), debug=True)
        r.read()
        assert r.sequence == list("ACGT")

    def test_sequence_and_enriched_copy_equal_after_read(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        r.read()
        assert r.sequence == r.sequence_enr

    def test_multiline_sequence_joined(self, tmp_path):
        p = tmp_path / "multi.fa"
        p.write_text(">chr1\nACGT\nACGT\n")
        r = GenomeReader(str(p), debug=True)
        r.read()
        assert len(r.sequence) == 8

    def test_empty_file_raises(self, tmp_path):
        p = tmp_path / "empty.fa"
        p.write_text("")
        with pytest.raises(GenomeReaderError):
            GenomeReader(str(p), debug=True).read()

    def test_missing_header_raises(self, tmp_path):
        p = tmp_path / "noheader.fa"
        p.write_text("ACGT\n")
        with pytest.raises(GenomeReaderError):
            GenomeReader(str(p), debug=True).read()

    def test_nonexistent_file_raises(self):
        with pytest.raises((GenomeReaderError, Exception)):
            GenomeReader("/nonexistent/path.fa", debug=True).read()

    def test_fname_property(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        assert r.fname == fasta_file

    def test_repr_before_read(self, fasta_file):
        assert "GenomeReader" in repr(GenomeReader(fasta_file, debug=False))


class TestGenomeReaderInsertSnp:
    def test_snp_modifies_enriched_not_reference(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        r.read()
        original = r.sequence[0]
        r.insert_snp("R", 0)
        assert r.sequence_enr[0] == "R"
        assert r.sequence[0] == original  # pristine unchanged

    def test_snp_at_arbitrary_position(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        r.read()
        r.insert_snp("Y", 50)
        assert r.sequence_enr[50] == "Y"


class TestGenomeReaderInsertIndel:
    def _long_reader(self, tmp_path, length=300):
        seq = "ACGT" * (length // 4)
        p = tmp_path / "long.fa"
        p.write_text(f">chrL\n{seq}\n")
        r = GenomeReader(str(p), debug=True)
        r.read()
        return r

    def test_indel_returns_indel_pair(self, tmp_path):
        r = self._long_reader(tmp_path)
        pair = r.insert_indel("ACG", pos=100, offset=1)
        assert hasattr(pair, "refseq")
        assert hasattr(pair, "indelseq")

    def test_refseq_length_uses_indeloffset(self, tmp_path):
        r = self._long_reader(tmp_path)
        offset = 2
        pair = r.insert_indel("AT", pos=100, offset=offset)
        expected_len = INDELOFFSET * 2 + offset
        assert len(pair.refseq) == expected_len

    def test_indelseq_splices_new_bases(self, tmp_path):
        r = self._long_reader(tmp_path)
        pair = r.insert_indel("TTT", pos=100, offset=1)
        # indel bases at positions [INDELOFFSET : INDELOFFSET+3]
        assert pair.indelseq[INDELOFFSET : INDELOFFSET + 3] == list("TTT")


class TestGenomeReaderToString:
    def test_to_string_joins_enriched_sequence(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        r.read()
        s = r.to_string()
        assert isinstance(s, str)
        assert len(s) == 100

    def test_to_string_reflects_snp_insertion(self, fasta_file):
        r = GenomeReader(fasta_file, debug=True)
        r.read()
        r.insert_snp("R", 0)
        assert r.to_string()[0] == "R"


# ===========================================================================
# GenomeWriter
# ===========================================================================


class TestGenomeWriter:
    def test_write_creates_file(self, tmp_path):
        p = tmp_path / "out.fa"
        GenomeWriter(str(p), debug=True).write("chr1", list("ACGT"))
        assert p.is_file()

    def test_written_header_prefixed(self, tmp_path):
        p = tmp_path / "out.fa"
        GenomeWriter(str(p), debug=True).write("chr2", list("ACGT"))
        assert p.read_text().startswith(">chr2\n")

    def test_written_sequence_correct(self, tmp_path):
        p = tmp_path / "out.fa"
        GenomeWriter(str(p), debug=True).write("chr1", list("GGCC"))
        lines = p.read_text().splitlines()
        assert lines[1] == "GGCC"

    def test_write_to_unwritable_path_raises(self):
        with pytest.raises((GenomeWriterError, Exception)):
            GenomeWriter("/nonexistent/dir/out.fa", debug=True).write("c", list("A"))

    def test_repr_and_str(self, tmp_path):
        w = GenomeWriter(str(tmp_path / "x.fa"), debug=False)
        assert "GenomeWriter" in repr(w)
        assert "GenomeWriter" in str(w)
