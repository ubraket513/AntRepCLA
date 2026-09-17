"""Loading and normalising IgBLAST output."""

import pytest

from antibody_repertoire.io.cache import load_cache, save_cache
from antibody_repertoire.io.igblast import load_igblast_results

HEADER = "v_call\td_call\tj_call\tcdr3\tcdr3_aa\tcdr3_start\tcdr3_end\n"


def write_tsv(path, rows):
    path.write_text(HEADER + "".join(rows))
    return path


def test_strips_allele_suffixes(tmp_path):
    """IGHV3-23*01 and IGHV3-23*02 are the same gene for lineage purposes."""
    tsv = write_tsv(tmp_path / "r.tsv", ["IGHV3-23*01\tIGHD1*01\tIGHJ4*02\tAAAA\tCAR\t1\t4\n"])
    df = load_igblast_results(tsv)

    assert df.loc[0, 'v_call'] == {'IGHV3-23'}
    assert df.loc[0, 'j_call'] == {'IGHJ4'}


def test_splits_ambiguous_calls_into_sets(tmp_path):
    tsv = write_tsv(tmp_path / "r.tsv", ["IGHV1-2*01,IGHV1-3*01\tIGHD1*01\tIGHJ4*02\tAAAA\tCAR\t1\t4\n"])
    df = load_igblast_results(tsv)

    assert df.loc[0, 'v_call'] == {'IGHV1-2', 'IGHV1-3'}


def test_cdr3_length_is_inclusive_of_both_endpoints(tmp_path):
    tsv = write_tsv(tmp_path / "r.tsv", ["IGHV1*01\tIGHD1*01\tIGHJ4*01\tAAAA\tCAR\t10\t13\n"])
    assert load_igblast_results(tsv).loc[0, 'cdr_length'] == 4


def test_coordinate_columns_are_dropped(tmp_path):
    tsv = write_tsv(tmp_path / "r.tsv", ["IGHV1*01\tIGHD1*01\tIGHJ4*01\tAAAA\tCAR\t1\t4\n"])
    df = load_igblast_results(tsv)

    assert 'cdr3_start' not in df.columns
    assert 'cdr3_end' not in df.columns


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_igblast_results(tmp_path / "absent.tsv")


def test_cache_round_trip(tmp_path):
    payload = {'a': 1, 'b': [2, 3]}
    path = tmp_path / "nested" / "c.pkl"
    save_cache(payload, path)
    assert load_cache(path) == payload


def test_missing_cache_returns_none(tmp_path):
    assert load_cache(tmp_path / "absent.pkl") is None


def test_corrupt_cache_returns_none_rather_than_raising(tmp_path):
    """An unreadable cache means recompute, not crash."""
    path = tmp_path / "bad.pkl"
    path.write_bytes(b"not a pickle")
    assert load_cache(path) is None
