from typing import List

from utils import getfilenames
from visual_peaks_AP.ap import match, VisualPeak, parse


def assert_matching(labels_file: str, peaks_file: str,
                    iouthr: List[float], extp: List[int], exfp: List[int], exfn: List[int]):
    peaks, labels = parse(peaks_file, labels_file, score_column=4, presorted=False)
    labels = [lbl for lbl in labels if isinstance(lbl, VisualPeak)]
    for iouthr, extp, exfp, exfn in zip(iouthr, extp, exfp, exfn):
        tp, fp, fn = match(peaks, labels, iouthr)
        assert len(fn) == 0 and exfn == len(labels) - sum(tp) or fn[-1] == exfn
        assert sum(tp) == extp and sum(fp) == exfp


def test_matching_chr1():
    assert_matching(
        *getfilenames('chr1'),
        iouthr=[0.001, 0.5, 0.75, 0.95],
        extp=[2, 1, 1, 1],
        exfp=[0, 1, 1, 1],
        exfn=[1, 2, 2, 2]
    )


def test_matching_chr2():
    assert_matching(
        *getfilenames('chr2'),
        iouthr=[0.001, 0.5, 0.75, 0.95],
        extp=[3, 2, 1, 0],
        exfp=[3, 4, 5, 6],
        exfn=[0, 1, 2, 3]
    )


def test_matching_chr3():
    assert_matching(
        *getfilenames('chr3'),
        iouthr=[0.001, 0.5, 0.75, 0.95],
        extp=[3, 3, 2, 0],
        exfp=[0, 0, 2, 4],
        exfn=[0, 0, 1, 3]
    )


def test_matching_chr4():
    assert_matching(
        *getfilenames('chr4'),
        iouthr=[0.001, 0.5, 0.6, 0.75],
        extp=[1, 0, 0, 0],
        exfp=[0, 0, 1, 2],
        exfn=[0, 1, 1, 1]
    )


def test_matching_chr5():
    assert_matching(
        *getfilenames('chr5'),
        iouthr=[0.001, 0.35, 0.6, 0.78],
        extp=[2, 1, 1, 0],
        exfp=[0, 0, 0, 2],
        exfn=[0, 1, 1, 2]
    )


def test_matching_chr6():
    assert_matching(
        *getfilenames('chr6'),
        iouthr=[0.001, 0.3, 0.99],
        extp=[2, 1, 1],
        exfp=[2, 3, 3],
        exfn=[0, 1, 1]
    )


def test_matching_chr7():
    assert_matching(
        *getfilenames('chr7'),
        iouthr=[0.001, 0.75, 0.85],
        extp=[2, 1, 0],
        exfp=[0, 1, 2],
        exfn=[1, 2, 3]
    )


def test_matching_chr8():
    assert_matching(
        *getfilenames('chr8'),
        iouthr=[0.001, 0.75],
        extp=[1, 0],
        exfp=[3, 5],
        exfn=[2, 3]
    )
