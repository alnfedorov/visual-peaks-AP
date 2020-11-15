import numpy as np

from utils import getfilenames
from visual_peaks_AP.AP import average_precision, VisualPeak, parse


IOU_THRESHOLDS = np.linspace(.5, 0.95, int(np.round((0.95 - .5) / .05)) + 1, endpoint=True)
RECALL_THRESHOLDS = np.linspace(.0, 1.00, int(np.round((1.00 - .0) / .01)) + 1, endpoint=True)


def assert_ap(labels_file: str, peaks_file: str, AP: float):
    eps = 1e-6

    peaks, labels = parse(peaks_file, labels_file, score_column=4, presorted=False)
    labels = [lbl for lbl in labels if isinstance(lbl, VisualPeak)]
    calc_AP = average_precision(peaks, labels, iou_thresholds=IOU_THRESHOLDS, recall_tresholds=RECALL_THRESHOLDS)
    assert abs(calc_AP - AP) < eps


def test_average_precision_chr1():
    # 34 * 0.5 / 101 => for iou threshold = 0.5, 34 recall points map to the precision 0.5, the rest are 0.
    # Average precision is averaged for total recall points(101) and IOU thresholds
    # IOU thresholds are not used here, PRC is the same for every threshold).
    assert_ap(
        *getfilenames('chr1'),
        AP=34 * 0.5 / 101
    )


def test_average_precision_chr2():
    assert_ap(
        *getfilenames('chr2'),
        AP=(
            # iou [0.5, 0.6]
            (67 * 1) * 3 +
            # iou [0.65, 0.8]
            (34 * 0.5) * 4
        ) / (101 * 10)
    )


def test_average_precision_chr3():
    assert_ap(
        *getfilenames('chr3'),
        AP=(
            # iou [0.5, 0.65]
            (101 * 1) * 4 +
            # iou 0.7
            (34 * 1 + 67 * 0.75) +
            # iou [0.75, 0.8]
            (67 * 0.5) * 2
        ) / (101 * 10)
    )


def test_average_precision_chr4():
    assert_ap(
        *getfilenames('chr4'),
        AP=0
    )


def test_average_precision_chr5():
    assert_ap(
        *getfilenames('chr5'),
        AP=(
            # iou [0.5, 0.6]
            (51 * 1) * 3 +
            # iou [0.65, 0.75]
            (51 * 0.5) * 3
        ) / (101 * 10)
    )


def test_average_precision_chr6():
    assert_ap(
        *getfilenames('chr6'),
        AP=(
            # iou [0.5, 0.95]
            (51 * 1) * 10
        ) / (101 * 10)
    )


def test_average_precision_chr7():
    assert_ap(
        *getfilenames('chr7'),
        AP=(
            # iou [0.5, 0.65]
            (67 * 1) * 4 +
            # iou [0.7, 0.8]
            (34 * 1) * 3
        ) / (101 * 10)
    )


def test_average_precision_chr8():
    assert_ap(
        *getfilenames('chr8'),
        AP=(
            # iou 0.5
            (34 * 1)
        ) / (101 * 10)
    )
