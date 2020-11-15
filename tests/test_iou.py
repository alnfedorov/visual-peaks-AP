import pytest

from visual_peaks_AP.AP import Interval, VisualPeak, iou


@pytest.fixture
def setup():
    chrom, startlbl, bodylbl, endlbl = 'chr1', 'peakStart', 'peaks', 'peakEnd'
    start = Interval(chrom, 90, 110, startlbl)
    body = Interval(chrom, 110, 190, bodylbl)
    end = Interval(chrom, 190, 240, endlbl)
    peak = VisualPeak(start, body, end)
    return peak


def assert_iou(interval: Interval, peak: VisualPeak, iouval: float):
    eps = 1e-6
    assert abs(iou(interval, peak) - iouval) < eps


def test_iou(setup):
    peak = setup
    chrom = peak.chrom
    workload = [
        # expected
        (Interval(chrom, 100, 200), 1.0),  # start/end is inside the corresponding regions
        (Interval(chrom, 100, 300), 80 / (80 + 60)),  # end is out of the peak
        (Interval(chrom, 10, 240), 80 / (80 + 80)),  # start is out of the peak
        # start/end is out of the peak
        (Interval(chrom, 0, 400), peak.ibody.length / (400 - peak.istart.length - peak.iend.length)),

        # no intersection
        (Interval(chrom, 10, 90), 0.0),  # to the left
        (Interval(chrom, 241, 262), 0.0),  # to the right
        # inside start
        (Interval(chrom, 90, 110), 0.0),
        (Interval(chrom, 95, 105), 0.0),
        # inside end
        (Interval(chrom, 190, 240), 0.0),
        (Interval(chrom, 200, 230), 0.0),

        # completely inside the peak
        (peak.ibody, 1.0),
        (Interval(chrom, 160, 170), 10 / peak.ibody.length),

        # interval starts inside the peak's body
        (Interval(chrom, 110, 400), peak.ibody.length / (290 - peak.iend.length)),
        (Interval(chrom, 160, 400), 30 / 240),
        (Interval(chrom, 190, 444), 0),

        # interval ends inside the peak's body
        (Interval(chrom, 0, 190), peak.ibody.length / (190 - peak.istart.length)),
        (Interval(chrom, 0, 160), 50 / 170),
        (Interval(chrom, 0, 110), 0),
    ]
    for inter, iouval in workload:
        assert_iou(inter, peak, iouval)
