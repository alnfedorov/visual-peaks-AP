import os
import tempfile

import pytest

from utils import getfilenames
from visual_peaks_AP.ap import Interval, VisualPeak, PredictedPeak, parse


# hand-parsed peaks/labels from the data/test_parsing folder
@pytest.fixture
def chr1():
    labels = [VisualPeak(
        Interval('chr1', 0, 20, 'peakStart'),
        Interval('chr1', 20, 60, 'peaks'),
        Interval('chr1', 60, 80, 'peakEnd')
    )]

    peaks = [
        PredictedPeak(((labels[0], 1),), 0, 1, Interval('chr1', 10, 70, '.', 1))
    ]

    lbl1 = VisualPeak(
        Interval('chr1', 150, 160, 'peakStart'),
        Interval('chr1', 160, 170, 'peaks'),
        Interval('chr1', 170, 180, 'peakEnd'),
    )
    lbl2 = VisualPeak(
        Interval('chr1', 190, 200, 'peakStart'),
        Interval('chr1', 200, 210, 'peaks'),
        Interval('chr1', 210, 220, 'peakEnd'),
    )
    labels.extend([
        Interval('chr1', 80, 150, 'noPeaks'),
        lbl1,
        Interval('chr1', 180, 190, 'noPeaks'),
        lbl2,
        Interval('chr1', 220, 240, 'noPeaks')
    ])

    peaks.append(PredictedPeak(
        ((lbl1, 10 / 70), (lbl2, 10 / 70)), 30 / 90, 60 / 90, Interval('chr1', 140, 230, '.', 2)
    ))

    return labels, peaks, getfilenames('chr1')


@pytest.fixture
def chr2():
    lbl = VisualPeak(
        Interval('chr2', 0, 10, 'peakStart'),
        Interval('chr2', 10, 100, 'peaks'),
        Interval('chr2', 100, 110, 'peakEnd'),
    )
    labels = [lbl]
    peaks = [
        PredictedPeak(((lbl, 30 / 90),), 0, 1, Interval('chr2', 10, 40, '.', 1)),
        PredictedPeak(((lbl, 30 / 90),), 0, 1, Interval('chr2', 50, 80, '.', 2)),
        PredictedPeak(((lbl, 10 / 90),), 0, 1, Interval('chr2', 90, 110, '.', 3)),
    ]

    labels.append(Interval('chr2', 110, 160, 'noPeaks'))
    peaks.append(PredictedPeak(tuple(), 1, 0, Interval('chr2', 120, 130, '.', 4)))

    lbl = VisualPeak(
        Interval('chr2', 160, 170, 'peakStart'),
        Interval('chr2', 170, 210, 'peaks'),
        Interval('chr2', 210, 220, 'peakEnd'),
    )
    peak = PredictedPeak(
        ((lbl, 4 / 5),), 1 / 6, 5 / 6, Interval('chr2', 150, 210, '.', 5)
    )
    labels.append(lbl)
    peaks.append(peak)

    lbl = VisualPeak(
        Interval('chr2', 250, 260, 'peakStart'),
        Interval('chr2', 260, 310, 'peaks'),
        Interval('chr2', 310, 320, 'peakEnd'),
    )
    labels.extend([lbl, Interval('chr2', 220, 250, 'noPeaks')])
    peaks.append(PredictedPeak(
        ((lbl, 3 / 5),), 0, 1, Interval('chr2', 270, 300, '.', 6)
    ))

    return labels, peaks, getfilenames('chr2')


@pytest.mark.parametrize('workload', [chr1, chr2])
def test_single_chr_parsing(workload, request):
    labels, peaks, (labels_file, peaks_file) = request.getfixturevalue(workload.__name__)

    parsed_peaks, parsed_labels = parse(peaks_file, labels_file, score_column=4, presorted=False)
    assert set(peaks) == set(parsed_peaks)
    assert set(labels) == set(parsed_labels)


def test_all_chr_parsing(chr1, chr2):
    labels, peaks = [], []
    labels_file, peaks_file = tempfile.mkstemp()[1], tempfile.mkstemp()[1]

    labels_stream, peaks_stream = open(labels_file, 'w'), open(peaks_file, 'w')
    for lbls, pks, (lfile, pfile) in [chr1, chr2]:
        labels.extend(lbls)
        peaks.extend(pks)
        for infile, outfile in [(lfile, labels_stream), (pfile, peaks_stream)]:
            with open(infile) as infile:
                outfile.write(infile.read() + "\n")
    labels_stream.close(), peaks_stream.close()

    parsed_peaks, parsed_labels = parse(peaks_file, labels_file, score_column=4, presorted=False)
    assert set(peaks) == set(parsed_peaks)
    assert set(labels) == set(parsed_labels)

    os.unlink(labels_file)
    os.unlink(peaks_file)
