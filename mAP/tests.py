import os
import tempfile

import numpy as np

from mAP import iou, Interval, VisualPeak, parse, PredictedPeak, match, average_precision

EPS = 1e-6

# Global test data
# it is ok to use """strings with tabulation""" here, bed parser strips spaces before processing
LABELS = {
    "chr1": """
    chr1	0	20  peakStart
    chr1	20	60	peaks
    chr1	60	80	peakEnd
    chr1	80	150	noPeaks
    chr1	150	160	peakStart
    chr1	160	170	peaks
    chr1	170	180	peakEnd
    chr1	180	190	noPeaks
    chr1	190	200	peakStart
    chr1	200	210	peaks
    chr1	210	220	peakEnd
    chr1	220	240	noPeaks
    """,

    "chr2": """
    chr2	0	10	peakStart
    chr2	10	100	peaks
    chr2	100	110	peakEnd
    chr2	110	160	noPeaks
    chr2	160	170	peakStart
    chr2	170	210	peaks
    chr2	210	220	peakEnd
    chr2	220	250	noPeaks
    chr2	250	260	peakStart
    chr2	260	310	peaks
    chr2	310	320	peakEnd
    """,

    "chr3": """
    chr3	0	10	peakStart
    chr3	10	60	peaks
    chr3	60	70	peakEnd
    chr3	70	100	noPeaks
    chr3	100	110	peakStart
    chr3	110	160	peaks
    chr3	160	170	peakEnd
    chr3	170	200	noPeaks
    chr3	240	280	noPeaks
    chr3	280	290	peakStart
    chr3	290	340	peaks
    chr3	340	350	peakEnd
    """,

    "chr4": """
    chr4	0	60	noPeaks
    chr4	90	140	noPeaks
    chr4	140 150	peakStart
    chr4	150	230	peaks
    chr4	230	240	peakEnd
    """,

    "chr5": """
    chr5	0   10	peakStart
    chr5	10  80	peaks
    chr5	80  90	peakEnd
    chr5	100 110 peakStart
    chr5	110 180 peaks
    chr5	180 190 peakEnd
    """,

    "chr6": """
    chr6	10  60  peakStart
    chr6	60  90  peaks
    chr6	90  100 peakEnd
    chr6	100 140 noPeaks
    chr6	140 170 peakStart
    chr6	170 220 peaks
    chr6	220 230 peakEnd
    """,

    "chr7": """
    chr7	10  20  peakStart
    chr7	20  80  peaks
    chr7	80  110 peakEnd
    chr7	110 140 peakStart
    chr7	140 200 peaks
    chr7	200 210 peakEnd
    chr7	330 340 peakStart
    chr7	340 360 peaks
    chr7	360 370 peakEnd
    """,

    "chr8": """
    chr8    0  70  peakStart
    chr8    70  80  peaks
    chr8    80  140 peakEnd
    chr8    150 180 noPeaks
    chr8    190 220 peakStart
    chr8    220 290 peaks
    chr8    290 320 peakEnd
    chr8    330 340 peakStart
    chr8    340 350 peaks
    chr8    350 360 peakEnd
    """,
}

PEAKS = {
    "chr1": """
    chr1	10	70	.	1
    chr1	140	230	.	2
    """,

    "chr2": """
    chr2	10	40	.	1
    chr2	50	80	.	2
    chr2	90	110	.	3
    chr2	120	130	.	4
    chr2	150	210	.	5
    chr2	270	300	.	6
    """,

    "chr3": """
    chr3	5   80	.	1
    chr3	90	170 .	2
    chr3	190 250	.	3
    chr3	290	370	.	4
    """,

    "chr4": """
    chr4	40  95	.	1
    chr4	220 280 .	2
    """,

    "chr5": """
    chr5	10  110	.	1
    chr5	120 310 .	2
    """,

    "chr6": """
    chr6	0   70  .	1
    chr6	130 150 .	2
    chr6	150 160 .	3
    chr6	160 230 .	4
    """,

    "chr7": """
    chr7	0   130  .	1
    chr7	130 190 .	2
    chr7	220 260 .	3
    chr7	270 320 .	4
    """,

    "chr8": """
    chr8	0   10  .	1
    chr8	20  40 .	2
    chr8	100 120 .	3
    chr8	130 180 .	4
    chr8	200 380 .	5
    """,
}


def to_file_and_parse(peaks: str, labels: str):
    fd, peaks_file = tempfile.mkstemp()
    with open(peaks_file, 'w') as file:
        file.write(peaks)
    os.close(fd)

    fd, labels_file = tempfile.mkstemp()
    with open(labels_file, 'w') as file:
        file.write(labels)
    os.close(fd)

    peaks, labels = parse(peaks_file, labels_file, score_column=4, presorted=False)
    os.remove(labels_file)
    os.remove(peaks_file)
    return peaks, labels


def test_iou():
    chrom, startlbl, bodylbl, endlbl = 'chr1', 'peakStart', 'peaks', 'peakEnd'

    start = Interval(chrom, 90, 110, startlbl)
    body = Interval(chrom, 110, 190, bodylbl)
    end = Interval(chrom, 190, 240, endlbl)
    peak = VisualPeak(start, body, end)

    # normal situations
    # ---start/end is inside the corresponding regions
    interval = Interval(chrom, 100, 200, '')
    assert abs(iou(interval, peak) - 1.0) < EPS

    # ---end is out of the corresponding region
    interval = Interval(chrom, 100, 300, '')
    assert abs(iou(interval, peak) - 80 / (80 + 60)) < EPS

    # ---start is out of the corresponding region
    interval = Interval(chrom, 10, 240, '')
    assert abs(iou(interval, peak) - 80 / (80 + 80)) < EPS

    # no intersection
    interval = Interval(chrom, 10, 90, '')
    assert iou(interval, peak) < EPS

    interval = Interval(chrom, 241, 262, '')
    assert iou(interval, peak) < EPS

    # inside start
    interval = Interval(chrom, 90, 110, '')
    assert iou(interval, peak) < EPS

    interval = Interval(chrom, 95, 105, '')
    assert iou(interval, peak) < EPS

    # inside end
    interval = Interval(chrom, 190, 240, '')
    assert iou(interval, peak) < EPS

    interval = Interval(chrom, 200, 230, '')
    assert iou(interval, peak) < EPS

    # peak inside the interval
    interval = Interval(chrom, 0, 400, '')
    assert abs(iou(interval, peak) - body.length / (interval.length - start.length - end.length)) < EPS

    interval = Interval(chrom, start.start, end.end, '')
    assert abs(iou(interval, peak) - 1.0) < EPS

    # interval starts inside the peak's body
    interval = Interval(chrom, 110, 400, '')
    assert abs(iou(interval, peak) - body.length / (interval.length - end.length)) < EPS

    interval = Interval(chrom, 160, 400, '')
    assert abs(iou(interval, peak) - 30 / 240) < EPS

    interval = Interval(chrom, 190, 400, '')
    assert iou(interval, peak) < EPS

    # interval ends inside the peak's body
    interval = Interval(chrom, 0, 190, '')
    assert abs(iou(interval, peak) - body.length / (interval.length - start.length)) < EPS

    interval = Interval(chrom, 0, 160, '')
    assert abs(iou(interval, peak) - 50 / 170) < EPS

    interval = Interval(chrom, 0, 110, '')
    assert iou(interval, peak) < EPS

    # interval is inside the peak
    interval = Interval(chrom, 110, 190, '')
    assert abs(iou(interval, peak) - 1) < EPS

    interval = Interval(chrom, 160, 170, '')
    assert abs(iou(interval, peak) - interval.length / body.length) < EPS


def test_parse():
    labels = LABELS['chr2'] + LABELS['chr1']
    peaks = PEAKS['chr2'] + PEAKS['chr1']

    peaks, labels = to_file_and_parse(peaks, labels)

    # recreate the above peaks here by hand
    ########################################
    handlabels, handpeaks = [], []

    label = VisualPeak(
        Interval('chr1', 0, 20, 'peakStart'),
        Interval('chr1', 20, 60, 'peaks'),
        Interval('chr1', 60, 80, 'peakEnd')
    )
    handlabels.append(label)

    handpeaks.append(
        PredictedPeak(((label, 1),), 0, 1, Interval('chr1', 10, 70, '.', 1))
    )

    ########################################

    handlabels.append(Interval('chr1', 80, 150, 'noPeaks'))
    label1 = VisualPeak(
        Interval('chr1', 150, 160, 'peakStart'),
        Interval('chr1', 160, 170, 'peaks'),
        Interval('chr1', 170, 180, 'peakEnd'),
    )
    label2 = VisualPeak(
        Interval('chr1', 190, 200, 'peakStart'),
        Interval('chr1', 200, 210, 'peaks'),
        Interval('chr1', 210, 220, 'peakEnd'),
    )
    handlabels.extend([
        label1,
        Interval('chr1', 180, 190, 'noPeaks'),
        label2,
        Interval('chr1', 220, 240, 'noPeaks')
    ])

    handpeaks.append(PredictedPeak(
        ((label1, 10 / 70), (label2, 10 / 70)), 30 / 90, 60 / 90, Interval('chr1', 140, 230, '.', 2)
    ))

    ########################################

    label = VisualPeak(
        Interval('chr2', 0, 10, 'peakStart'),
        Interval('chr2', 10, 100, 'peaks'),
        Interval('chr2', 100, 110, 'peakEnd'),
    )
    handlabels.append(label)

    handpeaks.extend([
        PredictedPeak(((label, 30 / 90),), 0, 1, Interval('chr2', 10, 40, '.', 1)),
        PredictedPeak(((label, 30 / 90),), 0, 1, Interval('chr2', 50, 80, '.', 2)),
        PredictedPeak(((label, 10 / 90),), 0, 1, Interval('chr2', 90, 110, '.', 3)),
    ])

    ########################################

    handlabels.append(Interval('chr2', 110, 160, 'noPeaks'))
    handpeaks.append(PredictedPeak(
        tuple(), 1, 0, Interval('chr2', 120, 130, '.', 4)
    ))

    ########################################

    label = VisualPeak(
        Interval('chr2', 160, 170, 'peakStart'),
        Interval('chr2', 170, 210, 'peaks'),
        Interval('chr2', 210, 220, 'peakEnd'),
    )
    handlabels.append(label)

    handpeaks.append(PredictedPeak(
        ((label, 4 / 5),), 1 / 6, 5 / 6, Interval('chr2', 150, 210, '.', 5)
    ))

    ########################################

    handlabels.append(Interval('chr2', 220, 250, 'noPeaks'))
    label = VisualPeak(
        Interval('chr2', 250, 260, 'peakStart'),
        Interval('chr2', 260, 310, 'peaks'),
        Interval('chr2', 310, 320, 'peakEnd'),
    )
    handlabels.append(label)

    handpeaks.append(PredictedPeak(
        ((label, 3 / 5),), 0, 1, Interval('chr2', 270, 300, '.', 6)
    ))

    ########################################

    # validate labels
    assert set(handlabels) == set(labels)

    # validate peaks
    assert set(handpeaks) == set(peaks)


def test_match():
    def assertions(peaks, labels, iouthr, extp, exfp, exfn):
        peaks, labels = to_file_and_parse(peaks, labels)
        labels = [lbl for lbl in labels if isinstance(lbl, VisualPeak)]
        for iouthr, extp, exfp, exfn in zip(iouthr, extp, exfp, exfn):
            tp, fp, fn = match(peaks, labels, iouthr)
            assert len(fn) == 0 and exfn == len(labels) - sum(tp) or fn[-1] == exfn
            assert sum(tp) == extp and sum(fp) == exfp

    # chr1
    assertions(
        PEAKS['chr1'], LABELS['chr1'],
        iouthr=[0.001, 0.5, 0.75, 0.95],
        extp=[2, 1, 1, 1],
        exfp=[0, 1, 1, 1],
        exfn=[1, 2, 2, 2]
    )

    # chr2
    assertions(
        PEAKS['chr2'], LABELS['chr2'],
        iouthr=[0.001, 0.5, 0.75, 0.95],
        extp=[3, 2, 1, 0],
        exfp=[3, 4, 5, 6],
        exfn=[0, 1, 2, 3]
    )

    # chr3
    assertions(
        PEAKS['chr3'], LABELS['chr3'],
        iouthr=[0.001, 0.5, 0.75, 0.95],
        extp=[3, 3, 2, 0],
        exfp=[0, 0, 2, 4],
        exfn=[0, 0, 1, 3]
    )

    # chr4
    assertions(
        PEAKS['chr4'], LABELS['chr4'],
        iouthr=[0.001, 0.5, 0.6, 0.75],
        extp=[1, 0, 0, 0],
        exfp=[0, 0, 1, 2],
        exfn=[0, 1, 1, 1]
    )

    # chr5
    assertions(
        PEAKS['chr5'], LABELS['chr5'],
        iouthr=[0.001, 0.35, 0.6, 0.78],
        extp=[2, 1, 1, 0],
        exfp=[0, 0, 0, 2],
        exfn=[0, 1, 1, 2]
    )

    # chr6
    assertions(
        PEAKS['chr6'], LABELS['chr6'],
        iouthr=[0.001, 0.3, 0.99],
        extp=[2, 1, 1],
        exfp=[2, 3, 3],
        exfn=[0, 1, 1]
    )

    # chr7
    assertions(
        PEAKS['chr7'], LABELS['chr7'],
        iouthr=[0.001, 0.75, 0.85],
        extp=[2, 1, 0],
        exfp=[0, 1, 2],
        exfn=[1, 2, 3]
    )

    # chr8
    assertions(
        PEAKS['chr8'], LABELS['chr8'],
        iouthr=[0.001, 0.75],
        extp=[1, 0],
        exfp=[3, 5],
        exfn=[2, 3]
    )


def test_average_precision():
    iou_thresholds = np.linspace(.5, 0.95, int(np.round((0.95 - .5) / .05)) + 1, endpoint=True)
    recall_thresholds = np.linspace(.0, 1.00, int(np.round((1.00 - .0) / .01)) + 1, endpoint=True)

    def assertions(peaks, labels, AP):
        peaks, labels = to_file_and_parse(peaks, labels)
        labels = [lbl for lbl in labels if isinstance(lbl, VisualPeak)]
        calc_AP = average_precision(peaks, labels, iou_thresholds=iou_thresholds, recall_tresholds=recall_thresholds)
        assert abs(calc_AP - AP) < EPS, f"{calc_AP} != expected {AP}"

    # chr1
    # 34 * 0.5 / 101 => 34 is a recall points that map to the precision 0.5. 
    # Everything is averaged for total recall points(101) and IOU(here of no use, PRC is the same for every threshold).
    assertions(PEAKS['chr1'], LABELS['chr1'], 34 * 0.5 / 101)

    # chr2
    assertions(
        PEAKS['chr2'], LABELS['chr2'],
        (
            # iou [0.5, 0.6]
            (67 * 1) * 3 +
            # iou [0.65, 0.8]
            (34 * 0.5) * 4
        ) / (101 * 10)
    )

    # chr3
    assertions(
        PEAKS['chr3'], LABELS['chr3'],
        (
            # iou [0.5, 0.65]
            (101 * 1) * 4 +
            # iou 0.7
            (34 * 1 + 67 * 0.75) +
            # iou [0.75, 0.8]
            (67 * 0.5) * 2
        ) / (101 * 10)
    )

    # chr4
    assertions(
        PEAKS['chr4'], LABELS['chr4'], 0
    )

    # chr5
    assertions(
        PEAKS['chr5'], LABELS['chr5'],
        (
            # iou [0.5, 0.6]
            (51 * 1) * 3 +
            # iou [0.65, 0.75]
            (51 * 0.5) * 3
        ) / (101 * 10)
    )

    # chr6
    assertions(
        PEAKS['chr6'], LABELS['chr6'],
        (
            # iou [0.5, 0.95]
            (51 * 1) * 10
        ) / (101 * 10)
    )

    # chr7
    assertions(
        PEAKS['chr7'], LABELS['chr7'],
        (
            # iou [0.5, 0.65]
            (67 * 1) * 4 +
            # iou [0.7, 0.8]
            (34 * 1) * 3
        ) / (101 * 10)
    )

    # chr8
    assertions(
        PEAKS['chr8'], LABELS['chr8'],
        (
            # iou 0.5
            (34 * 1)
        ) / (101 * 10)
    )


if __name__ == '__main__':
    test_iou()
    test_parse()
    test_match()
    test_average_precision()
