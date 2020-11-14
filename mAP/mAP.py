import logging
import os
from dataclasses import dataclass
from itertools import chain
from typing import Tuple, List, Union, Optional

import numpy as np

logger = logging.getLogger(__file__)
DEFAULT_IOU_THRESHOLDS = np.linspace(.5, 0.95, int(np.round((0.95 - .5) / .05)) + 1, endpoint=True)
DEFAULT_RECALL_THRESHOLDS = np.linspace(.0, 1.00, int(np.round((1.00 - .0) / .01)) + 1, endpoint=True)


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    name: str
    score: Optional[float] = None

    def __post_init__(self):
        assert 0 <= self.start < self.end, f"Incorrect interval: {self}"

    @staticmethod
    def from_bed_file(path: str, score_column: Optional[int] = None) -> List["Interval"]:
        intervals = []
        with open(path, 'r') as file:
            for line in file:
                line = line.strip().split()
                if len(line) == 0:
                    continue
                if score_column is not None and score_column < len(line):
                    score = float(line[score_column])
                else:
                    score = None
                try:
                    new_int = Interval(line[0], int(line[1]), int(line[2]), line[3], score)
                    intervals.append(new_int)
                except AssertionError:
                    print(f"Failed to parse {line}")
                    continue
        return intervals

    @property
    def length(self):
        return self.end - self.start


@dataclass(frozen=True)
class VisualPeak:
    istart: Interval
    ibody: Interval
    iend: Interval

    def __post_init__(self):
        assert self.istart.end == self.ibody.start and self.ibody.end == self.iend.start and \
               self.istart.chrom == self.ibody.chrom == self.iend.chrom and \
               (self.istart.name, self.ibody.name, self.iend.name) == ("peakStart", "peaks", "peakEnd"), \
            f"Incorrect labels are found, the following intervals were expected to form a visual peak: \n" \
            f"\tStart: {self.istart}\n" \
            f"\tBody:  {self.ibody}\n" \
            f"\tEnd:   {self.iend}\n"

    @property
    def start(self):
        return self.istart.start

    @property
    def end(self):
        return self.iend.end

    @property
    def chrom(self):
        return self.istart.chrom


@dataclass(frozen=True)
class PredictedPeak:
    iou: Tuple[Tuple[VisualPeak, float], ...]
    bckg_fraction: float  # noPeaks / peak area
    peaks_fraction: float  # peakStart+peaks+peakEnd area / peak area
    interval: Interval

    def __post_init__(self):
        assert len(self.iou) > 0 or self.bckg_fraction >= 0
        if self.iou:
            assert self.peaks_fraction > 0


def iou(peak: Interval, label: VisualPeak) -> float:
    assert peak.chrom == label.chrom
    # ignore start and end intervals during IOU calculation
    intersection = min(peak.end, label.ibody.end) - max(peak.start, label.ibody.start)
    union = max(peak.end, label.ibody.end) - min(peak.start, label.ibody.start)

    intersection_start = min(peak.end, label.istart.end) - max(peak.start, label.istart.start)
    intersection_end = min(peak.end, label.iend.end) - max(peak.start, label.iend.start)

    iou = max(0, intersection) / (union - max(0, intersection_start) - max(0, intersection_end))
    assert iou >= 0
    return iou


def parse(peaks: str, labels: str, score_column: int, presorted: bool):
    assert os.path.exists(peaks), f"Peaks file {peaks} is not found"
    assert os.path.exists(labels), f"Labels file {labels} is not found"

    # load in memory
    labels, peaks = Interval.from_bed_file(labels, score_column), Interval.from_bed_file(peaks, score_column)

    if not presorted:
        peaks, labels = sorted(peaks, key=lambda x: x.start), sorted(labels, key=lambda x: x.start)
        peaks, labels = sorted(peaks, key=lambda x: x.chrom), sorted(labels, key=lambda x: x.chrom)

    # parse and validate labels
    parsed = []
    lblind = 0
    while lblind < len(labels):
        lbl = labels[lblind]
        # ensure there are no intersections between labels, they are in order
        if lblind != 0:
            prev = labels[lblind - 1]
            assert prev.chrom != lbl.chrom or prev.end <= lbl.start
        if lbl.name == "noPeaks":
            parsed.append(lbl)
            lblind += 1
        elif lbl.name == "peakStart":
            assert lblind + 2 < len(labels), f"Not finished peak found: {lbl}"
            parsed.append(VisualPeak(lbl, labels[lblind + 1], labels[lblind + 2]))
            lblind += 3
        else:
            raise Exception(f"Unexpected/Unknown label type {lbl.name}")
    labels = parsed

    # parse and validate peaks
    parsed = []
    numlabels, numpeaks = len(labels), len(peaks)
    lblind, peakind = 0, 0

    while peakind < numpeaks:
        curpeak = peaks[peakind]

        if len(parsed) != 0:
            prev = parsed[-1].interval
            in_order = prev.chrom < curpeak.chrom or \
                       prev.chrom == curpeak.chrom and prev.end <= curpeak.start
            if not in_order:
                logger.error(f"Passed peaks are not sorted: \n"
                             f"prev: {prev}\n"
                             f"cur: {curpeak}.\n"
                             f"It also might mean that peaks intersect with each other.\n"
                             f"Skipping {curpeak}...")
                peakind += 1
                continue

        while lblind < numlabels and (
                # skip to the current chromosome
                labels[lblind].chrom < curpeak.chrom or
                # skip to the current position
                labels[lblind].chrom == curpeak.chrom and labels[lblind].end <= curpeak.start
        ):
            lblind += 1

        if lblind >= numlabels:
            break

        visual_peaks_iou = []
        bckg_coverage = 0
        peaks_coverage = 0
        while lblind < numlabels and labels[lblind].chrom == curpeak.chrom and labels[lblind].start < curpeak.end:
            curlbl = labels[lblind]
            intersection = min(curpeak.end, curlbl.end) - max(curpeak.start, curlbl.start)
            assert intersection > 0

            if isinstance(curlbl, VisualPeak):
                visual_peaks_iou.append((curlbl, iou(curpeak, curlbl)))
                peaks_coverage += intersection
            elif isinstance(curlbl, Interval):
                assert curlbl.name == "noPeaks"
                bckg_coverage += intersection
            else:
                raise NotImplementedError()

            lblind += 1

        if peaks_coverage != 0 or bckg_coverage != 0:
            curpeak = PredictedPeak(iou=tuple(visual_peaks_iou),
                                    peaks_fraction=peaks_coverage / curpeak.length,
                                    bckg_fraction=bckg_coverage / curpeak.length, interval=curpeak)
            parsed.append(curpeak)

        peakind += 1
        # Is it really necessary? It is.
        lblind = max(0, lblind - 1)

    peaks = parsed
    logger.info(f"Successfully parsed {len(labels)} labels, {len(peaks)} peaks were covered by them.")
    return peaks, labels


def match(peaks: List[PredictedPeak], labels: List[VisualPeak], iouthr: float):
    matched = {lbl: None for lbl in labels}
    totaltp = 0
    tp, fp, fn = [], [], []
    for peak in peaks:
        matchlbl = None
        bestiou = iouthr
        for (lbl, lbliou) in peak.iou:
            # label is already matched
            if matched[lbl] is not None:
                continue
            # current match is better or iou is not high enough
            if lbliou < bestiou:
                continue
            bestiou = lbliou
            matchlbl = lbl

        if matchlbl is not None:
            matched[matchlbl] = peak
            totaltp += 1
            tp.append(1)
            fp.append(0)
            fn.append(len(labels) - totaltp)
        elif 1 - peak.peaks_fraction - peak.bckg_fraction < iouthr:
            # Predicted peak might be covered by unlabeled visual peak with the best IOU = 1 - bckg - labeled peaks area
            # conservatively assign peak as fp only if there is no way it is covered by the unlabeled data
            tp.append(0)
            fp.append(1)
            fn.append(len(labels) - totaltp)
    return tp, fp, fn


def average_precision(peaks: List[PredictedPeak], labels: List[Union[VisualPeak, Interval]],
                      iou_thresholds: np.ndarray, recall_tresholds: np.ndarray):
    peaks = sorted(peaks, key=lambda x: x.interval.score, reverse=True)
    labels = [lbl for lbl in labels if isinstance(lbl, VisualPeak)]
    average_precision = []

    for iouthr in iou_thresholds:
        tp, fp, fn = match(peaks, labels, iouthr)
        tp, fp, fn = np.cumsum(tp), np.cumsum(fp), np.asarray(fn)
        assert np.all(tp + fn == len(labels))

        recall = tp / (tp + fn)
        precision = tp / (tp + fp + np.spacing(1))

        for i in range(precision.size - 1, 0, -1):
            precision[i - 1] = max(precision[i - 1], precision[i])
        # The best precision is picked for each recall value (in case there are several of them)
        inds = np.searchsorted(recall, recall_tresholds, side='left')
        inds = inds[inds < precision.size]
        average_precision.append(precision[inds].sum())

    average_precision = sum(average_precision) / (recall_tresholds.size * iou_thresholds.size)
    return average_precision


def mAP(peaks: List[str], labels: List[str], classes: List[str],
        iou_thresholds: np.ndarray = DEFAULT_IOU_THRESHOLDS, recall_thresholds: np.ndarray = DEFAULT_RECALL_THRESHOLDS,
        score_column: int = 4, presorted: bool = False):
    # preprocess the data
    preprocessed = {cls: [] for cls in classes}
    for peaks, labels, cls in zip(peaks, labels, classes):
        peaks, labels = parse(peaks, labels, score_column, presorted)
        preprocessed[cls].append(
            (peaks, labels)
        )

    # ap for each class
    ap = {}
    for cls, preprocessed in preprocessed.items():
        # flatten peaks and labels, merge IOU dictionaries
        peaks, labels = zip(*preprocessed)
        peaks, labels = list(chain(*peaks)), list(chain(*labels))
        ap[cls] = average_precision(peaks, labels, iou_thresholds, recall_thresholds)
    mAP = sum(ap.values()) / len(ap)
    logger.info(f"mAP value is {mAP:.8f}")
    return mAP


if __name__ == '__main__':
    peaks = [
        "/data/encode/H3K4me3/ENCSR930HLX/original/peak-calling/replica[0]-homer-hard.broadPeak",
        "/data/encode/H3K4me3/ENCSR930HLX/original/peak-calling/replica[1]-homer-hard.broadPeak",

        "/data/encode/H3K4me3/ENCSR956CTX/original/peak-calling/replica[0]-homer-hard.broadPeak",
        "/data/encode/H3K4me3/ENCSR956CTX/original/peak-calling/replica[1]-homer-hard.broadPeak",

        "/data/encode/H3K4me3/ENCSR361FWQ/original/peak-calling/replica[0]-homer-hard.broadPeak",
        "/data/encode/H3K4me3/ENCSR361FWQ/original/peak-calling/replica[1]-homer-hard.broadPeak",
    ]

    # epic2 - 0.2686
    # macs2 - 0.25

    labels = [
        "/data/encode/H3K4me3/ENCSR930HLX/labels/labels.bed",
        "/data/encode/H3K4me3/ENCSR930HLX/labels/labels.bed",

        "/data/encode/H3K4me3/ENCSR956CTX/labels/labels.bed",
        "/data/encode/H3K4me3/ENCSR956CTX/labels/labels.bed",

        "/data/encode/H3K4me3/ENCSR361FWQ/labels/labels.bed",
        "/data/encode/H3K4me3/ENCSR361FWQ/labels/labels.bed",
    ]
    logging.basicConfig(level=logging.NOTSET)
    logging.info(f"The following thresholds are used:\n"
                 f"\tRecall:\n{DEFAULT_RECALL_THRESHOLDS}\n"
                 f"\tIOU:\n{DEFAULT_IOU_THRESHOLDS}")
    mAP(peaks, labels, ["" for _ in peaks], score_column=6)

    # macs2 0.599
    # epic2 0.539
    # idr + macs2 = 0.4
    # INFO:/home/mAP.py:Successfully parsed 941 labels, 5858 peaks were covered by them.
    # INFO:/home/mAP.py:mAP value is 0.599
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", '-l', nargs="+", required=True)
    parser.add_argument("--peaks", '-p', nargs="+", required=True)
    parser.add_argument("--score-column", default=6)
    # args = parser.parse_args()
    # mAP([args.peaks], [args.labels], ['1'])
