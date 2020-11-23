import gzip
import logging
import os
from typing import List, Optional

from .dataclasses import Interval, PredictedPeak, VisualPeak

logger = logging.getLogger(__file__)


def iou(peak: Interval, label: VisualPeak) -> float:
    assert peak.chrom == label.chrom
    # ignore start and end intervals during IOU calculation
    intersection = min(peak.end, label.ibody.end) - max(peak.start, label.ibody.start)
    union = max(peak.end, label.ibody.end) - min(peak.start, label.ibody.start)

    ignore_start = min(peak.end, label.istart.end) - max(peak.start, label.istart.start)
    ignore_end = min(peak.end, label.iend.end) - max(peak.start, label.iend.start)

    iou = max(0, intersection) / (union - max(0, ignore_start) - max(0, ignore_end))
    assert iou >= 0
    return iou


def _parse_intervals(path: str, score_column: Optional[int] = None) -> List[Interval]:
    intervals = []
    if path.endswith(".gz"):
        logger.info(f"gzipped file detected: {path}")
        stream = gzip.open(path, 'rt', encoding='utf-8')
    else:
        logger.info(f"plain text file detected: {path}")
        stream = open(path, 'r')

    for line in stream:
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

    stream.close()

    return intervals


def _parse_labels(labels: List[Interval]):
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
            assert lblind + 2 < len(labels), f"Unfinished peak found: {lbl}"
            parsed.append(VisualPeak(lbl, labels[lblind + 1], labels[lblind + 2]))
            lblind += 3
        else:
            raise Exception(f"Unexpected/Unknown label type {lbl.name}")
    return parsed


def parse(peaks: str, labels: str, score_column: int, presorted: bool):
    assert os.path.exists(peaks), f"Peaks file {peaks} is not found"
    assert os.path.exists(labels), f"Labels file {labels} is not found"

    # load in memory
    labels, peaks = _parse_intervals(labels, score_column), _parse_intervals(peaks, score_column)

    if not presorted:
        peaks, labels = sorted(peaks, key=lambda x: x.start), sorted(labels, key=lambda x: x.start)
        peaks, labels = sorted(peaks, key=lambda x: x.chrom), sorted(labels, key=lambda x: x.chrom)

    # parse and validate labels
    labels = _parse_labels(labels)

    # parse and validate peaks
    parsed = []
    numlabels, numpeaks = len(labels), len(peaks)
    lblind, peakind = 0, 0

    while peakind < numpeaks:
        curpeak = peaks[peakind]

        if len(parsed) != 0:
            prev = parsed[-1].interval
            in_order = prev.chrom < curpeak.chrom or prev.chrom == curpeak.chrom and prev.end <= curpeak.start
            if not in_order:
                logger.error(f"Passed peaks are not sorted: \n"
                             f"previous: {prev}\n"
                             f"current:  {curpeak}.\n"
                             f"It also might mean that peaks intersect with each other.\n"
                             f"Skipping {curpeak}...")
                peakind += 1
                continue

        # skip labels
        while lblind < numlabels and (
                # to the current chromosome
                labels[lblind].chrom < curpeak.chrom or
                # to the current position
                labels[lblind].chrom == curpeak.chrom and labels[lblind].end <= curpeak.start
        ):
            lblind += 1

        if lblind >= numlabels:
            break

        # match all labels to the peak
        matches = []
        bckg_coverage = 0
        peaks_coverage = 0
        while lblind < numlabels and labels[lblind].chrom == curpeak.chrom and labels[lblind].start < curpeak.end:
            curlbl = labels[lblind]
            intersection = min(curpeak.end, curlbl.end) - max(curpeak.start, curlbl.start)
            assert intersection > 0

            if isinstance(curlbl, VisualPeak):
                matches.append((curlbl, iou(curpeak, curlbl)))
                peaks_coverage += intersection
            elif isinstance(curlbl, Interval):
                assert curlbl.name == "noPeaks"
                bckg_coverage += intersection
            else:
                raise NotImplementedError()

            lblind += 1

        if peaks_coverage != 0 or bckg_coverage != 0:
            curpeak = PredictedPeak(iou=tuple(matches),
                                    peaks_fraction=peaks_coverage / curpeak.length,
                                    bckg_fraction=bckg_coverage / curpeak.length, interval=curpeak)
            parsed.append(curpeak)

        peakind += 1
        # Rollback to the previous labels
        lblind = max(0, lblind - 2)

    peaks = parsed
    logger.info(f"Successfully parsed {len(labels)} labels, {len(peaks)} peaks are covered by them.")
    return peaks, labels
