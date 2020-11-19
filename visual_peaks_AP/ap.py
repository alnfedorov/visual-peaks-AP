import logging
from itertools import chain
from typing import List, Union, Optional

import numpy as np

from .dataclasses import Interval, PredictedPeak, VisualPeak
from .parsing import parse

logger = logging.getLogger(__file__)
DEFAULT_IOU_THRESHOLDS = np.linspace(.5, 0.95, int(np.round((0.95 - .5) / .05)) + 1, endpoint=True)
DEFAULT_RECALL_THRESHOLDS = np.linspace(.0, 1.00, int(np.round((1.00 - .0) / .01)) + 1, endpoint=True)


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
        elif (1 - peak.peaks_fraction - peak.bckg_fraction) < iouthr:
            # Predicted peak might be covered by unlabeled visual peak with the max IOU = 1 - bckg - labeled peaks area
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

        # Pick the best precision for each recall value
        for i in range(precision.size - 1, 0, -1):
            precision[i - 1] = max(precision[i - 1], precision[i])

        inds = np.searchsorted(recall, recall_tresholds, side='left')
        # It is not mandatory for the recall to achieve 1 here
        inds = inds[inds < precision.size]
        average_precision.append(precision[inds].sum())

    average_precision = sum(average_precision) / (recall_tresholds.size * iou_thresholds.size)
    return average_precision


def AP(peaks: List[str], labels: List[str], classes: Optional[List[str]] = None,
       iou_thresholds: np.ndarray = DEFAULT_IOU_THRESHOLDS, recall_thresholds: np.ndarray = DEFAULT_RECALL_THRESHOLDS,
       score_column: int = 4, presorted: bool = False):

    assert len(peaks) == len(labels), \
        f"Number of peak/labels files must be the same, got {len(peaks)} peaks and {len(labels)} labels"

    classes = classes if classes else [""] * len(peaks)
    assert len(classes) == len(peaks), \
        f"Each peak/label pair must have an associated class, got {len(peaks)} pairs and {len(classes)} classes"

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
    ap = sum(ap.values()) / len(ap)
    logger.info(f"AP value is {ap:.8f}")
    return ap
