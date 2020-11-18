from typing import Tuple, List, Optional, Sequence


class Interval:
    __slots__ = ['chrom', 'start', 'end', 'name', 'score']

    def __init__(self, chrom: str, start: int, end: int, name: str = ".", score: Optional[float] = None):
        assert 0 <= start < end, f"Incorrect interval: [{start}, {end})"
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score

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

    def __repr__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.score}"

    @property
    def length(self):
        return self.end - self.start

    def __eq__(self, other):
        return isinstance(other, Interval) and self.chrom == other.chrom and self.start == other.start and \
               self.end == other.end and self.name == other.name and self.score == other.score

    def __hash__(self):
        return sum(hash(getattr(self, x)) for x in Interval.__slots__)


class VisualPeak:
    __slots__ = ['istart', 'ibody', 'iend']

    def __init__(self, istart: Interval, ibody: Interval, iend: Interval):
        assert istart.end == ibody.start and ibody.end == iend.start and \
               istart.chrom == ibody.chrom == iend.chrom and \
               (istart.name, ibody.name, iend.name) == ("peakStart", "peaks", "peakEnd"), \
            f"Incorrect labels are found, the following intervals were expected to form a visual peak: \n" \
            f"\tStart: {istart}\n" \
            f"\tBody:  {ibody}\n" \
            f"\tEnd:   {iend}\n"
        self.istart = istart
        self.ibody = ibody
        self.iend = iend

    @property
    def start(self):
        return self.istart.start

    @property
    def end(self):
        return self.iend.end

    @property
    def chrom(self):
        return self.istart.chrom

    def __eq__(self, other):
        return isinstance(other, VisualPeak) and \
               self.istart == other.istart and \
               self.ibody == other.ibody and \
               self.iend == other.iend

    def __hash__(self):
        return hash(self.istart) + hash(self.ibody) + hash(self.iend)


class PredictedPeak:
    __slots__ = ["iou", "bckg_fraction", "peaks_fraction", "interval"]
    Matches = Sequence[Tuple[VisualPeak, float]]

    def __init__(self, iou: Matches, bckg_fraction: float, peaks_fraction: float,
                 interval: Interval):
        assert len(iou) > 0 or bckg_fraction >= 0
        if iou:
            assert peaks_fraction > 0
        self.iou = iou
        self.bckg_fraction = bckg_fraction
        self.peaks_fraction = peaks_fraction
        self.interval = interval

    def __eq__(self, other):
        return isinstance(other, PredictedPeak) and \
               self.iou == other.iou and \
               self.bckg_fraction == other.bckg_fraction and \
               self.peaks_fraction == other.peaks_fraction and \
               self.interval == other.interval

    def __hash__(self):
        return hash(self.iou) + hash(self.bckg_fraction) + hash(self.peaks_fraction) + hash(self.interval)
