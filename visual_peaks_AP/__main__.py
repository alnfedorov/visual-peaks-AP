#
# if __name__ == '__main__':
#     peaks = [
#         "/data/encode/H3K4me3/ENCSR930HLX/original/peak-calling/replica[0]-homer-hard.broadPeak",
#         "/data/encode/H3K4me3/ENCSR930HLX/original/peak-calling/replica[1]-homer-hard.broadPeak",
#
#         "/data/encode/H3K4me3/ENCSR956CTX/original/peak-calling/replica[0]-homer-hard.broadPeak",
#         "/data/encode/H3K4me3/ENCSR956CTX/original/peak-calling/replica[1]-homer-hard.broadPeak",
#
#         "/data/encode/H3K4me3/ENCSR361FWQ/original/peak-calling/replica[0]-homer-hard.broadPeak",
#         "/data/encode/H3K4me3/ENCSR361FWQ/original/peak-calling/replica[1]-homer-hard.broadPeak",
#     ]
#
#     # epic2 - 0.2686
#     # macs2 - 0.25
#
#     labels = [
#         "/data/encode/H3K4me3/ENCSR930HLX/labels/labels.bed",
#         "/data/encode/H3K4me3/ENCSR930HLX/labels/labels.bed",
#
#         "/data/encode/H3K4me3/ENCSR956CTX/labels/labels.bed",
#         "/data/encode/H3K4me3/ENCSR956CTX/labels/labels.bed",
#
#         "/data/encode/H3K4me3/ENCSR361FWQ/labels/labels.bed",
#         "/data/encode/H3K4me3/ENCSR361FWQ/labels/labels.bed",
#     ]
#     logging.basicConfig(level=logging.NOTSET)
#     logging.info(f"The following thresholds are used:\n"
#                  f"\tRecall:\n{DEFAULT_RECALL_THRESHOLDS}\n"
#                  f"\tIOU:\n{DEFAULT_IOU_THRESHOLDS}")
#     AP(peaks, labels, ["" for _ in peaks], score_column=6)
#
#     # macs2 0.599
#     # epic2 0.539
#     # idr + macs2 = 0.4
#     # INFO:/home/mAP.py:Successfully parsed 941 labels, 5858 peaks were covered by them.
#     # INFO:/home/mAP.py:mAP value is 0.599
#     import argparse
#
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--labels", '-l', nargs="+", required=True)
#     parser.add_argument("--peaks", '-p', nargs="+", required=True)
#     parser.add_argument("--score-column", default=6)
#     # args = parser.parse_args()
#     # mAP([args.peaks], [args.labels], ['1'])
