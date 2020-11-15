import os


def getfilenames(chr: str):
    file = os.path.abspath(__file__)
    folder = os.path.join(os.path.split(file)[0], "data")
    return f"{folder}/labels/{chr}.bed", f"{folder}/peaks/{chr}.bed"
