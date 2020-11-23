import os


def getfilenames(chr: str, gz: bool = False):
    file = os.path.abspath(__file__)
    folder = os.path.join(os.path.split(file)[0], "data")
    postfix = "-gz" if gz else ""
    extension = ".bed.gz" if gz else ".bed"
    return f"{folder}/labels{postfix}/{chr}{extension}", f"{folder}/peaks{postfix}/{chr}{extension}"
