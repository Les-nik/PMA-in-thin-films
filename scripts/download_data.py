"""Download the MuMax3 result files used in Fig. 1(e, f) and Fig. 3(b, c).

The files (about 53 MB) are attached to the GitHub release ``data-v1`` of this
repository and are saved to ``data/release/``.

    python scripts/download_data.py
"""

from pathlib import Path
from urllib.request import urlretrieve

BASE = "https://github.com/Les-nik/PMA-in-thin-films/releases/download/data-v1/"
FILES = [
    "IFC_1500MHz.npz",
    "bireflection_1500MHz_25deg.npz",
    "thick_20nm_beam3_ABL_65.5deg_130rad_9GHz_v2.npz",
]
TARGET = Path(__file__).resolve().parents[1] / "data" / "release"


def main():
    TARGET.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        dest = TARGET / name
        if dest.exists():
            print(f"{name}: already there")
            continue
        print(f"{name}: downloading")
        urlretrieve(BASE + name, dest)


if __name__ == "__main__":
    main()
