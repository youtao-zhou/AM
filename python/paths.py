"""Project paths for the upstream analysis scripts."""

import os
from pathlib import Path

REPO = Path(os.environ.get("AM_REPO_DIR", ".")).resolve()
DATA = Path(os.environ.get("AM_DATA_DIR", REPO / "data")).resolve()


def data_path(*parts):
    return DATA.joinpath(*parts)


def code_path(*parts):
    return REPO.joinpath(*parts)
