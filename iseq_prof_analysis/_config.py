import configparser
import errno
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# from enum import Enum

__all__ = [
    "load_config",
    "config",
    "Config",
    "ConfigBaseline",
    "ConfigChlamydia",
    "Label",
]


# class Per(Enum):
#     organism = 1
#     profile = 2
#     clan = 3


@dataclass(frozen=True)
class Label:
    recall = "recall (sensitivity, true positive rate)"
    precision = "precision (positive predictive value)"
    auc = "AUC"
    hits = "# HMMER3 hits"
    f1score = "F1 score"
    hmmer_hits = "# HMMER3 hits"
    log_evalue = "-log10(e-value)"


@dataclass
class ConfigBaseline:
    root_dir: Path
    cache_dir: Path
    downsample: int
    min_hmmer_hits: int


@dataclass
class ConfigChlamydia:
    root_dir: Path
    hybrid_consensus: Path


@dataclass
class Config:
    label = Label()
    baseline: Optional[ConfigBaseline] = None
    chlamydia: Optional[ConfigChlamydia] = None


config = Config()


def load_config(filepath: Optional[Path] = None, verbose=False):
    if filepath is None:
        root = Path(os.environ["XDG_CONFIG_HOME"])
        filepath = root / "iseq-prof-analysis" / "config.cfg"

    if verbose:
        print(f"Loading {filepath}.")

    _load_config(filepath)

    if verbose:
        print(config.baseline)
        print(config.chlamydia)


def _load_config(filepath: Path):
    global config

    cfg = configparser.ConfigParser()
    cfg.read(filepath)

    root_dir = Path(cfg["baseline"]["root_dir"])
    cache_dir = Path(cfg["baseline"]["cache_dir"])
    downsample = int(cfg["baseline"]["downsample"])
    min_hmmer_hits = int(cfg["baseline"]["min_hmmer_hits"])

    if not Path(root_dir).exists():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), root_dir)

    config.baseline = ConfigBaseline(root_dir, cache_dir, downsample, min_hmmer_hits)

    root_dir = Path(cfg["chlamydia"]["root_dir"])
    hybrid_consensus = Path(cfg["chlamydia"]["hybrid_consensus"])

    if not Path(root_dir).exists():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), root_dir)

    config.chlamydia = ConfigChlamydia(root_dir, hybrid_consensus)
