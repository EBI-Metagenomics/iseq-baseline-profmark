import configparser
import errno
import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

__all__ = ["root_dir", "cache_dir", "downsample", "label", "min_hmmer_hits", "Per"]


class Per(Enum):
    organism = 1
    profile = 2
    clan = 3


config = configparser.ConfigParser()
config.read("config.cfg")

root_dir = config["iseq_baseline_profmark"]["root_dir"]
cache_dir = Path(config["iseq_baseline_profmark"]["cache_dir"])
downsample = int(config["iseq_baseline_profmark"]["downsample"])
min_hmmer_hits = int(config["iseq_baseline_profmark"]["min_hmmer_hits"])


@dataclass(frozen=True)
class Label:
    recall = "recall (sensitivity, true positive rate)"
    precision = "precision (positive predictive value)"
    auc = "AUC"
    hits = "# HMMER3 hits"
    f1score = "F1 score"
    hmmer_hits = "# HMMER3 hits"
    log_evalue = "-log10(e-value)"


label = Label()


if not Path(root_dir).exists():
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), root_dir)
