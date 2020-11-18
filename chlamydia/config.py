import configparser
import errno
import os
from pathlib import Path

__all__ = ["root_dir", "hybrid_consensus", "orig_subdir"]


config = configparser.ConfigParser()
config.read("config.cfg")

root_dir = config["chlamydia"]["root_dir"]
orig_subdir = config["chlamydia"]["orig_subdir"]
hybrid_consensus = config["chlamydia"]["hybrid_consensus"]

if not Path(root_dir).exists():
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), root_dir)
