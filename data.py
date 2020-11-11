from pathlib import Path
from typing import List, Optional, Tuple

import lttb
import numpy as np
import pandas as pd
from diskcache import FanoutCache
from iseq_prof import ConfusionMatrix, OrganismResult, Profiling, pfam
from tqdm import tqdm

import config

size_limit = 250 * 1024 ** 3
timeout = 10 * 60
cache = FanoutCache(directory=config.cache_dir, size_limit=size_limit, timeout=timeout)

EVALUES = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
EVALUES += [1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18]
EVALUES += [1e-19, 1e-20]


def progress(iterable, **kwargs):
    return tqdm(iterable, mininterval=5, ascii=True, **kwargs)


@cache.memoize()
def _fetch_organism(
    root_dir: str, organism: str
) -> Tuple[OrganismResult, ConfusionMatrix]:
    prof = Profiling(Path(root_dir).resolve())
    result = prof.read_organism_result(organism)
    result.accession.domain
    return result, result.confusion_matrix()


@cache.memoize()
def create_per_organism_df(root_dir: str, organisms: List[str], downsample: int):
    dfs = []
    for organism in progress(organisms):

        result, cm = _fetch_organism(root_dir, organism)

        pr = cm.pr_curve
        evalues = cm.sample_scores
        assert len(pr.recall) == len(evalues)

        x = pr.recall
        x, idx = np.unique(x, return_index=True)
        y = pr.precision[idx]
        evalues = evalues[idx]

        matrix = lttb.downsample(np.stack((x, y), axis=1), n_out=downsample)
        idx = np.searchsorted(x, matrix[:, 0])
        x = matrix[:, 0]
        y = matrix[:, 1]
        evalues = evalues[idx]

        df = pd.DataFrame()
        df[config.label.recall] = x
        df[config.label.precision] = y
        df["organism"] = organism
        df[config.label.auc] = pr.auc
        df[config.label.hmmer_hits] = cm.P
        df["domain"] = result.accession.domain
        df["e-value"] = evalues

        dfs.append(df)
    return pd.concat(dfs).reset_index(drop=True)


@cache.memoize()
def create_per_organism_df_evalue(root_dir: str, organisms: List[str]):
    dfs = []
    for organism in progress(organisms):

        result, cm = _fetch_organism(root_dir, organism)

        df = pd.DataFrame()
        df["e-value"] = EVALUES
        df[config.label.log_evalue] = -np.log10(EVALUES)
        df[config.label.f1score] = [cm.f1score[cm.cutpoint(e)] for e in EVALUES]

        recall = [cm.recall[cm.cutpoint(e)] for e in EVALUES]
        df[config.label.recall] = recall

        precision = [cm.precision[cm.cutpoint(e)] for e in EVALUES]
        df[config.label.precision] = precision

        df["organism"] = organism
        df["domain"] = result.accession.domain
        df[config.label.hmmer_hits] = cm.P
        df[config.label.auc] = cm.pr_curve.auc

        dfs.append(df)
    return pd.concat(dfs).reset_index(drop=True)


@cache.memoize(typed=True)
def _prof_confusion_matrices(root_dir: str, organisms: List[str], clan_wise: bool):
    prof = Profiling(Path(root_dir).resolve())
    return prof.confusion_matrix(organisms, clan_wise=clan_wise)


@cache.memoize()
def create_per_clan_df(root_dir: str, organisms: List[str]):
    dfs = []
    confusion_matrices = _prof_confusion_matrices(root_dir, organisms, True)
    for clan, cm in progress(confusion_matrices.items()):

        pr = cm.pr_curve
        if cm.P <= config.min_hmmer_hits:
            continue

        evalues = cm.sample_scores
        assert len(pr.recall) == len(evalues)

        x = pr.recall
        x, idx = np.unique(x, return_index=True)
        y = pr.precision[idx]
        evalues = evalues[idx]

        n_out = min(config.downsample, len(x))
        matrix = lttb.downsample(np.stack((x, y), axis=1), n_out=n_out)
        idx = np.searchsorted(x, matrix[:, 0])
        x = matrix[:, 0]
        y = matrix[:, 1]
        evalues = evalues[idx]

        df = pd.DataFrame()
        df[config.label.recall] = x
        df[config.label.precision] = y
        df["clan"] = clan
        df[config.label.auc] = pr.auc
        df[config.label.hmmer_hits] = cm.P
        df["e-value"] = evalues
        dfs.append(df)

    df = pd.concat(dfs)
    df.sort_values(
        [
            "clan",
            config.label.recall,
            config.label.precision,
        ],
        inplace=True,
    )
    df = df.reset_index(drop=True)

    return df


def _clan_name(name: Optional[str]) -> str:
    if name is None:
        return "Unclassified"
    return name


@cache.memoize()
def create_per_profile_df(root_dir: str, organisms: List[str]):
    dfs = []
    confusion_matrices = _prof_confusion_matrices(root_dir, organisms, False)
    clans = pfam.Clans()
    for profile, cm in progress(confusion_matrices.items()):

        pr = cm.pr_curve
        if cm.P <= config.min_hmmer_hits:
            continue

        evalues = cm.sample_scores
        assert len(pr.recall) == len(evalues)

        x = pr.recall
        x, idx = np.unique(x, return_index=True)
        y = pr.precision[idx]
        evalues = evalues[idx]

        n_out = min(config.downsample, len(x))
        matrix = lttb.downsample(np.stack((x, y), axis=1), n_out=n_out)
        idx = np.searchsorted(x, matrix[:, 0])
        x = matrix[:, 0]
        y = matrix[:, 1]
        evalues = evalues[idx]

        df = pd.DataFrame()
        df[config.label.recall] = x
        df[config.label.precision] = y
        df["profile"] = profile
        df["clan"] = _clan_name(clans.get(profile))
        df[config.label.auc] = pr.auc
        df[config.label.hmmer_hits] = cm.P
        df["e-value"] = evalues
        dfs.append(df)

    df = pd.concat(dfs)
    df.sort_values(
        [
            "profile",
            "clan",
            config.label.recall,
            config.label.precision,
        ],
        inplace=True,
    )
    df = df.reset_index(drop=True)

    return df


@cache.memoize()
def create_per_clan_df_evalue(root_dir: str, organisms: List[str]):
    dfs = []
    confusion_matrices = _prof_confusion_matrices(root_dir, organisms, True)
    for clan, cm in progress(confusion_matrices.items()):

        if cm.P <= config.min_hmmer_hits:
            continue

        df = pd.DataFrame()
        df["e-value"] = EVALUES
        df[config.label.log_evalue] = -np.log10(EVALUES)
        df[config.label.f1score] = [cm.f1score[cm.cutpoint(e)] for e in EVALUES]

        recall = [cm.recall[cm.cutpoint(e)] for e in EVALUES]
        df[config.label.recall] = recall

        precision = [cm.precision[cm.cutpoint(e)] for e in EVALUES]
        df[config.label.precision] = precision

        df[config.label.hmmer_hits] = cm.P
        df["clan"] = clan
        df[config.label.auc] = cm.pr_curve.auc

        dfs.append(df)
    return pd.concat(dfs).reset_index(drop=True)


@cache.memoize()
def create_per_profile_df_evalue(root_dir: str, organisms: List[str]):
    dfs = []
    confusion_matrices = _prof_confusion_matrices(root_dir, organisms, False)
    clans = pfam.Clans()
    for profile, cm in progress(confusion_matrices.items()):

        if cm.P <= config.min_hmmer_hits:
            continue

        df = pd.DataFrame()
        df["e-value"] = EVALUES
        df[config.label.log_evalue] = -np.log10(EVALUES)
        df[config.label.f1score] = [cm.f1score[cm.cutpoint(e)] for e in EVALUES]

        recall = [cm.recall[cm.cutpoint(e)] for e in EVALUES]
        df[config.label.recall] = recall

        precision = [cm.precision[cm.cutpoint(e)] for e in EVALUES]
        df[config.label.precision] = precision

        df[config.label.hmmer_hits] = cm.P
        df["profile"] = profile
        df["clan"] = _clan_name(clans.get(profile))
        df[config.label.auc] = cm.pr_curve.auc

        dfs.append(df)
    return pd.concat(dfs).reset_index(drop=True)
