from typing import Optional

import pandas as pd
import plotly.express as px
from pandas import DataFrame
from plotly.graph_objs import Figure

import config

__all__ = [
    "plot_auc_dist",
    "plot_auc_scatter",
    "plot_pr_curve_per_organism",
    "plot_pr_scatter",
    "plot_score_boxplot",
]

XRANGE = [-0.05, 1.05]
YRANGE = [-0.05, 1.05]
NBINS = 100


def plot_pr_curve_per_organism(df: DataFrame) -> Figure:
    # df = df[df[config.label.precision] > 0].reset_index(drop=True)
    fig = px.line(
        df,
        x=config.label.recall,
        y=config.label.precision,
        color="domain",
        line_group="organism",
        hover_data=[config.label.auc, "domain", config.label.hmmer_hits, "e-value"],
        hover_name="organism",
        title="Precision-Recall curves, organism-wise",
    )
    fig.update_xaxes(range=XRANGE)
    fig.update_yaxes(range=YRANGE)
    return fig


def plot_auc_dist(df: DataFrame, per: config.Per) -> Figure:
    """
    Plot AUC distribution.

    Parameters
    ----------
    df
        DataFrame.
    per
        Per what?
    """
    if per == config.Per.organism:
        color: Optional[str] = "domain"
        fields = ["organism", "domain"]
        hover_data = [config.label.auc, config.label.hmmer_hits]
    else:
        color = None
        hover_data = [per.name, "clan", config.label.auc, config.label.hmmer_hits]
        hover_data = list(set(hover_data))
        fields = [per.name, "clan", config.label.auc, config.label.hmmer_hits]
        fields = list(set(fields))

    title = f"AUC distribution, {per.name}-wise"
    df = df.drop_duplicates(fields)
    fig = px.histogram(
        df,
        x=config.label.auc,
        color=color,
        title=title,
        marginal="rug",
        hover_data=hover_data,
        hover_name=per.name,
        nbins=NBINS,
    )
    fig.update_xaxes(range=XRANGE)
    return fig


def plot_auc_scatter(df: DataFrame, per: config.Per) -> Figure:
    """
    Plot AUC scatter.

    Parameters
    ----------
    df
        DataFrame.
    per
        Per what?
    """
    if per == config.Per.organism:
        color: Optional[str] = "domain"
        hover_name = per.name
        hover_data = ["domain", config.label.auc, config.label.hmmer_hits]
        fields = ["organism", "domain"]
    else:
        color = None
        hover_name = per.name
        hover_data = [per.name, "clan", config.label.auc, config.label.hmmer_hits]
        fields = [per.name, "clan", config.label.auc, config.label.hmmer_hits]
        fields = list(set(fields))

    df = df.drop_duplicates(fields)
    df = df.sort_values(config.label.auc)
    title = f"AUC scatter, {per.name}-wise"
    fig = px.scatter(
        df,
        x=per.name,
        y=config.label.auc,
        color=color,
        title=title,
        hover_data=hover_data,
        hover_name=hover_name,
        marginal_y="violin",
    )
    fig.update_layout(showlegend=False)
    fig.update_yaxes(range=YRANGE)
    return fig


def plot_score_boxplot(dfe: DataFrame, per: config.Per) -> Figure:
    """
    Plot score bloxplot.

    Parameters
    ----------
    df
        DataFrame.
    per
        Per what?
    """
    precision = dfe.copy()
    precision["value"] = precision[config.label.precision]
    del precision[config.label.precision]
    precision["score"] = "precision"

    recall = dfe.copy()
    recall["value"] = recall[config.label.recall]
    del recall[config.label.recall]
    recall["score"] = "recall"

    f1score = dfe.copy()
    f1score["value"] = f1score[config.label.f1score]
    del f1score[config.label.f1score]
    f1score["score"] = config.label.f1score

    if per == config.Per.organism:
        fields = ["organism", "domain"]
        hover_data = [
            "organism",
            "domain",
            config.label.auc,
            config.label.hmmer_hits,
        ]
    else:
        hover_data = [
            per.name,
            "clan",
            config.label.auc,
            config.label.hmmer_hits,
        ]
        hover_data = list(set(hover_data))
        fields = [per.name, "clan", config.label.auc, config.label.hmmer_hits]
        fields = list(set(fields))

    dfe = pd.concat([precision, recall, f1score])

    title = f"Score boxplot, {per.name}-wise"
    fig = px.box(
        dfe,
        x="-log10(e-value)",
        color="score",
        y="value",
        title=title,
        hover_name=per.name,
        hover_data=hover_data,
    )
    fig.update_yaxes(range=YRANGE)

    return fig


def plot_pr_scatter(dfe: DataFrame, per: config.Per, size_max: int) -> Figure:
    """
    Plot Precision-Recall scatter.

    Parameters
    ----------
    df
        DataFrame.
    per
        Per what?
    """
    if per == config.Per.organism:
        color: Optional[str] = "domain"
        hover_data = [
            "organism",
            "domain",
            config.label.f1score,
            config.label.auc,
            config.label.hmmer_hits,
        ]
    else:
        color = None
        hover_data = [
            per.name,
            "clan",
            config.label.f1score,
            config.label.auc,
            config.label.hmmer_hits,
        ]
    title = f"Precision vs Recall, {per.name}-wise"
    fig = px.scatter(
        dfe,
        x=config.label.recall,
        y=config.label.precision,
        animation_frame="e-value",
        color=color,
        title=title,
        size=config.label.hmmer_hits,
        size_max=size_max,
        hover_data=hover_data,
        hover_name=per.name,
    )
    fig.update_layout(showlegend=False)
    fig.update_xaxes(range=XRANGE)
    fig.update_yaxes(range=YRANGE)
    fig.show()
