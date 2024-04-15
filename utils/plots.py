import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from typing import Optional
from rdkit import Chem
from rdkit.Chem import Draw
import bisect


def normal_heat(df, x, y, **kwargs):
    title = kwargs.get("title", "seriously, no title ?")
    fig = px.density_heatmap(
        df,
        x=x,
        y=y,
        #    nbinsx=20, nbinsy=20,
        color_continuous_scale="RdPu",
        title=title,
        # range_x=[-1,12],range_y=[-1,12],
        height=650,
        width=650,
    )  # text_auto=True
    fig.show()


def heat_log(df: pd.DataFrame, x: "str", y: str, **kwargs):
    title = kwargs.get("title", "seriously, no title ?")
    xaxis_title = kwargs.get("xaxis_title", "")
    yaxis_title = kwargs.get("yaxis_title", "")
    tickvals = kwargs.get(
        "tickvals",
        [
            0,
            20,
            40,
            60,
            80,
            100,
            120,
            140,
            160,
            180,
            200,
            220,
            240,
            260,
            280,
            300,
        ],
    )
    z = pd.crosstab(index=df[y], columns=df[x])
    fig = go.Figure(
        go.Heatmap(
            z=z,
            y=z.index,
            x=z.columns,
            colorscale=[
                [0, "rgb(250, 250, 250)"],  # 0
                [1.0 / 10000, "rgb(200, 200, 200)"],  # 10
                [1.0 / 1000, "rgb(150, 150, 150)"],  # 100
                [1.0 / 100, "rgb(100, 100, 100)"],  # 1000
                [1.0 / 10, "rgb(50, 50, 50)"],  # 10000
                [1.0, "rgb(0, 0, 0)"],  # 100000
            ],
            colorbar=dict(tick0=0, tickmode="array", tickvals=tickvals),
        )
    )
    fig.update_layout(
        autosize=False,
        width=650,
        height=650,
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
    )
    fig.show()


def nice_heat_log(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str = "",
    fig_width: int = 650,
    fig_height: int = 650,
    xaxis_title: str = "",
    yaxis_title: str = "",
    idx_max: int = 27,
    col_max: int = 27,
    base_colour: list = [0, 0, 0],
    logbase: int = 10,
    tickvals: list = [
        0,
        1,
        2,
        3,
        4,
        5,
        20,
        40,
        60,
        80,
        100,
        120,
        140,
        160,
        180,
        200,
        220,
        240,
        260,
        280,
        300,
    ],
    family_font: str = "Arial",
    xaxis_range: list = None,
    yaxis_range: list = None,
) -> go.Figure:
    # tickvals = tickvals[0]
    lw = 1

    z = pd.crosstab(index=df[y], columns=df[x])

    a = np.linspace(base_colour, 255, endpoint=True, num=6, dtype=None, axis=0)
    a = a.round().tolist()
    a_sumed = [
        f'rgb({", ".join(i)})' for i in [[str(int(i)) for i in j] for j in a]
    ]
    val_list = [1.0 / (logbase**i) for i in range(0, 5)] + [0.0]
    colorscale_value = [[i, j] for i, j in zip(val_list, a_sumed)]
    colorscale_value.reverse()

    # fill from 1
    if z.index.min() > 1:
        d = pd.DataFrame(
            0, index=list(range(1, int(z.index.min()))), columns=z.columns
        )
        z = pd.concat([d, z], axis=0)

    if z.columns.min() > 1:
        d = pd.DataFrame(
            0, columns=list(range(1, int(z.columns.min()))), index=z.index
        )
        z = pd.concat([d, z], axis=1)

    # first add row
    if z.index.max() < idx_max:
        d = pd.DataFrame(
            0,
            index=list(range(int(z.index.max()) + 1, idx_max + 1)),
            columns=list(z.columns),
        )
        z = pd.concat([z, d], axis=0)

    # add column
    if z.columns.max() < col_max:
        d = pd.DataFrame(
            0,
            columns=list(range(int(z.columns.max()) + 1, col_max + 1)),
            index=list(z.index),
        )
        z = pd.concat([z, d], axis=1)

    fig = go.Figure(
        go.Heatmap(
            z=z,
            y=z.index,
            x=z.columns,
            colorscale=colorscale_value,
            colorbar=dict(
                tick0=0,
                tickmode="array",
                tickvals=tickvals,
                borderwidth=0,
                title="count",
            ),
        )
    )
    fig.update_layout(
        autosize=False,
        width=fig_width,
        height=fig_height,
        title=dict(
            text=title, font=dict(family=family_font, size=18, color="#000000")
        ),
        font=dict(family=family_font, size=14, color="black"),
        yaxis=go.layout.YAxis(
            linecolor="black",
            linewidth=lw,
            mirror=True,
            ticklen=7,
            ticks="outside",
        ),
        xaxis=go.layout.XAxis(
            linecolor="black",
            linewidth=lw,
            mirror=True,
            ticklen=7,
            ticks="outside",
        ),
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        xaxis_range=xaxis_range,
        yaxis_range=yaxis_range,
    )

    return fig


def heat_log_with_binning(
    df: pd.DataFrame,
    x: str,
    y: str,
    bins: list,
    title: str = "",
    xaxis_title: str = "",
    yaxis_title: str = "",
) -> go.Figure:
    df["MolBin"] = df[x].apply(
        lambda i: bins[int(bisect.bisect_left(bins, i)) - 1]
    )
    z = pd.crosstab(index=df[y], columns=df["MolBin"])
    fig = go.Figure(
        go.Heatmap(
            z=z,
            y=z.index,
            x=z.columns,
            colorscale=[
                [0, "rgb(250, 250, 250)"],
                [1.0 / 10000, "rgb(200, 200, 200)"],
                [1.0 / 1000, "rgb(150, 150, 150)"],
                [1.0 / 100, "rgb(100, 100, 100)"],
                [1.0 / 10, "rgb(50, 50, 50)"],
                [1.0, "rgb(0, 0, 0)"],
            ],
            colorbar=dict(
                tick0=0,
                tickmode="array",
            ),
        ),
    )
    # fig.update_traces(customdata=pyLogo)

    fig.update_layout(
        autosize=False,
        width=650,
        height=650,
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
    )
    return fig


def show_molecules(
    df: pd.DataFrame,
    inchi_tag: str = "inchi",  # name of the inchi cotaining column
    assembly_index: int = "assembly_index",  # name of the inchi containing MA
    range_of_indexes: list[int] = None,
    how_many_structures: Optional[int] = 50,
):
    """Creates a grid plot of all molecules in the dataframe.
    It requires column of InChI. If it has assembly index and
    number of peaks, it can print it to each compound.

    Args:
        df (Pandas DataFrame): _description_
        inchi_tag (str, optional): _description_. Defaults to 'inchi'.
        assembly_index (str, optional): _description_. Defaults to 'assembly_index'.
        peaks_tag (str, optional): _description_. Defaults to 'filtered_peaks'.

    Returns:
        _type_: _description_
    """

    from rdkit.Chem import Draw

    ms = []
    labels = []
    if range_of_indexes:
        df2plot = df.iloc[range_of_indexes[0] : range_of_indexes[1], ::]
    else:
        if len(df) > how_many_structures:
            df2plot = df.sample(n=how_many_structures, random_state=100)

    for i in df2plot.iterrows():
        inchi = df.loc[i[0], inchi_tag]
        try:
            ma = int(df.loc[i[0], assembly_index])
        except:
            ma = "n.a."

        m = Chem.MolFromInchi(inchi)
        ms.append(m)
        labels.append(f"MA = {ma}")

    return Draw.MolsToGridImage(
        ms, molsPerRow=5, legends=labels, subImgSize=(275, 180), maxMols=100
    )
