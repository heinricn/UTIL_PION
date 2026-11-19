import json
import pathlib
import numpy as np
from collections import OrderedDict
from typing import Union


def load_json(fname: Union[str, pathlib.Path]) -> OrderedDict:
    """
    Load JSON file and return as an OrderedDict.
    fname : str | pathlib.Path
        Path to the JSON file.
    Returns
    -------
    OrderedDict
        Contents of the JSON file as an OrderedDict.
    """
    fname = pathlib.Path(fname)
    with fname.open("rt") as handle:
        return json.load(handle, object_hook=OrderedDict)


def point_in_polygon(
    x: np.ndarray,
    y: np.ndarray,
    polygon: np.ndarray,
) -> np.ndarray:
    """

    Parameters
    ----------
    x : np.ndarray
        Array of x coordinates of the points to test.
    y : np.ndarray
        Array of y coordinates of the points to test.
    polygon : np.ndarray
        Array of shape (N, 2) representing the vertices of the polygon.

    Returns
    -------
    np.ndarray
        Boolean array indicating whether each point is inside the polygon.

    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    polygon = np.asarray(polygon, dtype=float)

    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape.")

    if polygon.ndim != 2 or polygon.shape[1] != 2:
        raise ValueError("polygon must have shape (N, 2).")
    if polygon.shape[0] < 3:
        raise ValueError("polygon must contain at least three points.")

    inside = np.zeros_like(x, dtype=bool)
    x1, y1 = polygon[-1]
    for x2, y2 in polygon:
        cond = (y1 > y) != (y2 > y)
        denom = y2 - y1
        denom = np.where(denom == 0.0, np.nan, denom)
        xinters = (x2 - x1) * (y - y1) / denom + x1
        toggle = cond & (x < xinters)
        inside ^= np.nan_to_num(toggle, nan=False)
        x1, y1 = x2, y2
    return inside


def load_data(path: Union[str, pathlib.Path]):
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:,0], data[:,1]

# example use

def main():

    hgcer_x, hgcer_y = load_data("P_hgcer.csv")
    points = load_json("Q4p4_W2p74_left_lowe/03.json")

    # e.g. [True, True, True, False]
    is_inside = point_in_polygon(
        x=np.array(hgcer_x),
        y=np.array(hgcer_y),
        polygon=np.array(points),
    )
    make_plot(hgcer_x, hgcer_y, path="no_cut.png")
    make_plot(hgcer_x[~is_inside], hgcer_y[~is_inside], path="cut.png")


def make_plot(x, y, path):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "figure.dpi": 120,
        "figure.figsize": (4, 3.5),
        "figure.facecolor": "white",
        "xtick.top": True,
        "xtick.direction": "in",
        "xtick.minor.visible": True,
        "ytick.right": True,
        "ytick.direction": "in",
        "ytick.minor.visible": True,
    })

    fig, ax = plt.subplots(figsize=(5, 4), dpi=300, constrained_layout=True)

    hist, xedges, yedges = np.histogram2d(
        x,
        y,
        bins=[120,120],
        range=[[-60.0, 60.0], [-60.0, 60.0]],
    )
    x = 0.5 * (xedges[:-1] + xedges[1:])
    y = 0.5 * (yedges[:-1] + yedges[1:])
    ax.pcolormesh(x, y, hist.T, cmap="rainbow", norm=mpl.colors.LogNorm())
    ax.set_xlabel(r"$\mathrm{SHMS\ HGCER\ x [cm]}$")
    ax.set_ylabel(r"$\mathrm{SHMS\ HGCER\ y [cm]}$")
    fig.savefig(path)
    plt.close(fig)



if __name__ == "__main__":
    main()