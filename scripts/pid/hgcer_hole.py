import pandas as pd
import scipy
import numpy as np


def define_hgcer_hole_boundary(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    x0: float = 0.0,
    y0: float = 0.0,
    thres:float =0.1,
    max_dr:float=10,
    n_angles:int=360,
    smooth_sigma:float=2,
):
    """
    Given a 2D histogram defined by x, y coordinates and content values, determine the low-threshold boundary. The algorithm searches radially from the specified center and identifies the boundary where the content exceeds the threshold.

    Parameters:
    -----------
    x : np.ndarray
        x coordinates of the given histogram
    y : np.ndarray
        y coordinates of the given histogram
    z : np.ndarray
        contents of the histogram (efficiency for determining hgcer hole)
    x0 : float
        x coordinate of the hole center
    y0 : float
        y coordinate of the hole center
    thres : float
        Threshold value to determine hole boundary
    max_dr : float
        Maximum radial distance to search
    n_angles : int
        Number of angular divisions for polar coordinate search
    smooth_sigma : float
        Sigma for Gaussian smoothing of the boundary points

    Returns:
    --------
    Tuple of np.ndarray: arrays containing x and y coordinates of boundary points
    """

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # Polar coordinate search from center
    theta_arr = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    r_pts = []
    x_pts = []
    y_pts = []

    for theta in theta_arr:
        r_boundary = 0.0

        # Radial search with fine step size
        r = 0.0
        dr = min(dx, dy) / 2  # Fine resolution for accurate boundary

        while r <= max_dr:
            x_pt = x0 + r * np.cos(theta)
            y_pt = y0 + r * np.sin(theta)

            # Find bin indices using searchsorted for accuracy
            i = np.searchsorted(x, x_pt) - 1
            j = np.searchsorted(y, y_pt) - 1

            # Check bounds
            if 0 <= i < len(x) - 1 and 0 <= j < len(y) - 1:
                # Use bilinear interpolation for smoother boundary
                z_interp = bilinear_interp(x, y, z, x_pt, y_pt)

                if z_interp > thres:
                    r_boundary = r
                    break
            else:
                # Out of histogram bounds
                break

            r += dr

        # Store boundary point
        if r_boundary > 0:
            x_pt = x0 + r_boundary * np.cos(theta)
            y_pt = y0 + r_boundary * np.sin(theta)
            r_pts.append(r_boundary)
            x_pts.append(x_pt)
            y_pts.append(y_pt)

    if len(x_pts) == 0:
        print("Warning: No boundary points found")
        return np.array([]), np.array([])

    # Convert to arrays
    r_pts = np.array(r_pts)
    x_pts = np.array(x_pts)
    y_pts = np.array(y_pts)

    # Smooth the boundary using Gaussian filter on radial coordinates
    if smooth_sigma > 0 and len(r_pts) > 10:
        # Smooth in polar coordinates to preserve shape
        r_smooth = scipy.ndimage.gaussian_filter1d(
            r_pts, sigma=smooth_sigma, mode='wrap'
        )

        theta_smooth = np.linspace(0, 2 * np.pi, len(r_smooth), endpoint=False)
        x_pts = x0 + r_smooth * np.cos(theta_smooth)
        y_pts = y0 + r_smooth * np.sin(theta_smooth)

    return np.array(x_pts), np.array(y_pts)


def bilinear_interp(x, y, z, x_pt, y_pt):
    """
    Perform bilinear interpolation on 2D histogram.

    Parameters:
    -----------
    x, y : array
        Bin edges
    z : 2D array
        Histogram values
    x_pt, y_pt : float
        Point to interpolate

    Returns:
    --------
    float : interpolated value
    """
    # Find bin indices
    i = np.searchsorted(x, x_pt) - 1
    j = np.searchsorted(y, y_pt) - 1

    # Boundary check
    if i < 0 or i >= len(x) - 1 or j < 0 or j >= len(y) - 1:
        return 0.0

    # Get bin edges for interpolation
    x1, x2 = x[i], x[i + 1]
    y1, y2 = y[j], y[j + 1]

    # Normalized coordinates within bin
    wx = (x_pt - x1) / (x2 - x1) if x2 != x1 else 0.5
    wy = (y_pt - y1) / (y2 - y1) if y2 != y1 else 0.5
    wx = np.clip(wx, 0, 1)
    wy = np.clip(wy, 0, 1)

    # Ensure indices are within bounds
    i_max = z.shape[0] - 1
    j_max = z.shape[1] - 1

    # Get corner values
    z00 = z[min(i, i_max), min(j, j_max)]
    z10 = z[min(i + 1, i_max), min(j, j_max)]
    z01 = z[min(i, i_max), min(j + 1, j_max)]
    z11 = z[min(i + 1, i_max), min(j + 1, j_max)]

    # Bilinear interpolation
    z_interp = (
        (1 - wx) * (1 - wy) * z00
        + wx * (1 - wy) * z10
        + (1 - wx) * wy * z01
        + wx * wy * z11
    )

    return z_interp
