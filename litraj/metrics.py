import numpy as np
from scipy.stats import linregress
from ase import geometry


def mae(y, y_hat):
    """mean absolute error"""
    return abs(np.array(y) - y_hat).mean()



def rmse(y, y_hat):
    """root mean squared error"""
    return np.sqrt(np.mean((np.array(y)-y_hat)**2))



def r2_score(y, y_hat):
    """
    coefficient of determination
    """
    tss = np.sum((y - np.mean(y)) ** 2)
    rss = np.sum((y - y_hat) ** 2)
    r2 = 1 - (rss / tss)
    return r2



def get_metrics(y, y_hat):

    """
    calculate MAE, RMSE, Rp, slope, and R2 metrics
    """

    res = linregress(y, y_hat)
    metrics = {
        'mae': mae(y, y_hat),
        'rmse': rmse(y, y_hat),
        'Rp': res.rvalue,
        'slope': res.slope,
        'R2': r2_score(y, y_hat)
    }
    return metrics



def PA_MAE(trajectory_true, trajectory_surrogate):
    """
    Path averaged mean absolute error of geometry prediction
    """
    cumulative_path_displacement = 0
    for im_true, im_surrogate in zip(trajectory_true, trajectory_surrogate):
        mean_displacement = 0
        for p1, p2 in zip(im_true.positions, im_surrogate.positions):
            mean_displacement += geometry.get_distances(p1, p2, cell = im_true.cell, pbc = True)[1][0][0]
        mean_displacement /= len(im_true)
        cumulative_path_displacement += mean_displacement
    return cumulative_path_displacement / len(trajectory_true)
    


