# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
##############################################################################
Production envelope module
Author: Vitor Pereira
##############################################################################
"""
import numpy as np
from mewpy.simulation import get_simulator


def flux_envelope(model, r_x, r_y, steps=10, constraints=None, x_range=None, tolerance=0):
    """ Calculate the flux envelope for a pair of reactions.
        Adapted from REFRAMED to be compatible both with REFRAMED and COBRApy.


    :param model : The model or simulator.
    :param str r_x: Reaction on x-axis.
    :param str r_y: Reaction on y-axis.
    :param int steps: Number of steps to compute (default: 10).
    :param dict constraints: Custom constraints to the FBA problem.
    :param dict envcond: Environmental conditions.
    :param tuple range: x value range. Default None.
    :returns:  x values, y_min values, y_max values

    """

    try:
        simul = get_simulator(model)
    except Exception:
        raise ValueError(
            'The model should be an instance of model or simulator')

    obj_frac = 0
    # if r_x in simul.get_objective():
    #    obj_frac = 0.0

    xrange = simul.FVA(obj_frac=obj_frac, reactions=[r_x], constraints=constraints)
    xmin, xmax = xrange[r_x]
    xmax -= tolerance
    if x_range:
        if x_range[0] > xmin:
            xmin = x_range[0]
        if x_range[1] < xmax:
            xmax = x_range[1]
    xvals = np.linspace(xmin, xmax, steps)
    ymins, ymaxs = np.zeros(steps), np.zeros(steps)

    if constraints is None:
        _constraints = {}
    else:
        _constraints = {}
        _constraints.update(constraints)

    for i, xval in enumerate(xvals):
        _constraints[r_x] = xval
        y_range = simul.FVA(obj_frac=obj_frac, reactions=[r_y], constraints=_constraints)
        ymins[i], ymaxs[i] = y_range[r_y]

    return xvals, ymins, ymaxs


def plot_flux_envelope(model, r_x, r_y, steps=10, substrate=None, constraints=None,
                       label_x=None, label_y=None, flip_x=False, flip_y=False,
                       plot_kwargs=None, fill_kwargs=None, ax=None, x_range=None):
    """ Plots the flux envelope for a pair of reactions.
        Adapted from REFRAMED.

    :param model: The model or simulator.
    :param str r_x: Reaction on x-axis.
    :param str r_y: Reaction on y-axis.
    :param int steps: Number of steps to compute (default: 20).
    :param str substrate: Compute yields for given substrate instead of rates (optional).
    :param dict constraints: Additional simulation constraints.
    :param str label_x: x label (optional, uses reaction name by default).
    :param str label_y: y label (optional, uses reaction name by default).
    :param bool flip_x: Flip direction of r_x (default: False).
    :param dict flip_y: Flip direction of r_y (default: False).
    :param dict plot_kwargs: Additional parameters to *pyplot.plot* (optional).
    :param dict fill_kwargs: Additional parameters to *pyplot.fill_between* (optional).
    :param matplotlib.Axes ax: Plot over existing axes (optional).
    :param tuple range: x value range. Default None.
    :returns:  matplotlib.Axes: Axes object.
    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise RuntimeError("Matplotlib is not installed.")

    try:
        simul = get_simulator(model)
    except Exception:
        raise ValueError(
            'model should be an instance of model or simulator')

    offset = 0.03

    if ax is None:
        _, ax = plt.subplots()

    if not plot_kwargs:
        plot_kwargs = {'color': 'k'}

    if not fill_kwargs:
        fill_kwargs = {'color': 'k', 'alpha': 0.1}

    xvals, ymins, ymaxs = flux_envelope(model, r_x, r_y, steps, constraints, x_range=x_range)

    if flip_x:
        xvals, ymins, ymaxs = -xvals, ymins[::-1], ymaxs[::-1]

    if flip_y:
        ymins, ymaxs = -ymaxs, -ymins

    if substrate:
        sol = simul.simulate()
        uptk = abs(sol.fluxes[substrate])
        xvals, ymins, ymaxs = xvals / uptk, ymins / uptk, ymaxs / uptk

    ax.plot(xvals, ymins, **plot_kwargs)
    ax.plot(xvals, ymaxs, **plot_kwargs)
    ax.plot([xvals[0], xvals[0]], [ymins[0], ymaxs[0]], **plot_kwargs)
    ax.plot([xvals[-1], xvals[-1]], [ymins[-1], ymaxs[-1]], **plot_kwargs)

    ax.fill_between(xvals, ymins, ymaxs, **fill_kwargs)

    ax.set_xlabel(label_x) if label_x else ax.set_xlabel(r_x)
    ax.set_ylabel(label_y) if label_y else ax.set_ylabel(r_y)

    xmin, xmax = min(xvals), max(xvals)
    dx = offset * (xmax - xmin)
    ax.set_xlim((xmin - dx, xmax + dx))

    ymin, ymax = min(ymins), max(ymaxs)
    dy = offset * (ymax - ymin)
    ax.set_ylim((ymin - dy, ymax + dy))

    return ax
