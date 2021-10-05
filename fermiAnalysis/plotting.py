import numpy as np
import matplotlib.pyplot as plt
from .utils import interp_2d_likelihood

def plot_model(model_flux, **kwargs):

    ax = kwargs.pop('ax', plt.gca())

    color = kwargs.pop('color', 'k')
    noband = kwargs.pop('noband', False)
    nomargin = kwargs.pop('nomargin', True)

    e2 = 10 ** (2 * model_flux['log_energies'])

    ax.plot(10 ** model_flux['log_energies'],
        model_flux['dnde'] * e2,
        color=color, zorder = -10,
        label = '$Fermi$-LAT best fit')

    if not nomargin:
        ax.plot(10 ** model_flux['log_energies'],
            model_flux['dnde_lo'] * e2,
            color=color,
            linestyle='--', zorder = -10)
        ax.plot(10 ** model_flux['log_energies'],
            model_flux['dnde_hi'] * e2,
            color=color,
            linestyle='--', zorder = -10)

    if not noband:
        ax.fill_between(10 ** model_flux['log_energies'],
            model_flux['dnde_lo'] * e2,
            model_flux['dnde_hi'] * e2,
            alpha=0.5, color=color, zorder=-10)

def plot_data(sed, **kwargs):
    # FERMI
    ax = kwargs.pop('ax', plt.gca())
    m = (sed['ts'] > 4.) 

    plt.errorbar(sed['e_ref'][m], sed['e2dnde'][m], 
        yerr = np.array([sed['e2dnde_err_lo'][m],sed['e2dnde_err_hi'][m]]),
        xerr = np.array([sed['e_ref'] - sed['e_min'],sed['e_max'] - sed['e_ref']])[:,m],
        marker = 's', ls = 'None', color = '0.5', 
        label = "$Fermi$ LAT")

    plt.errorbar(sed['e_ref'][~m], sed['e2dnde_ul'][~m], 
        yerr = np.array([sed['e2dnde_ul'][~m] / 3., np.zeros(np.sum(~m))]),
        xerr = np.array([sed['e_ref'] - sed['e_min'],sed['e_max'] - sed['e_ref']])[:,~m],
        uplims = True,
        marker = 's', ls = 'None', color = '0.5')

                    
    plot_model(sed['model_flux'], color = "0.5")

    plt.xlabel("Energy [MeV]")
    plt.ylabel("$E^2dN/dE\,[\mathrm{MeV}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$")
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.legend()

def plot_2d_delta_likelihood(x, y, logl, cbar=True, fig=None, ax=None, log_x=True, log_y=True, intp_kwargs={}, **kwargs):
    """Plot a 2d likelihood surface"""
    cp = plt.cm.get_cmap(kwargs.pop('cmap', 'PuBuGn_r'))

    if fig is None:
        fig = plt.figure(figsize=(6,4))
    if ax is None:
        ax = fig.add_subplot(111)

    xx, yy, ll = interp_2d_likelihood(x, y, logl, log_x=log_x, log_y=log_y, **intp_kwargs)

    im = ax.pcolormesh(xx, yy, ll - ll.max(), cmap=cp, **kwargs)
    if cbar:
        plt.colorbar(im, label='$\Delta\ln\mathcal{L}$')

    if log_x:
        ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')

    return fig, ax, im
