import numpy as np
import matplotlib.pyplot as plt
from .utils import interp_2d_likelihood
from fermipy.plotting import ROIPlotter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patheffects import withStroke
from astropy.visualization import wcsaxes

# Create my own colormap
colors = [
            (0, 0, 0),
            plt.cm.tab10(0.),
            plt.cm.tab10(0.1),
            (1, 1, 1)
         ] 
cmap_name = "my_cmap"
cm = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=200)
plt.register_cmap(cmap=cm)
effect = dict(path_effects=[withStroke(foreground="w", linewidth=2)])
effect_k = dict(path_effects=[withStroke(foreground="k", linewidth=2)])

class MyROIPlotter(ROIPlotter):
    def __init__(self, *args, **kwargs):
        super(MyROIPlotter, self).__init__(*args, **kwargs)

    def myplot(self, **kwargs):

        zoom = kwargs.get('zoom', None)
        graticule_radii = kwargs.get('graticule_radii',
                                     self.config['graticule_radii'])
        label_ts_threshold = kwargs.get('label_ts_threshold',
                                        self.config['label_ts_threshold'])

        im_kwargs = dict(cmap=self.config['cmap'],
                         interpolation='nearest', transform=None,
                         vmin=None, vmax=None, clip=False, levels=None,
                         zscale='lin', subplot=111, colors=['k'])

        cb_kwargs = dict(orientation='vertical', shrink=1.0, pad=0.1,
                         extend='max',
                         fraction=0.1, cb_label=None)

        im_kwargs = merge_dict(im_kwargs, kwargs)
        cb_kwargs = merge_dict(cb_kwargs, kwargs)

        im, ax = self._implot.plot(**im_kwargs)

        self._ax = ax
        for c in self._catalogs:
            self.plot_catalog(c)

        if self._roi is not None:
            self.plot_roi(self._roi,
                          label_ts_threshold=label_ts_threshold)

        self._extent = im.get_extent()
        ax.set_xlim(self._extent[0], self._extent[1])
        ax.set_ylim(self._extent[2], self._extent[3])

        self.zoom(zoom)

        cb_label = cb_kwargs.pop('cb_label', None)
        cb = plt.colorbar(im, **cb_kwargs)

        if cb_label:
            cb.set_label(cb_label)
        for r in graticule_radii:
            self.draw_circle(r)
        return ax, cb

def plot_r68(gta, ax, psf_file, **kwargs):
    """
    add the r68 PSF containment radius to a sky map
    """
    kwargs.setdefault("transform", ax.get_transform('fk5'))
    kwargs.setdefault("facecolor", "none")
    kwargs.setdefault("edgecolor", "k")
    kwargs.setdefault("lw", 2)
    kwargs.pop("radius", None)

    # get the PSF 68% containment radius, use the last component
    logging.info("using PSF {0:s}".format(os.path.basename(psf_file)))
    psf = PSFReader(psf_file)
    r68 = psf.containment_radius(gta.energies[0], fraction=0.68)
    logging.info("r68 at {0:.2f} MeV: {1:.2f}".format(gta.energies[0], r68))

    c_psf = wcsaxes.SphericalCircle([gta.roi.skydir.ra, gta.roi.skydir.dec],
                                   radius=r68 * units.deg,
                                   **kwargs)
    ax.add_patch(c_psf)

def tsmap_plot_nice(gta, state):
    """make a nice ts map plot"""

    tsmap_file = os.path.join(gta.config['fileio']['workdir'], state + "_pointsource_powerlaw_2.00_tsmap.npy")
    tsmap = np.load(tsmap_file, allow_pickle=True).flat[0]
    sigma_levels = [3, 5, 7] + list(np.logspace(1, 3, 17))

    p = MyROIPlotter(tsmap['sqrt_ts'], roi=gta.roi)
    p.config['cmap'] = 'my_cmap'
    ax, cb = p.myplot(vmin=0, vmax=4, levels=sigma_levels,
                      cb_label='$\sqrt{\mathrm{TS}}$', interpolation='bicubic', cmap='my_cmap',
                      extend='max',
                      colors=['w'],
                      zoom=2)
    ax.coords.grid(color='w', ls=":", lw=0.5)
    ax.coords['ra'].set_ticks(color="k")
    ax.coords['dec'].set_ticks(color="k")
    ax.coords['ra'].ticks.set_tick_out(True)
    ax.coords['dec'].ticks.set_tick_out(True)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")

    psf_file = os.path.join(gta.config['fileio']['workdir'], "psf_01.fits")
    plot_r68(gta, ax, psf_file, lw=2, ls='-', edgecolor='w')

    cen_src = gta.get_sources()[0]
    #ax.annotate("$\mathbf{{{0:s}}}$".format(cen_src.assoc["ASSOC_TEV"]), xy=(0.05,0.95), va='top', color='k', xycoords="axes fraction", **effect)
    ax.annotate("{0:s}".format(cen_src.assoc["ASSOC_TEV"]), xy=(0.05,0.95), va='top', color='k', xycoords="axes fraction", **effect)

    plt.subplots_adjust(left=0.15, bottom=0.1, top=0.95, right=0.9)
    for plot_format in ['png', 'pdf']:
        plt.savefig(tsmap_file.replace(".npy", "_new_{1:s}.{0:s}".format(plot_format, cen_src.assoc["ASSOC_TEV"].replace(" ",""))), dpi=120)
    plt.close("all")


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
