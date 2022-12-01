import logging
import subprocess
import shlex
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from . import setup
from fermipy.gtanalysis import GTAnalysis
from fermipy.utils import merge_dict
from gammapy.maps import WcsNDMap
from astropy.table import Table
from astropy import units
from astropy.utils import lazyproperty
from gammapy.maps import MapAxis
from scipy.interpolate import UnivariateSpline, RectBivariateSpline
from scipy.optimize import curve_fit
from .plotting import effect_k, plot_r68
import fermiAnalysis


def generate_psmap(gta, state, overwrite=False, prob_epsilon=1e-7, make_plots=False):
    """
    Generate the model maps and ps maps
    """

    # export the model map for each component
    gta.write_model_map(state)

    cen_src = gta.get_sources()[0]

    script = os.path.join(os.path.dirname(fermiAnalysis.__file__), 'scripts/', 'gtpsmap.py')
    
    # run gtpsmap for each component
    ccubes = glob.glob(os.path.join(gta.config['fileio']['workdir'], "ccube*"))
    mcubes = glob.glob(os.path.join(gta.config['fileio']['workdir'], "mcube*"))
    
    ccubes = sorted(ccubes, key=lambda f: os.path.basename(f))
    mcubes = sorted(mcubes, key=lambda f: os.path.basename(f))
    psmaps = []

    cp = cm
    script = 
    for i, c in enumerate(ccubes):

        logging.info("Running gtpsmap for {0:s} and {1:s}".format(os.path.basename(c), os.path.basename(mcubes[i])))

        # check if we're dealing with a component
        component = os.path.basename(c).rstrip(".fits").split("_")
        if len(component) > 1:
            out_suffix = "_" + component[1] + ".fits"
        else:
            out_suffix = ".fits"

        # check if psf file exists
        if len(component) > 1:
            psf_file = os.path.join(gta.config['fileio']['workdir'], "psf" + out_suffix)
        else:
            # use zeroth component
            psf_file = os.path.join(gta.config['fileio']['workdir'], "psf_00.fits")

        if not os.path.exists(psf_file):
            gta.compute_psf()

        logging.info("using PSF {0:s}".format(os.path.basename(psf_file)))
        # read in the psf
        psf = PSFReader(psf_file)
        # get the 68% confidence radius in deg
        r68 = psf.containment_radius(gta.energies[0], fraction=0.68)
        logging.info("r68 at {0:.2f} MeV: {1:.2f} deg".format(gta.energies[0], r68))
        
        # perfrom fit to containent radius
        psf_params = psf.fit_containment(fraction=0.68)
        logging.info("PSF params at 100 MeV: {0}".format(psf_params))

        outfile = os.path.join(gta.config['fileio']['workdir'], state + "_" + "psmap" + out_suffix)
        logging.info("outfile: {0:s}".format(os.path.basename(outfile)))

        if not os.path.exists(outfile) or overwrite:
            cmd = "python {4:s} --cmap {0:s} --mmap {1:s} --outfile {2:s} --prob_epsilon {3:.3e}".format(ccubes[0],
                mcubes[0], outfile, prob_epsilon, script)

            cmd += " --psfpar0 {0[0]:.5f} --psfpar1 100.  --psfpar2 {0[1]:.5f} --psfpar3 {0[2]:.5f}".format(psf_params)
            logging.info("running command: {0:s}".format(cmd))
            subprocess.call(shlex.split(cmd))
        else:
            logging.info("{0:s} exists, skipping gtpsmap".format(outfile))
        psmaps.append(outfile)

        if make_plots:
            ps_map = WcsNDMap.read(outfile)
            p = MyROIPlotter(ps_map, roi=gta.roi)
            p.config['cmap'] = 'coolwarm'
            ax, cb = p.myplot(vmin=-4, vmax=4, extend='both', interpolation='bicubic',zoom=2, cb_label="$\log_{10}(\mathrm{PS})$", cmap='coolwarm')

            ax.coords.grid(color='k', ls=":", lw=0.5)
            ax.coords['ra'].set_ticks(color="k")
            ax.coords['dec'].set_ticks(color="k")
            ax.coords['ra'].ticks.set_tick_out(True)
            ax.coords['dec'].ticks.set_tick_out(True)
            ax.set_xlabel("Right Ascension")
            ax.set_ylabel("Declination")

            # add PSF circle
            plot_r68(gta, ax, psf_file, edgecolor='k', lw=2, ls='-', facecolor='none')

            ax.annotate("{0:s}".format(cen_src.assoc["ASSOC_TEV"]), xy=(0.05,0.95), va='top', color='w', xycoords="axes fraction", **effect_k)

            plt.subplots_adjust(left=0.15, bottom=0.12, top=0.95, right=0.9)
            for plot_format in ['png', 'pdf']:
                plt.savefig(outfile.replace(".fits", "_{1:s}.{0:s}".format(plot_format, cen_src.assoc["ASSOC_TEV"].replace(" ",""))), dpi=120)
            plt.close("all")

    return psmaps


def psf_fit(E, *p):
    """The PSF fitting function for gtpsmap for E=100 MeV"""
    return np.sqrt(p[0]**2. * np.power(E / 100., -2. * p[1]) + p[2]**2.)

p0 = [4.,0.9,0.1]  # the initial parameters

class PSFReader(object):
    """
    Class to read in Fermi PSF and calculate containment radius
    inspired from gammapy EnergyDependentIRFTable

    Written for gammapy 0.9...
    """
    def __init__(self, psf_fits_file):
        """
        Init the class

        Parameters
        ----------
        psf_fits_file: str
            path the psf fits file generated with fermipy
        """
        psf_table = Table.read(psf_fits_file, hdu="PSF")
        theta = Table.read("psf_01.fits", hdu="THETA")['Theta'].data

        self._rad_axis = MapAxis.from_nodes(theta)
        self._energy = psf_table['Energy'].data
        self._log_energy = np.log10(psf_table['Energy'].data)
        self._psf = psf_table['Psf'].data

    @property
    def energy(self):
        return self._energy

    @lazyproperty
    def _interpolate_containment(self):
        """Perform a 2D interpolation of CDF of integrated PSF over energy and theta angle"""
        psf_int = 2. * np.pi * self._rad_axis.center * self._psf 
        psf_int *= (self._rad_axis.edges[1:] - self._rad_axis.edges[:-1])

        cdf = np.cumsum(psf_int, axis=1)  
        cdf = ((cdf.T - cdf[:,0]) / (cdf[:,-1] - cdf[:,0])).T

        return RectBivariateSpline(self._log_energy, self._rad_axis.center, cdf)

    def containment(self, theta):
        """Get the cointainment angle from 2D interpolation"""
        theta = np.atleast_1d(theta)
        return self._interpolate_containment(self._log_energy, theta)

    def containment_radius(self, energy, fraction=0.68):
        """
        Calculate containment radius as a function of energy

        Parameters
        ----------
        energy: array-like
            energies to calculate containment radius in MeV

        fraction: float
            Fraction of PSF containment (default: 0.68)

        Returns 
        -------
        Containment radii for each energy in deg
        """

        theta_upsample = np.linspace(0., self._rad_axis.center[-1], 10 * self._rad_axis.nbin)
        containment = self.containment(theta=theta_upsample)

        fraction_idx = np.argmin(np.abs(containment - fraction), axis=1)
        
        spline = UnivariateSpline(self._log_energy, theta_upsample[fraction_idx], k=1, s=0)

        return spline(np.log10(energy))

    def fit_containment(self, fraction=0.68):
        """
        Fit the containment radius as function of energy

        Parameters
        ----------
        fraction: float
            Fraction of PSF containment (default: 0.68)

        Returns 
        -------
        best-fit parameters
        """
        containment_radius = self.containment_radius(self._energy, fraction = fraction)

        # perform the fit 
        p_final, _ = curve_fit(psf_fit, self.energy, containment_radius, p0=p0)

        return p_final
