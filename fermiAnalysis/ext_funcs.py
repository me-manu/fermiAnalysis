"""
Functions for extended source analysis

Adapted from Matthew Wood's extpipe package

see https://github.com/fermiPy/extpipe/tree/master/extpipe
and in particular https://github.com/fermiPy/extpipe/blob/master/extpipe/fit_funcs.py
"""

import os
import glob
import yaml
import sys
import copy
import itertools
import logging
from astropy.io import fits
from astropy import units as u
from scipy.interpolate import RectBivariateSpline as RBSpline
from scipy.interpolate import UnivariateSpline as USpline

import numpy as np
from astropy.table import Column, Table
from fermipy.utils import get_parameter_limits
from fermipy.gtanalysis import GTAnalysis
from fermiAnalysis.utils import set_src_spec_plexpcut, set_src_spec_pl, add_ebl_atten
import fermiAnalysis as fa
#from LikelihoodState import LikelihoodState
from collections import OrderedDict
from itertools import repeat
from multiprocessing import Pool

halo_source_dict_default = {
    'SpectrumType' : 'PowerLaw', 
    'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 1.0, 'max' : 4.5 }, 
    'Scale' : 1000,
    'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13, 'min' : 1E-5, 'max' : 1E4 }
    }

def stack_files(files, outfile, new_cols=None):

    h = fits.open(files[0])

    tables = []
    for hdu in h:
        if isinstance(hdu,fits.BinTableHDU):
            tables += [stack_tables(files,hdu.name,new_cols=new_cols)]

    hdus = [fits.PrimaryHDU()]
    hdus += [fits.table_to_hdu(t) for t in tables]
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfile,overwrite=True)
    

def stack_tables(files, hdu=None, new_cols=None):
    
    tables = []
    for f in sorted(files):
        tables += [Table.read(f,hdu=hdu)]

    cols = []
    for c in tables[0].colnames:

        col = tables[0][c]
        cols += [Column(name=col.name, unit=col.unit, shape=col.shape,
                        dtype=col.dtype)]

    tab = Table(cols,meta=tables[0].meta)

    for t in tables:
        row = [ t[c] for c in tables[0].colnames ]
        tab.add_row(row)

    if new_cols is not None:

        for col in new_cols:
            tab.add_column(col)
                    
    return tab


def fit_region(gta, modelname, src_name, loge_bounds=None,
               skip_opt=[], shape_ts_threshold=9.,
               force_ps=False,
               create_maps=True,
               create_sed=True,
               free_radius_sed=1.0,
               crate_ts_maps=True,
               distance_free_norm=1.5,
               distance_free_shape=1.0,
               **kwargs):
    """
    Fit an ROI and determine best fit of extended models

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        the gt analysis object

    modelname: str
        a name for the model

    src_name: str
        name of the source for which extension will be tested

    loge_bounds array-like or None
        if not None, use these energy bounds for the analysis. 
        Otherwise, use the energy bounds of the gta object

    {options}

    create_maps: bool
        compute ts and residual maps for different test sources
        prior to extension fitting.
        Default: True

    create_sed: bool
        compute SEDs
        prior to extension fitting.
        Default: True

    skip_opt: list
        list of source names that will be skipped in gta.optimize
       
    shape_ts_threshold: float
        above this value, shape parameters are left free in 
        gta.optimize. 
        Default = 9.

    force_ps: bool
        If true, force the source of interest to be modeled 
        as a point source for the base model

    free_radius_sed: float
        radius in degree in which normalizations are left free
        for SED.
        Default: 1.

    kwargs: dict
        additional key word arguments passed to gta.extension
    """
    kwargs.setdefault('fit_ebin', False)
    kwargs.setdefault('free_radius', 1.0)
    kwargs.setdefault('make_tsmap', False)
    kwargs.setdefault('fit_position', True)
    kwargs.setdefault('free_background', True)
    kwargs.setdefault('make_plots', True)
    kwargs.setdefault('update', True)

    gta.logger.info('Starting Region Fit {0:s}'.format(modelname))
    lnl0 = -gta.like()    
    gta.logger.info('{0:s} Model Likelihood: {1:f}'.format(modelname,lnl0))
    
    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0],loge_bounds[1])
    
    # point source models for TS maps
    model_pl20 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model_pl27 = { 'SpatialModel' : 'PointSource', 'Index' : 2.7 }

    # model for residual map
    model3 = { 'SpatialModel' : 'Gaussian', 'Index' : 2.0, 'SpatialWidth' : 0.1 }
    #model4 = { 'SpatialModel' : 'RadialDisk', 'Index' : 2.0,
               #'SpatialWidth' : 0.1 * 0.8246211251235321 }

    # force the central point source to be a point source
    if force_ps and not gta.roi[src_name]['SpatialModel'] == 'PointSource':
        gta.roi[src_name]['SpatialModel'] = 'PointSource'

    gta.print_roi()
    gta.print_params()
    
    # optimization. Optimize shape parameters with ts > ts_threshold (9 is default)
    gta.optimize(skip=skip_opt, shape_ts_threshold=shape_ts_threshold)

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]

    skydir = gta.roi[src_name].skydir

    # freeze all sources
    gta.free_sources(False)

    # free normalization within 1.5 deg
    gta.free_sources(skydir=skydir, distance=distance_free_norm, pars='norm')

    # free shape parameters within 1 deg except for diff sources
    gta.free_sources(skydir=skydir, distance=distance_free_shape, pars='shape', exclude=diff_sources)

    # free central source
    gta.free_source(src_name)
    gta.fit()
    gta.update_source(src_name,reoptimize=True)
    gta.write_roi(modelname + '_roi', make_plots=True)

    gta.print_roi()
    gta.print_params()

    lnl1 = -gta.like()
    gta.logger.info('{0:s} Model Likelihood: {1:f}'.format(modelname, lnl1))
    gta.logger.info('{0:s} Model Likelihood Delta: {1:f}'.format(modelname, lnl1-lnl0))
    
    # TS and Residual Maps
    if create_maps:
        # resid map broken :-(
        # hopefully not anymore?
        gta.residmap(modelname, model=model3,
                     loge_bounds=loge_bounds,
                     make_plots=True)
        gta.tsmap(modelname, model=model_pl20,
                                    loge_bounds=loge_bounds, make_plots=True)
        gta.tsmap(modelname, model=model_pl27,
                  loge_bounds=loge_bounds, make_plots=True)

        gta.tsmap('{0:s}_nosource'.format(modelname),
                  model=model_pl20, exclude=[src_name],
                  loge_bounds=loge_bounds, make_plots=True)
        gta.tsmap('{0:s}_nosource'.format(modelname),
                  model=model_pl27, exclude=[src_name],
                  loge_bounds=loge_bounds, make_plots=True)
        #maps_model4_nosource = gta.tsmap('%s_nosource'%modelname,
        #                                 model=model4, exclude=[src_name],
        #                                 loge_bounds=loge_bounds, make_plots=True)    
                                                                        
    # SED Analysis
    if create_sed:
        gta.sed(src_name, outfile=modelname + '_sed_fixed',
                prefix=modelname + '_fixed',
                make_plots=True)
    
        gta.sed(src_name, outfile=modelname + '_sed',
                prefix=modelname,
                free_radius=1.0,
                make_plots=True)    

        gta.sed(src_name,outfile=modelname + '_sed_bin4',
                prefix=modelname + '_bin4',
                loge_bins=gta.log_energies[::2],
                free_radius=1.0,
                make_plots=True)
    
    psf_syst_scale = np.array([0.05,0.05,0.2])
    psf_fnlo = ([3.0,4.0,5.5],list(-1.0*psf_syst_scale))
    psf_fnhi = ([3.0,4.0,5.5],list(1.0*psf_syst_scale))

    # -----------------------------------------------------------------
    # Gaussian Analysis
    # -----------------------------------------------------------------
    kw = dict(spatial_model='RadialGaussian',
              outfile=modelname + '_ext_gauss_ext',
              prefix=modelname + '_gauss',
              psf_scale_fn=None
              )

    kwargs.update(kw)

    gta.logger.info("Starting extension analysis for radial Gaussian")

    gta.logger.info("Testing extension for nominal psf")
    gta.extension(src_name, **kwargs)

    gta.logger.info("Testing extension for psf low")

    #kwargs.update(dict(
    #                   outfile=modelname + '_ext_gauss_ext_psflo',
    #                   prefix=modelname + '_gauss_psflo',
    #                   psf_scale_fn=psf_fnlo
    #                   )
    #)
    #gta.extension(src_name, **kwargs)
    # Matt's way
    gta.extension(src_name,
                  outfile=modelname + '_ext_gauss_ext_psflo',
                  prefix=modelname + '_gauss_psflo',
                  psf_scale_fn=psf_fnlo,
                  spatial_model='RadialGaussian',
                  free_radius=1.0,
                  fit_position=False,
                  free_background=False,
                  update=False,
                  make_tsmap=False)

    gta.logger.info("Testing extension for nominal hi")

    #kwargs.update(dict(
    #                   outfile=modelname + '_ext_gauss_ext_psfhi',
    #                   prefix=modelname + '_gauss_psfhi',
    #                   psf_scale_fn=psf_fnhi
    #                   )
    #)

    # Matt's way
    gta.extension(src_name,
    #              **kwargs)
                  outfile=modelname + '_ext_gauss_ext_psfhi',
                  prefix=modelname + '_gauss_psfhi',
                  psf_scale_fn=psf_fnhi,
                  spatial_model='RadialGaussian',
                  free_radius=1.0,
                  fit_position=False,
                  free_background=False,
                  update=False,
                  make_tsmap=False)

    gta.free_source(src_name)
    gta.fit()
    gta.update_source(src_name,reoptimize=True)
    gta.print_roi()
    gta.print_params()
    
    gta.sed(src_name,outfile=modelname + '_ext_gauss_sed',
            prefix=modelname + '_gauss',
            free_radius=free_radius_sed,
            make_plots=True)
    gta.sed(src_name,outfile=modelname + '_ext_gauss_sed_bin4',
            prefix=modelname + '_gauss_bin4', loge_bins=gta.log_energies[::2],
            free_radius=free_radius_sed,
            make_plots=True)
    gta.write_roi(modelname + '_ext_gauss_roi')

    gta.tsmap(modelname + '_ext_gauss', model=model_pl20,
              loge_bounds=loge_bounds, make_plots=True)
    gta.tsmap(modelname + '_ext_gauss', model=model_pl27,
              loge_bounds=loge_bounds, make_plots=True)

    # -----------------------------------------------------------------
    # Disk Analysis
    # -----------------------------------------------------------------
    gta.logger.info("Starting extension analysis for radial disk")
    gta.load_roi(modelname + '_roi')
    gta.reload_source(src_name)

    kw = dict(spatial_model='RadialDisk',
              outfile=modelname + '_ext_disk_ext',
              prefix=modelname + '_disk',
              psf_scale_fn=None
              )

    kwargs.update(kw)
    
    gta.extension(src_name, **kwargs)

    gta.logger.info("Testing extension for psf low")

#    kwargs.update(dict(
#                       outfile=modelname + '_ext_disk_ext_psflo',
#                       prefix=modelname + '_disk_psflo',
#                       psf_scale_fn=psf_fnlo
#                       )
#    )
#    gta.extension(src_name, **kwargs)
#
#    gta.logger.info("Testing extension for nominal hi")
#
#    kwargs.update(dict(
#                       outfile=modelname + '_ext_disk_ext_psfhi',
#                       prefix=modelname + '_disk_psfhi',
#                       psf_scale_fn=psf_fnhi
#                       )
#    )
#    gta.extension(src_name, **kwargs)

    gta.extension(src_name,
                  outfile=modelname + '_ext_disk_ext_psflo',
                  prefix=modelname + '_disk_psflo',
                  psf_scale_fn=psf_fnlo,
                  spatial_model='RadialDisk',
                  free_radius=1.0,
                  fit_position=False,
                  free_background=False,
                  update=False,
                  make_tsmap=False)

    gta.logger.info("Testing extension for nominal hi")

    # Matt's way
    gta.extension(src_name,
    #              **kwargs)
                  outfile=modelname + '_ext_disk_ext_psfhi',
                  prefix=modelname + '_disk_psfhi',
                  psf_scale_fn=psf_fnhi,
                  spatial_model='RadialDisk',
                  free_radius=1.0,
                  fit_position=False,
                  free_background=False,
                  update=False,
                  make_tsmap=False)
 
    gta.free_source(src_name)
    gta.fit()
    gta.update_source(src_name,reoptimize=True)
    gta.print_roi()
    gta.print_params()
    
    gta.sed(src_name,outfile=modelname + '_ext_disk_sed',
            prefix=modelname + '_disk',
            free_radius=free_radius_sed,
            make_plots=True)
    gta.sed(src_name,outfile=modelname + '_ext_disk_sed_bin4',
            prefix=modelname + '_disk_bin4', loge_bins=gta.log_energies[::2],
            free_radius=free_radius_sed, make_plots=True)
    gta.write_roi(modelname + '_ext_disk_roi')
    
    gta.load_roi(modelname + '_roi')
    gta.reload_source(src_name)    
    gta.logger.info('Finished Region Fit {0:s}'.format(modelname))

def plot_lnlscan(ext, **kwargs):

    if np.all(np.isnan(ext['ebin_loglike'])):
        raise ValueError("All likelihood values are nan, analysis wasn't carried out with ebin extension option")

    from fermipy.my_plotting import truncate_colormap
    from fermipy.utils import init_matplotlib_backend
    from scipy import interpolate
    init_matplotlib_backend()
    import matplotlib.pyplot as plt

    ax = kwargs.pop('ax', plt.gca())
    llhcut = kwargs.pop('llhcut', -2.70)
    cmap = kwargs.pop('cmap', 'BuGn')
    cmap_trunc_lo = kwargs.pop('cmap_trunc_lo', None)
    cmap_trunc_hi = kwargs.pop('cmap_trunc_hi', None)

    ylim = kwargs.pop('ylim', None)

    m = ext['width'] > 0.
    if ylim is None:
        wmin, wmax = ext['width'][m].min(), ext['width'][m].max()
    else:
        wmin, wmax = ylim

    widthM = np.arange(np.log10(wmin), np.log10(wmax), 0.01)
    wbins = len(widthM)
    llhMatrix = np.zeros((len(ext['ebin_e_ctr']), wbins))

    # loop over energy bins
    for i in range(len(ext['ebin_e_ctr'])):
        logl = ext['ebin_loglike'][i][m]
        logl -= np.max(logl)
        try:
            fn = interpolate.interp1d(np.log10(ext['width'][m]), logl, fill_value='extrapolate')
            logli = fn(widthM)
        except:
            logli = np.interp(np.log10(widthM), ext['width'][m], logl)
        llhMatrix[i, :] = logli

    cmap = copy.deepcopy(plt.cm.get_cmap(cmap))
    # cmap.set_under('w')

    if cmap_trunc_lo is not None or cmap_trunc_hi is not None:
        cmap = truncate_colormap(cmap, cmap_trunc_lo, cmap_trunc_hi, 1024)

    xedge = np.insert(ext['ebin_e_max'], 0, ext['ebin_e_min'][0])
    yedge = np.logspace(np.log10(wmin), np.log10(wmax), wbins)
    xedge, yedge = np.meshgrid(xedge, yedge)
    im = ax.pcolormesh(xedge, yedge, llhMatrix.T,
                       vmin=llhcut,
                       vmax=0,
                       cmap=cmap,
                       linewidth=0)

    ax.set_ylim(wmin, wmax)
    #plt.gca().set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(ext['ebin_e_min'][0], ext['ebin_e_max'][-1])

    if kwargs.pop('plotcb', True):
        cb = plt.colorbar(im)
        cb.set_label('Delta LogLikelihood')

    return ax, im

def rebin_ebin_logl(ext_old, rebin=2):
# TODO rebin > 2 not working, loglike and e bins need to be calculated differently

    ext = copy.deepcopy(ext_old)
    ext['ebin_e_min'] = ext['ebin_e_min'][::rebin]  
    ext['ebin_e_max'] = ext['ebin_e_max'][rebin-1::rebin] 

    if rebin < 3:
        ext['ebin_loglike'] = ext['ebin_loglike'][1::rebin,:] + ext['ebin_loglike'][::rebin,:]
    else:
        ebin_loglike = np.zeros((ext['ebin_loglike'].shape[0] / rebin, ext['ebin_loglike'].shape[1]))
        for i in range(ebin_loglike.shape[0]):
            ebin_loglike[i] = np.sum(ext['ebin_loglike'][i * rebin : (i + 1) * rebin], axis=0)
        ext['ebin_loglike'] = ebin_loglike
    ext['ebin_e_ctr'] = np.sqrt(ext['ebin_e_min'] * ext['ebin_e_max'])    
    ext['ebin_ts_ext'] = -2. * (ext['ebin_loglike'][:,0] - ext['ebin_loglike'].max(axis=1))       
    ext['ebin_ext'] = ext['width'][np.argmax(ext['ebin_loglike'], axis=1)]
    dloglike = (ext['ebin_loglike'].T - ext['ebin_loglike'].max(axis=1)).T

    # recalculate errors for each energy bin
    m = ext['width'] > 0.
    wmin, wmax = ext['width'][m].min(), ext['width'][m].max()
    widthM = np.arange(np.log10(wmin), np.log10(wmax), 0.01)
    ext['ebin_ext_err'] = np.zeros(ext['ebin_ext'].size)
    ext['ebin_ext_err_hi'] = np.zeros(ext['ebin_ext'].size)
    ext['ebin_ext_err_lo'] = np.zeros(ext['ebin_ext'].size)
    ext['ebin_ext_ul95'] = np.zeros(ext['ebin_ext'].size)
    for i in range(len(ext['ebin_e_ctr'])):

        for j, sl in enumerate([slice(1, np.argmax(dloglike[i])), slice(np.argmax(dloglike[i]),dloglike[i].size)]):

            if not j:
                widthM = np.arange(np.log10(wmin), np.log10(ext['width'][m][np.argmax(dloglike[i])]), 0.01)
            else:
                widthM = np.arange(np.log10(ext['width'][np.argmax(dloglike[i])]), np.log10(wmax), 0.01)

            try:
                fn = USpline(np.log10(ext['width'][sl]), dloglike[i][sl], ext='extrapolate', s=0, k=2)
                logli = fn(widthM)
            except:
                logli = np.interp(np.log10(widthM), ext['width'][sl], dloglike[i][sl])

            if not j:
                ext['ebin_ext_err_lo'][i] = ext['ebin_ext'][i] - 10.**widthM[np.argmin(np.abs(logli + 1. / 2.))]
            else:
                ext['ebin_ext_err_hi'][i] = 10.**widthM[np.argmin(np.abs(logli + 1. / 2.))] - ext['ebin_ext'][i]
                ext['ebin_ext_ul95'][i] = 10.**widthM[np.argmin(np.abs(logli + 2.71 / 2.))]
            #if j and i:
                #1. / 0.
    ext['ebin_ext_err'] = 0.5 * (ext['ebin_ext_err_hi'] + ext['ebin_ext_err_lo'])
    return ext



def plot_ext_ebin_logl(ext, ts_thr=4., fig=None, ax=None, plot_lnl=True, lnl_kwargs={}):
    from fermipy.utils import init_matplotlib_backend
    init_matplotlib_backend()
    import matplotlib.pyplot as plt

    m = ext['ebin_ts_ext'] > ts_thr

    if fig is None:
        fig = plt.figure()

    ectr = ext['ebin_e_ctr']
    delo = ext['ebin_e_ctr'] - ext['ebin_e_min']
    dehi = ext['ebin_e_max'] - ext['ebin_e_ctr']
    xerr0 = np.vstack((delo[m], dehi[m]))
    xerr1 = np.vstack((delo[~m], dehi[~m]))

    if ax is None:
        ax = plt.gca()

    ax.errorbar(ectr[m], ext['ebin_ext'][m], xerr=xerr0,
                yerr=(ext['ebin_ext_err_lo'][m],
                      ext['ebin_ext_err_hi'][m]),
                color='k',
                linestyle='None',
                marker='o')

    ax.errorbar(ectr[~m], ext['ebin_ext_ul95'][~m], xerr=xerr1,
                yerr=0.2 * ext['ebin_ext_ul95'][~m],
                uplims=True,
                color='k',
                linestyle='None', marker='o')
    ax.set_xlabel('Energy [log$_{10}$(E/MeV)]')
    ax.set_ylabel('Extension [deg]')

    ymin = min(10**-1.5, 0.8 * ext['ext_ul95'])
    ymax = max(10**-0.5, 1.2 * ext['ext_ul95'])
    if np.any(np.isfinite(ext['ebin_ext_ul95'])):
        ymin = min(ymin, 0.8 * np.nanmin(ext['ebin_ext_ul95']))
        ymax = max(ymax, 1.2 * np.nanmax(ext['ebin_ext_ul95']))

    if plot_lnl:
        lnl_kwargs['ylim'] = (ymin, ymax)
        _, im = plot_lnlscan(ext, **lnl_kwargs)
    else:
        im = None

    if ext['ts_ext'] > ts_thr:
        ax.axhline(ext['ext'], color='k')
        ext_lo = ext['ext'] - ext['ext_err_lo']
        ext_hi = ext['ext'] + ext['ext_err_hi']
        ax.fill_between([ext['ebin_e_min'][0], ext['ebin_e_max'][-1]],
                        [ext_lo, ext_lo], [ext_hi, ext_hi],
                        alpha=0.5, color='k')

        ymin = min(ymin, 0.8 * (ext['ext'] - ext['ext_err_lo']))
        ymax = max(ymax, 1.2 * (ext['ext'] + ext['ext_err_hi']))

    else:
        ax.axhline(ext['ext_ul95'], color='k', linestyle='--')

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(ext['ebin_e_min'][0], ext['ebin_e_max'][-1])
    ax.set_xscale('log')
    ax.set_yscale('log')

    return fig, ax, im


def fit_halo_sed(gta, modelname, src_name,
                 halo_width=np.logspace(-1.5, 0.25, 15),
                 spatial_model='RadialGaussian',
                 halo_source_dict = halo_source_dict_default,
                 cov_scale=5.,
                 distance_free_norm=1.,
                 loge_bounds=None):
    """
    Fit a halo with free spectral shape on top of a point source
    and extract SED. Do this for an array of assumed halo widths.

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        the gt analysis object

    modelname: str
        base name of the model

    src_name: str
        name of the source for which extension will be tested

    halo_width: array-like
        assumed with of the halo. The elements of the array 
        will be passed to the SpatialWidth parameter of the 
        chosen spatial model.

    {options}

    halo_source_dict:
        halo source dictionary for the spectral parameters

    cov_scale: float
        the cov_scale used for the SED analysis. 
        Default: 5.

    distance_free_norm: float
        distance in deg up to which normalizations are left 
        free in fit. Default: 1.

    spatial_model: str
        name of spatial model assumed for the halo. Default: RadialGausssian

    loge_bounds array-like or None
        if not None, use these energy bounds for the analysis. 
        Otherwise, use the energy bounds of the gta object
        Default: None
    """

    gta.logger.info('Starting Halo SED Fit {0:s}'.format(modelname))
    
    halo_source_name = 'halo_' + spatial_model
    halo_source_dict['SpatialModel'] = spatial_model
    halo_source_dict['SpatialWidth'] = 1.

    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

    #gta.load_roi(modelname)
    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
        
    gta.free_sources(False)

    # define the base model, 
    # only let normalizations of innermost sources free
    # and exclude background models
    gta.free_sources(distance=distance_free_norm,
                     pars='norm',
                     exclude=diff_sources)
    gta.write_xml(modelname + '_base')
    
    # loop over halo width
    for i, w in enumerate(halo_width):

        halo_source_dict['SpatialWidth'] = w        
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name, halo_source_dict, free=True)

        # Do one fit with index free
        gta.set_parameter(halo_source_name,'Index',-2.0,
                          update_source=False)

        gta.fit()
        
        # SED w/ Index = 2.0
        gta.sed(halo_source_name,
                prefix='{0:s}_{1:02n}'.format(modelname,i),
                fix_background=False, cov_scale=cov_scale)

        gta.write_roi('{0:s}_halo_gauss_sed_{1:02n}'.format(modelname,i),
                      make_plots=False)

    gta.logger.info('Finished Halo SED Fit %s'%(modelname))

def fit_halo_scan(gta, modelname,
                  src_name,
                  halo_width=np.logspace(-1.5, 0.25, 15),
                  halo_index=np.arange(1., 4.25, 0.25),
                  halo_source_dict=halo_source_dict_default,
                  spatial_model='RadialGaussian',
                  index_par_name='Index',
                  cov_scale=5.,
                  free_radius_sed=1.,
                  distance_free_norm=1.,
                  loge_bounds=None, optimizer='NEWTON'):
    """
    Fit a halo with free spectral shape on top of a point source.
    Extract the SED over fixed values of the halo width and spectral 
    index.

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        the gt analysis object

    modelname: str
        name of the model that will be reloaded

    src_name: str
        name of the source for which extension will be tested

    halo_width: array-like
        assumed width of the halo. The elements of the array 
        will be passed to the SpatialWidth parameter of the 
        chosen spatial model.

    halo_index: array-like
        assumed spectral index of the halo. The elements of the array 
        will be passed to the parameter of the 
        chosen spectral model given by the index_par_name parameter

    {options}

    halo_source_dict:
        halo source dictionary for the spectral parameters

    cov_scale: float
        the cov_scale used for the SED analysis. 
        Default: 5.

    index_par_name: str
        name of the parameter that will be changed by halo_index values
        Default: Index

    spatial_model: str
        name of spatial model assumed for the halo. Default: RadialGausssian

    loge_bounds array-like or None
        if not None, use these energy bounds for the analysis. 
        Otherwise, use the energy bounds of the gta object
        Default: None

    distance_free_norm: float
        distance in deg up to which normalizations are left 
        free in fit. Default: 1.

    free_radius_sed: float
        distance in deg up to which normalizations are left 
        free in SED extraction. Default: 1.

    optimizer: str
        name of optimizer to be used. Default: NEWTON
    """

    gta.logger.info('Starting Halo Scan %s'%(modelname))

    halo_source_name = 'halo_' + spatial_model
    halo_source_dict['SpatialModel'] = spatial_model
    halo_source_dict['SpatialWidth'] = 1.

    outprefix = '_'.join([modelname, halo_source_name])
    
    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

    #gta.load_roi(modelname)
    #if loge_bounds is not None:
    #    gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    skydir = gta.roi[src_name].skydir
    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(skydir=skydir,
                     distance=distance_free_norm,
                     pars='norm',
                     exclude=diff_sources)
    gta.write_xml(modelname + '_base')

    halo_tab = gta.roi.create_table([])
    halo_tab_idx_free = gta.roi.create_table([])
    halo_data = []
    halo_data_idx_free = []
        
    # loop over halo width
    for i, w in enumerate(halo_width):

        gta.logger.info('Fitting Halo Width {0:.3f}'.format(w))
        
        halo_source_dict['SpatialWidth'] = w
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name, halo_source_dict, free=True)

        # Free Index
        gta.free_norm(halo_source_name)
        gta.fit(optimizer=optimizer)

        # SED should be indendepnt of Index
        gta.sed(halo_source_name,
                prefix='{0:s}_cov_{1:02n}'.format(modelname, i),
                outfile='{0:s}_cov_{1:02n}'.format(outprefix, i),
                free_radius=free_radius_sed,
                cov_scale=cov_scale,
                optimizer={'optimizer' : 'MINUIT'},
                make_plots=False)

        gta.sed(halo_source_name,
                prefix='{0:s}_covx2_{1:02n}'.format(modelname, i),
                outfile='{0:s}_covx2_{1:02n}'.format(outprefix, i),
                free_radius=free_radius_sed,
                cov_scale=2. * cov_scale,
                optimizer={'optimizer' : 'MINUIT'},
                make_plots=False)
        
        gta.free_parameter(halo_source_name, index_par_name)
        gta.fit(optimizer=optimizer)
        gta.free_parameter(halo_source_name, index_par_name, False)
        gta.update_source(halo_source_name,reoptimize=True,
                          optimizer={'optimizer' : optimizer})

        halo_data_idx_free += [copy.deepcopy(gta.roi[halo_source_name].data)]
        gta.roi[halo_source_name].add_to_table(halo_tab_idx_free)
        gta.write_roi('{0:s}_{1:02n}'.format(outprefix,i),
                      make_plots=False)

        gta.print_params(loglevel=logging.DEBUG)
        
        # Scan over fixed index
        for j, idx in enumerate(halo_index):

            gta.logger.info('Fitting Halo Index {0:.3f}'.format(idx))
            
            model_idx = i*len(halo_index) + j
            gta.set_norm(halo_source_name, 0.1, update_source=False)            
            gta.set_parameter(halo_source_name, index_par_name, -1.0*idx,
                              update_source=False)
            
            gta.fit(update=False, optimizer=optimizer)

            gta.print_params(loglevel=logging.DEBUG)
            
            gta.update_source(halo_source_name,reoptimize=True,
                              optimizer={'optimizer' : optimizer})

            ul_flux = get_parameter_limits(gta.roi[halo_source_name]['flux_scan'],
                                           gta.roi[halo_source_name]['loglike_scan'])
            ul_eflux = get_parameter_limits(gta.roi[halo_source_name]['eflux_scan'],
                                            gta.roi[halo_source_name]['loglike_scan'])

            gta.roi[halo_source_name]['flux_err'] = ul_flux['err']
            gta.roi[halo_source_name]['eflux_err'] = ul_eflux['err']
            
            gta.logger.info('%{0:s} Halo Width: {1:6.3f} Index: {2:6.2f} TS: {3:6.2f} Flux: {4:8.4g}'.format(
                                modelname,w,idx,
                                gta.roi[halo_source_name]['ts'],
                                gta.roi[halo_source_name]['flux'])
            )
    
            #gta.write_roi('%s_%02i_%02i'%(outprefix,i,j),make_plots=False)
            halo_data += [copy.deepcopy(gta.roi[halo_source_name].data)]
            gta.roi[halo_source_name].add_to_table(halo_tab)
            
        gta.delete_source(halo_source_name,save_template=False) 

    np.save(os.path.join(gta.workdir, '{0:s}_data.npy'.format(outprefix)), halo_data)
    np.save(os.path.join(gta.workdir, '{0:s}_data_idx_free.npy'.format(outprefix)),
            halo_data_idx_free)

    tab_halo_width, tab_halo_index = np.meshgrid(halo_width,halo_index,indexing='ij')
    halo_tab['halo_width'] = np.ravel(tab_halo_width)
    halo_tab['halo_index'] = np.ravel(tab_halo_index)
    halo_tab_idx_free['halo_width'] = halo_width

    stack_files(sorted(glob.glob(os.path.join(gta.workdir,'{0:s}*cov_*fits'.format(outprefix)))),
                os.path.join(gta.workdir,'{0:s}_cov_sed.fits'.format(outprefix)),
                new_cols=[Column(name='halo_width',data=halo_width, unit='deg')])

    stack_files(sorted(glob.glob(os.path.join(gta.workdir,'{0:s}*covx2_*fits'.format(outprefix)))),
                os.path.join(gta.workdir,'{0:s}_covx2_sed.fits'.format(outprefix)),
                new_cols=[Column(name='halo_width',data=halo_width, unit='deg')])
    
    halo_tab.write(os.path.join(gta.workdir,'{0:s}_data.fits'.format(outprefix)),
                   overwrite=True)
    halo_tab_idx_free.write(os.path.join(gta.workdir,'{0:s}_data_idx_free.fits'.format(outprefix)),
                            overwrite=True)
    gta.logger.info('Finished Halo Scan {0:s}'.format(modelname))
    
    
def fit_halo(gta, modelname, src_name,
             spatial_model='RadialGaussian',
             halo_source_dict=halo_source_dict_default,
             loge_bounds=None,
             distance_free_norm=1.,
             free_ext_radius=1., 
             n_step=1,
             dlog_thr=0.5,
             index_par_name="Index",
             optimizer='NEWTON'):
    """
    Fit a halo with free spectral shape on top of a point source.

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        the gt analysis object

    modelname: str
        name of the model that will be reloaded

    src_name: str
        name of the source for which extension will be tested

    {options}

    halo_source_dict:
        halo source dictionary for the spectral parameters

    spatial_model: str
        name of spatial model assumed for the halo. Default: RadialGausssian

    loge_bounds array-like or None
        if not None, use these energy bounds for the analysis. 
        Otherwise, use the energy bounds of the gta object
        Default: None

    index_par_name: str
        name of the parameter that will be changed by halo_index values
        Default: Index

    optimizer: str
        name of optimizer to be used. Default: NEWTON

    distance_free_norm: float
        distance in deg up to which normalizations are left 
        free in fit. Default: 1.

    free_ext_radius: float
        distance in deg up to which normalizations are left 
        free in extension fit. Default: 1.

    n_step: int
        number of times spectral and extension fits are repeated
        default: 1

    dlog_thr: float
        dloglike threshold below which fit loop is stopped.
        Default: 0.5
    """

    gta.logger.info('Starting Halo Fit %s'%(modelname))
    halo_source_dict['SpatialModel'] = spatial_model
    halo_source_dict['SpatialWidth'] = 1.
    
    halo_source_name = 'halo_' + spatial_model

    outprefix = '_'.join([modelname, halo_source_name])
    
    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

#    gta.load_roi(modelname)
#    if loge_bounds is not None:
#        gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(distance=distance_free_norm,
                     pars='norm',
                     exclude=diff_sources)

    # Find best-fit halo model
    halo_source_dict['SpatialWidth'] = 0.1
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_norm(halo_source_name)

    i_step = 0
    dlog = 0.
    lnl0 = 0.

    while dlog > dlog_thr and i_step <= nstep:
        gta.extension(halo_source_name,update=True,
                      optimizer={'optimizer' : optimizer},
                      free_radius=free_ext_radius)

        # Fit spectrum
        gta.free_parameter(halo_source_name, index_par_name)
        gta.fit()

        dlog = np.abs(gta.like() - lnl0)
        i_step += 1

        gta.logger("Finished with fit step {0:n}/{1:n}, |dloglike| = {2:.2f}".format(i_step, n_step, dlog))

    gta.update_source(halo_source_name,reoptimize=True,
                      optimizer={'optimizer' : optimizer})

    gta.print_params()
    
    gta.write_roi(outprefix, make_plots=False)    
    np.save(os.path.join(gta.workdir,'{0:s}_data.npy'.format(outprefix)),
            copy.deepcopy(gta.roi[halo_source_name].data))
    gta.delete_source(halo_source_name,save_template=False)
    
    gta.logger.info('Finished Halo Fit {0:s}'.format(modelname))

#def plot_halo_loglike_profile(loglike, halo_width, halo_index): 
def interp_2dloglike_profile(x,y,z, oversample=10, **spline_kw):
    spline_kw.setdefault('kx', 1)
    spline_kw.setdefault('ky', 1)
    spline_kw.setdefault('s', 0.)

    spline = RBSpline(x[:,0], y[0,:], z, **spline_kw)

    x_interp = np.linspace(x[0,0], x[-1,0], x.shape[0] * oversample)
    y_interp = np.linspace(y[0,0], y[0,-1], y.shape[1] * oversample)
    xx, yy = np.meshgrid(x_interp, y_interp, indexing='ij')
    zz = spline(x_interp, y_interp)
    return xx, yy, zz

def get_halo_results(halo_data_fits_file, oversample=1, ax=None, vmin=None, vmax=0., cmap=None):
    import matplotlib.pyplot as plt
    from fermipy.utils import init_matplotlib_backend
    init_matplotlib_backend()

    if cmap is None:
        cmap= plt.cm.PuBuGn
        cmap.set_under('1.')
        cmap.set_over(cmap(1.))

    t = Table.read(halo_data_fits_file)
    t_idx_free = Table.read(halo_data_fits_file.replace('data', 'data_idx_free'))

    # generate the arrays for halo with and index
    halo_index_shape = np.unique(t['halo_index']).size
    halo_width_shape = np.unique(t['halo_width']).size

    halo_index = t['halo_index'].data.reshape((halo_width_shape, halo_index_shape))
    halo_width = t['halo_width'].data.reshape((halo_width_shape, halo_index_shape))

    idx = np.argmax(t['loglike'])
    idx_idx_free = np.argmax(t_idx_free['loglike'])
    loglike = t['loglike'].data.reshape((halo_width_shape, halo_index_shape))
    logL0 = loglike - 0.5 * t['ts'].data.reshape((halo_width_shape, halo_index_shape))

    if oversample > 1:
        _, _, logL0 = interp_2dloglike_profile(halo_width, halo_index, logL0,
                                               oversample=oversample,
                                               kx=2, ky=2)
        halo_width, halo_index, loglike = interp_2dloglike_profile(halo_width, halo_index, loglike,
                                                                   oversample=oversample,
                                                                   kx=2, ky=2)
    # store best fit results in a dict
    res = {'loglike': loglike,
           'logL0': logL0,
           'halo_index': halo_index,
           'halo_width': halo_width,
           'max_ts': t['ts'][idx], 
           'max_ts_halo_width': t['halo_width'][idx], 
           'max_ts_halo_index': t['halo_index'][idx], 
           'flux': t['flux1000'][idx], 
           'flux_err': t['flux1000_err'][idx], 
           'flux_ul95': t['flux1000_ul95'][idx], 
           'max_ts_idx_free': t_idx_free['ts'][idx_idx_free], 
           'max_ts_halo_width': t_idx_free['halo_width'][idx_idx_free],
           'max_ts_halo_index': t_idx_free['param_values'][idx_idx_free][0], 
           'max_ts_halo_index_err': t_idx_free['param_errors'][idx_idx_free][1]
    }
    if oversample > 1:
        res['max_ts_halo_width_interp'] = res['halo_width'][:,0][np.where(res['loglike'] == res['loglike'].max())[0][0]]
        res['max_ts_halo_index_interp'] = res['halo_index'][0,:][np.where(res['loglike'] == res['loglike'].max())[1][0]]
    else:
        res['max_ts_halo_width_interp'] = res['max_ts_halo_width']
        res['max_ts_halo_index_interp'] = res['max_ts_halo_index']

    # plot 2D loglike profile
    if ax is None:
        ax = plt.gca()

    # contour levels for two-sided 68%, 90%, 95% C.L. for two degrees of freedom
    levels = np.array([2.30, 4.61, 5.99]) / 2. * -1.

    if t['ts'][idx] > 4.61:
        # contour levels for two-sided 68%, 90%, 95% C.L. for two degrees of freedom
        if vmin is None:
            vmin = -5.99 / 2.

    im = ax.pcolormesh(halo_width, halo_index, loglike - loglike.max(),
                        vmax=vmax, vmin=vmin, 
                        cmap=cmap
                        )
    plt.colorbar(im, label='$\Delta\log\mathcal{L}$')


    c = ax.contour(halo_width, halo_index, loglike - loglike.max(),
                   levels[::-1],
                   colors='k',
                   linestyles='-',
                   linewidths=1.,
                   )

    ax.plot(res['max_ts_halo_width_interp'],
            res['max_ts_halo_index_interp'],
            marker='x', color='r')
    ax.clabel(c, fmt="%.2f", ha='center', va='center')
    ax.set_ylabel("Index")
    ax.set_xlabel("$R_{68}$ (deg)")
    return res, ax

def profile_norm_tied_parameter(gta, name, parName,
                                tied_src_name, tied_par_name,
                                logemin=None, logemax=None,
                                reoptimize=False,
                                xvals=None, npts=None, savestate=True, **kwargs):
    """
    Helper function for fit_igmf_halo_scan.
    The central source normalization is profiled over 
    and a tied parameter (the halo template prefactor)
    is rescaled accordingly.

    Parameters
    ----------

    name : str
       Source name.

    parName : str
       Parameter name.

    tied_src_name: str
        Source name with tied parameter

    tied_par_name: str
        Name of tied parameter

    reoptimize : bool
       Re-fit nuisance parameters at each step in the scan.  Note
       that enabling this option will only re-fit parameters that
       were free when the method was executed.

    Returns
    -------

    lnlprofile : dict
       Dictionary containing results of likelihood scan.

    """
    # Find the source
    name = gta.roi.get_source_by_name(name).name

    par = gta.like.normPar(name)
    parName = gta.like.normPar(name).getName()
    idx = gta.like.par_index(name, parName)

    # find the tied source
    tied_src_name = gta.roi.get_source_by_name(tied_src_name).name
    tied_par = gta.like.normPar(tied_src_name)
    tied_par_name = gta.like.normPar(tied_src_name).getName()
    tied_idx = gta.like.par_index(tied_src_name, tied_par_name)
    tied_bounds = gta.like.model[tied_idx].getBounds()
    tied_value = gta.like.model[tied_idx].getValue()

    bounds = gta.like.model[idx].getBounds()
    value = gta.like.model[idx].getValue()
    loge_bounds = gta.loge_bounds

    optimizer = kwargs.get('optimizer', gta.config['optimizer'])

    if savestate:
        saved_state = gta._latch_state()

    # If parameter is fixed temporarily free it
    par.setFree(True)
    if optimizer['optimizer'] == 'NEWTON':
        gta._create_fitcache()

    if logemin is not None or logemax is not None:
        loge_bounds = gta.set_energy_range(logemin, logemax)
    else:
        loge_bounds = gta.loge_bounds

    loglike0 = -gta.like()

    if xvals is None:

        err = par.error()
        val = par.getValue()
        if err <= 0 or val <= 3 * err:
            xvals = 10 ** np.linspace(-2.0, 2.0, 51)
            if val < xvals[0]:
                xvals = np.insert(xvals, val, 0)
        else:
            xvals = np.linspace(0, 1, 25)
            xvals = np.concatenate((-1.0 * xvals[1:][::-1], xvals))
            xvals = val * 10 ** xvals

    if np.isnan(xvals).any():
        raise RuntimeError(
            "Parameter scan points for %s::%s include infinite value." % (name, parName))

    tied_xvals = tied_value * xvals / value

    # Update parameter bounds to encompass scan range
    try:
        gta.like[idx].setBounds(min(min(xvals), value, bounds[0]),
                                max(max(xvals), value, bounds[1]))

        gta.like[tied_idx].setBounds(min(min(tied_xvals), tied_value, tied_bounds[0]),
                                     max(max(tied_xvals), tied_value, tied_bounds[1]))
    except RuntimeError:
        gta.logger.warning(
            "Caught failure on setBounds for %s::%s." % (name, parName))

    o = {'xvals': xvals,
         'npred': np.zeros(len(xvals)),
         'npred_wt': np.zeros(len(xvals)),
         'dnde': np.zeros(len(xvals)),
         'flux': np.zeros(len(xvals)),
         'eflux': np.zeros(len(xvals)),
         'dloglike': np.zeros(len(xvals)),
         'loglike': np.zeros(len(xvals)),
         'xvals_tied': tied_xvals,
         'npred_tied': np.zeros(len(xvals)),
         'npred_wt_tied': np.zeros(len(xvals)),
         'dnde_tied': np.zeros(len(xvals)),
         'flux_tied': np.zeros(len(xvals)),
         'eflux_tied': np.zeros(len(xvals))
         }

    if reoptimize and hasattr(gta.like.components[0].logLike,
                              'setUpdateFixedWeights'):

        for c in gta.components:
            c.like.logLike.setUpdateFixedWeights(False)

    for i, x in enumerate(xvals):

        try:
            gta.like[idx] = x
            gta.like[tied_idx] = tied_xvals[i]
        except RuntimeError:
            gta.logger.warning(
                "Caught failure on set for %s::%s: %.2f" % (name, parName, x))

        if gta.like.nFreeParams() > 1 and reoptimize:
            # Only reoptimize if not all frozen
            gta.like.freeze(idx)
            gta.like.freeze(tied_idx)
            fit_output = gta._fit(errors=False, **optimizer)
            loglike1 = fit_output['loglike']
            gta.like.thaw(idx)
        else:
            loglike1 = -gta.like()

        o['flux'][i] = gta.like[name].flux(10 ** loge_bounds[0],
                                             10 ** loge_bounds[1])
        o['eflux'][i] = gta.like[name].energyFlux(10 ** loge_bounds[0],
                                                  10 ** loge_bounds[1])

        o['flux_tied'][i] = gta.like[tied_src_name].flux(10 ** loge_bounds[0],
                                                         10 ** loge_bounds[1])
        o['eflux_tied'][i] = gta.like[tied_src_name].energyFlux(10 ** loge_bounds[0],
                                                                10 ** loge_bounds[1])

        o['dloglike'][i] = loglike1 - loglike0
        o['loglike'][i] = loglike1
        o['dnde'][i] = gta.like[idx].getTrueValue()
        o['dnde_tied'][i] = gta.like[tied_idx].getTrueValue()

        cs = gta.model_counts_spectrum(name,
                                       loge_bounds[0],
                                       loge_bounds[1], summed=True)
        o['npred'][i] = np.sum(cs)
        cs_wt = gta.model_counts_spectrum(name,
                                          loge_bounds[0],
                                          loge_bounds[1], summed=True, weighted=True)
        o['npred_wt'][i] = np.sum(cs_wt)

        tied_cs = gta.model_counts_spectrum(tied_src_name,
                                            loge_bounds[0],
                                            loge_bounds[1], summed=True)

        o['npred_tied'][i] = np.sum(tied_cs)
        tied_cs_wt = gta.model_counts_spectrum(tied_src_name,
                                          loge_bounds[0],
                                          loge_bounds[1], summed=True, weighted=True)
        o['npred_wt_tied'][i] = np.sum(tied_cs_wt)

    gta.like[idx] = value
    gta.like[tied_idx] = tied_value

    if reoptimize and hasattr(gta.like.components[0].logLike,
                              'setUpdateFixedWeights'):

        for c in gta.components:
            c.like.logLike.setUpdateFixedWeights(True)

    # Restore model parameters to original values
    if savestate:
        saved_state.restore()

    gta.like[idx].setBounds(*bounds)
    gta.like[tied_idx].setBounds(*tied_bounds)
    if logemin is not None or logemax is not None:
        gta.set_energy_range(*loge_bounds)

    return o

def profile_halo_src_norms(gta, src_name, halo_source_name,
                           scale0, prefactor0,
                           index=1., 
                           src_norm_name='Prefactor',
                           halo_norm_name='Prefactor',
                           src_scale_name='Scale',
                           reoptimize=False,
                           savestate=True,
                           sigma=3.,
                           xsteps=30,
                           ysteps=31):
    # save state
    if savestate:
        #saved_state = LikelihoodState(gta.like)
        saved_state = gta._latch_state()

    s = gta.roi.get_source_by_name(src_name)
    halo = gta.roi.get_source_by_name(halo_source_name)

    loglike = -gta.like()

    # get the right index for the norm parameters of the central source 
    # and igmf halo
    name = gta.roi.get_source_by_name(src_name).name
    halo_name = gta.roi.get_source_by_name(halo_source_name).name
    idx_norm = gta.like.par_index(name, src_norm_name)
    idx_norm_halo = gta.like.par_index(halo_name, halo_norm_name)

    bounds = gta.like.model[idx_norm].getBounds()
    bounds_halo = gta.like.model[idx_norm_halo].getBounds()
    value = gta.like.model[idx_norm].getValue()
    loge_bounds = gta.loge_bounds

    #norm = np.logspace(np.log10(s.spectral_pars['Prefactor']['error'] * sigma/ \
                            #s.spectral_pars['Prefactor']['value']),
                        #np.log10(s.spectral_pars['Prefactor']['value'] / sigma \
                            #/ s.spectral_pars['Prefactor']['error']),
                        #ysteps)

    # norm steps for central source
    norm = np.logspace(np.log10(s.spectral_pars['Prefactor']['value'] / sigma / 2.), \
                       np.log10(s.spectral_pars['Prefactor']['value'] * sigma * 2.), \
                       xsteps)

    gta.logger.info("Performing 2D profile")
    gta.logger.info("norm steps: {0}".format(norm))
    
    m = (norm >= s.spectral_pars['Prefactor']['min']) & \
        (norm <= s.spectral_pars['Prefactor']['max']) 
    if not m.sum():
        raise ValueError("No Prefactor scan points!")
    norm= norm[m]

    output = {'norm_halo': []}

    try:
        gta.like[idx_norm].setBounds(min(min(norm), value, bounds[0]),
                                     max(max(norm), value, bounds[1]))

    except RuntimeError:
        gta.logger.warning(
            "Caught failure on setBounds for %s::%s." % (name, src_norm_name))

    for k in ['flux_src', 'eflux_src', 'dnde_src', 'npred_src']:
        output[k] = np.zeros_like(norm)

    # do the profiling
    for i,x in enumerate(norm):
        if i > 0:
            saved_state.restore()
#            s = gta.roi.get_source_by_name(target)
        #logging.info("{0:n}: {1}".format(i, s.spectral_pars))
        #s.spectral_pars['Index']['value'] = x
        #s.spectral_pars['Index']['free'] = False
        #logging.info("{0:n} set: {1}".format(i, s.spectral_pars))
        #s.set_spectral_pars(s.spectral_pars)
        #gta.roi[target].set_spectral_pars(s.spectral_pars)
        
        try:
            gta.like[idx_norm] = x
            gta.like.freeze(idx_norm)
        except RuntimeError:
            logging.warning(
                "Caught failure on set for %s::%s: %.2f" % (name, parName, x))

        # compute the new norm values for the halo
        max_halo_norm = set_halo_normalization(gta, src_name, halo_source_name,
                                               prefactor0, scale0,
                                               update_halo=True,
                                               injection_norm_name=src_norm_name,
                                               injection_scale_name=src_scale_name,
                                               index=index)
        xvals = np.logspace(-2., 0., ysteps) * max_halo_norm

        o = gta.profile_norm(halo_source_name, xvals=xvals, reoptimize=reoptimize)

        for k in ['loglike', 'flux', 'eflux', 'dnde', 'npred']:
            key = k if k == 'loglike' else k + '_halo'
            if not key in output:
                output[key] = []
            output[key].append(o[k])
        output['norm_halo'].append(xvals)

        output['flux_src'][i] = gta.like[name].flux(10 ** loge_bounds[0],
                                             10 ** loge_bounds[1])
        output['eflux_src'][i] = gta.like[name].energyFlux(10 ** loge_bounds[0],
                                                  10 ** loge_bounds[1])
        output['dnde_src'][i] = gta.like[idx_norm].getTrueValue()

        cs = gta.model_counts_spectrum(name,
                                       loge_bounds[0],
                                       loge_bounds[1], summed=True)
        output['npred_src'][i] = np.sum(cs)

        gta.like.thaw(idx_norm)

    for k in ['loglike', 'flux', 'eflux', 'dnde', 'norm', 'npred']:
        key = k if k == 'loglike' else k + '_halo'
        output[key] = np.array(output[key])

    output['loglike0'] = loglike
    output['dloglike'] = output['loglike'] - loglike
    output['norm_src'] = norm
    output['src_name'] = src_name
    output['halo_name'] = halo_name
    output['src_par_name'] = src_norm_name
    output['halo_par_name'] = halo_norm_name

    gta.logger.info("Finished")
    #np.save(os.path.join(gta.workdir, prefix + '_{0:s}_profile2d'.format(target.replace(' ','').lower())),
            #output)
    # Restore model parameters to original values
    if savestate:
        saved_state.restore()

    gta.like[idx_norm].setBounds(*bounds)
    gta.like[idx_norm_halo].setBounds(*bounds_halo)
    return output

def plot_2dprofile_scan(output, gta=None, fig=None, ax=None, ykey='dloglike', xkey='norm', **kwargs):
    """
    plot the output of the 2d profile scan
    """
    import matplotlib.pyplot as plt
    from scipy import interpolate

    cmap = copy.deepcopy(plt.cm.get_cmap(kwargs.pop('cmap', 'PuBuGn')))
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)

    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111)

    norm_halo = np.logspace(np.log10(output[xkey + '_halo'].min()), 
                            np.log10(output[xkey + '_halo'].max()),
                            100)
    xx, yy = np.meshgrid(output[xkey + '_src'], norm_halo, indexing='ij')

    logl_matrix = np.zeros_like(xx)

    for i in range(output[ykey].shape[0]):
        logl = output[ykey][i]
        norm = output[xkey + '_halo'][i]

        #fn = interpolate.interp1d(np.log10(norm), logl, fill_value='extrapolate')
        fn = USpline(np.log10(norm), logl, ext='extrapolate', s=0, k=2)
        logl_matrix[i] = fn(np.log10(norm_halo))
        m = norm_halo > output[xkey + '_halo'][i].max()
        logl_matrix[i][m] = np.nan

    # compute the likelihood difference yourself
    if ykey == 'loglike':
        m = np.isfinite(logl_matrix)
        logl_matrix[m] -= logl_matrix[m].max()
        print (logl_matrix[m].min(), logl_matrix[m].max(), vmin if vmin is not None else output[ykey].min())
        print (logl_matrix[m].min(), logl_matrix[m].max(), vmax if vmax is not None else output[ykey].max())

    cmap.set_under("1.")

    im = ax.pcolormesh(xx, yy, logl_matrix,
                       vmin=vmin if vmin is not None else output[ykey].min(),
                       vmax=vmax if vmax is not None else output[ykey].max(),
                       cmap=cmap,
                       linewidth=0)

    ax.plot(output[xkey + '_src'],
            output[xkey + '_halo'][:,-1],
            label="$s_\mathrm{halo} = 1$", marker='o', color='red')

    # plot best fit
    if gta is not None and xkey == 'norm' or xkey == 'dnde':
        s = gta.roi.get_source_by_name(output['src_name'])
        x = s.spectral_pars[output['src_par_name']]['value']
        xerr = s.spectral_pars[output['src_par_name']]['error']

        h = gta.roi.get_source_by_name(output['halo_name'])
        y = h.spectral_pars[output['halo_par_name']]['value']
        yerr = h.spectral_pars[output['halo_par_name']]['error']

        ax.errorbar(x,y,xerr=xerr, yerr=yerr, label='Best fit', marker='o', color='k')

    ax.legend()

    cb = plt.colorbar(im)
    cb.set_label('Delta LogLikelihood')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Source Normalization")
    ax.set_ylabel("Halo Normalization")

    return fig, ax

def set_halo_scale(gta, src_name, halo_source_name, injection_scale_name='Scale'):
    """
    Set the energy scale of the halo to the same one as the central source
    """
    idx_scale = gta.like.par_index(src_name, injection_scale_name) 
    idx_scale_halo = gta.like.par_index(halo_source_name, "Scale") 
    if not gta.like[idx_scale].value() == gta.like[idx_scale_halo].value():
        gta.like[idx_scale_halo].setBounds(gta.like[idx_scale].value() / 10., gta.like[idx_scale].value() * 10.)
        gta.logger.info("setting halo scale energy to {0:.3e}".format(gta.like[idx_scale].value()))
        gta.like[idx_scale_halo] = gta.like[idx_scale].value()

def set_halo_normalization(gta, src_name, halo_source_name,
                           prefactor0, scale0,
                           update_halo=True,
                           injection_norm_name='Prefactor',
                           injection_scale_name='Scale',
                           index=1.):
    """
    Given the normalization of the central source, recompute 
    expected normalization of the halo
    """
    idx_norm = gta.like.par_index(src_name, injection_norm_name) 
    idx_scale = gta.like.par_index(src_name, injection_scale_name) 
    idx_norm_halo = gta.like.par_index(halo_source_name, "Prefactor") 
    idx_scale_halo = gta.like.par_index(halo_source_name, "Scale") 

    new_norm = (gta.like.model[idx_norm].getValue() * gta.like.model[idx_norm].getScale())
    new_norm /= prefactor0 * (scale0 / (gta.like.model[idx_scale].getValue() * gta.like.model[idx_scale].getScale())) ** index

    if update_halo:
        set_halo_scale(gta, src_name, halo_source_name, injection_scale_name=injection_scale_name)
        gta.logger.info("setting halo prefactor to {0:.3e}".format(new_norm))
        gta.set_norm(halo_source_name, new_norm,
                     update_source=False)            
    return new_norm

def build_halo_profile_table(halo_profile_tied, index_par_name='Index', injection_par2_name='Cutoff'):
    """Build a fits table for the logl profile results for the IGMF halo fit"""
    cols = OrderedDict()
    # add B field parameters
    cols = OrderedDict({k: [d['Bfield'][k] for d in halo_profile_tied] for k in ['B', 'maxTurbScale']})
    cols.update(OrderedDict({k: [d['Source'][k] for d in halo_profile_tied] for k in ['th_jet', 'z']}))
    # add spectral parameters
    cols.update(OrderedDict({k: [d[k] if not 'Cutoff' in k else u.Quantity(d[k]).to('TeV').value for d in halo_profile_tied] \
        for k in [index_par_name, injection_par2_name]}))
    # add profile likelhood values
    keys = []
    for k in halo_profile_tied[0].keys():
        if 'name' in k:
            continue
        if 'src' in k or 'halo' in k and not 'name' in k:
            keys.append(k)
    keys.append('loglike')
    keys.append('loglike0')

    cols.update(OrderedDict({k: [d[k] for d in halo_profile_tied] for k in keys}))

    t = Table(cols)
    t['B'].unit='G'
    t['maxTurbScale'].unit='Mpc'
    t['th_jet'].unit='deg'
    t['Cutoff'].unit='TeV'
    t['flux_src'].unit='cm-2 s-1'
    t['flux_halo'].unit='cm-2 s-1'
    t['eflux_halo'].unit='MeV cm-2 s-1'
    t['eflux_src'].unit='MeV cm-2 s-1'
    t['dnde_src'].unit='MeV-1 cm-2 s-1'
    t['dnde_halo'].unit='MeV-1 cm-2 s-1'

    return t


def fit_igmf_halo_scan(gta, modelname,
                       src_name,
                       halo_template_dir,
                       halo_template_suffix='',
                       injection_spectrum='PLSuperExpCutoff',
                       injection_par2_name='Cutoff',
                       injection_norm_name='Prefactor',
                       injection_scale_name='Scale',
                       index_par_name='Index',
                       cov_scale=5.,
                       free_radius_sed=1.,
                       free_bkgs=False,
                       generate_maps=True,
                       generate_seds=False,
                       distance_free_norm=1.,
                       loge_bounds=None,
                       z=None, 
                       ebl_model_name='dominguez',
                       pool_size=1,
                       optimizer=None):
    """
    Fit a cascade halo template on top of a point source.
    The point source is modeled with a given injection spectrum and a loop over the 
    spectrum parameters of this spectrum is performed. In each step of the loop, 
    the correspoding halo template is loaded. 

    The fits template should have contain a META keyword in the header, which in turn
    contains a dictionary with the key 'spectral_parameters'. This key contains a sub dictionary
    with the spectral parameters used to generate the template. The parameter names should match 
    the parameter names of the chosen spectral function. 

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        the gt analysis object

    modelname: str
        name of the model that will be reloaded

    src_name: str
        name of the source for which extension will be tested

    halo_template_dir: str
        path to directory with the halo templates for each tested injection spectrum

    {options}

    halo_template_suffix: str
        an additional file suffix for the halo template.
        Files will be searched for in [halo_template_dir]/[halo_template_suffix]*.fits

    injection_spectrum: str
        name of the injection spectrum of the central point source. 
        Default: PLExpCutoff

    injection_par2_name: str
        name of second spectral parameter of injection spectrum. 
        Default: Cutoff

    injection_norm_name: str
        name of normalization parameter of injection spectrum. 
        Default: Prefactor 

    injection_scale_name: str
        name of energy scale parameter of injection spectrum. 
        Default: Scale

    cov_scale: float
        the cov_scale used for the SED analysis. 
        Default: 5.

    generate_maps: bool
        Generate residual and TS maps

    generate_seds: bool
        Generate seds

    loge_bounds array-like or None
        if not None, use these energy bounds for the analysis. 
        Otherwise, use the energy bounds of the gta object
        Default: None

    free_bkgs: bool
        if True, free normalizations of diffuse backgrounds
        Default: False

    distance_free_norm: float
        distance in deg up to which normalizations are left 
        free in fit. Default: 1.

    free_radius_sed: float
        distance in deg up to which normalizations are left 
        free in SED extraction. Default: 1.

    optimizer: str or None.
        name of optimizer to be used. Default: same optimizer as gta object
    """
    #TODO save and restore fit

    gta.logger.info('Starting IGMF Halo Scan %s'%(modelname))

    if optimizer is None:
        optimizer = gta.config['optimizer']['optimizer']

    # point source models for TS maps
    model_pl20 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    # model for residual map
    model3 = { 'SpatialModel' : 'Gaussian', 'Index' : 2.0, 'SpatialWidth' : 0.1 }

    halo_source_name = 'igmf_halo'

    # don't define any spectral paramameters, 
    # as they will be set to the default values by fermipy which reflect the parameters
    # of the Galactic diffuse model.
    # the spatial file name will be set in each step of the loop
    halo_source_dict =  {'Spatial_Filename' : '',
                         'SpatialModel' : 'MapCubeFunction' }
    
    # get spectral index and normalizaiton of original model
    index = gta.get_src_model(src_name)['spectral_pars'][index_par_name]['scale'] * \
            gta.get_src_model(src_name)['spectral_pars'][index_par_name]['value']
    norm = gta.get_norm(src_name)

    # change the central source spectrum to the desired injection 
    # spectrum
    src = gta.roi.get_source_by_name(src_name)
    if not src['SpectrumType'].split('::')[-1] == injection_spectrum:
        if injection_spectrum == 'PowerLaw':
            gta, index_par_name = set_src_spec_pl(gta, gta.get_source_name(src_name), 
                                                  e0=None,
                                                  e0_free=False)
        elif injection_spectrum == 'PLSuperExpCutoff' or injection_spectrum == 'PLExpCutoff':
            gta, index_par_name = set_src_spec_plexpcut(gta, gta.get_source_name(src_name),
                                                        e0=None,
                                                        e0_free=False,
                                                        index2=1.,
                                                        index2_free=False)
        else:
            raise ValueError("Injection spectrum {0:s} not supported".format(injection_spectrum))
        src = gta.roi.get_source_by_name(src_name)

    if z is not None and 'Ebl' not in src['SpectrumType']:
        gta = add_ebl_atten(gta, src_name, z, eblmodel=ebl_model_name)

    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0], loge_bounds[1])

    # set spectral parameters to original values
    gta.set_norm(src_name, norm)
    gta.set_parameter(src_name, index_par_name, index)

    skydir = gta.roi[src_name].skydir
    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(False, pars=fa.allidx)
    gta.free_sources(skydir=skydir,
                     distance=distance_free_norm,
                     pars=fa.allnorm,
                     exclude=None if free_bkgs else diff_sources)
    # leave index of central source free for first fit
    gta.free_source(gta.get_source_name(src_name), pars=["Index"])

    fit_result_base = gta.fit(optimizer=optimizer)

    gta.write_xml('base')
    gta.write_roi('base')

    if generate_maps:
        gta.tsmap('base', model=model_pl20,
                  #loge_bounds=loge_bounds,
                  make_plots=True)
        gta.residmap('base', model=model3,
                     #loge_bounds=loge_bounds,
                     make_plots=True)

    halo_data = []
    halo_fit_results = []
    halo_profile_tied = []

    gta.logger.info("Base model:")
    gta.print_roi()
    gta.print_model()
    gta.print_params()

    filenames = os.path.join(halo_template_dir, halo_template_suffix + "*.fits")
    files = glob.glob(filenames)
    gta.logger.info("Found {0:n} template files in {1:s}".format(len(files), filenames))

    if not len(files):
        raise ValueError("No template files found")

        # with multiprocessing

    if pool_size > 1:
        with Pool(processes=pool_size) as pool:
            result = pool.starmap(extract_halo_likelihood,
                                 zip(list(idx), repeat(analysis)))

    else:
        # loop over halo template files
        for i, f in enumerate(files):
            o, hd, outprefix, p_name = extract_halo_likelihood(copy.deepcopy(gta.config),
                                                               i,
                                                               files,
                                                               src_name,
                                                               halo_source_name,
                                                               halo_source_dict,
                                                               halo_template_suffix,
                                                               injection_par2_name=injection_par2_name,
                                                               injection_norm_name=injection_norm_name,
                                                               injection_scale_name=injection_scale_name,
                                                               index_par_name=index_par_name,
                                                               cov_scale=cov_scale,
                                                               free_radius_sed=free_radius_sed,
                                                               generate_maps=generate_maps,
                                                               generate_seds=generate_seds,
                                                               optimizer=optimizer)
        halo_profile_tied.append(o)
        halo_fit_results.append(hd)

    np.save(os.path.join(gta.workdir, '{0:s}_data.npy'.format(outprefix)), halo_data)
    np.save(os.path.join(gta.workdir, '{0:s}_fit_results.npy'.format(outprefix)), halo_fit_results)
    np.save(os.path.join(gta.workdir, '{0:s}_profile_tied.npy'.format(outprefix)), halo_profile_tied)

    # build a fits table from halo_profile_tied
    t = build_halo_profile_table(halo_profile_tied,
                                 index_par_name=p_name,
                                 injection_par2_name=injection_par2_name)
    t.write(os.path.join(gta.workdir, '{0:s}_profiles.fits'.format(outprefix)), overwrite=True)

    gta.logger.info('Finished IGMF Halo Scan {0:s}'.format(modelname))
    return halo_profile_tied


def extract_halo_likelihood(config,
                            i,
                            halo_template_files,
                            src_name,
                            halo_source_name,
                            halo_source_dict,
                            halo_template_suffix,
                            injection_par2_name='Cutoff',
                            injection_norm_name='Prefactor',
                            injection_scale_name='Scale',
                            index_par_name='Index',
                            cov_scale=5.,
                            free_radius_sed=1.,
                            generate_maps=True,
                            generate_seds=False,
                            model_pl20={'SpatialModel': 'PointSource', 'Index': 2.0 },
                            model3={'SpatialModel': 'Gaussian', 'Index': 2.0, 'SpatialWidth': 0.1 },
                            optimizer=None):
    """Extract halo likelihood for one particular halo file"""
    # get the used spectral parameters
    template = fits.open(halo_template_files[i])
    spectral_parameters = yaml.safe_load(template[0].header["META"])['spectral_parameters']
    # get the simulation parameters
    sim_parameters = yaml.safe_load(template[0].header["META"])['sim_config']
    # spectral index
    if index_par_name == 'Index1' and index_par_name not in spectral_parameters.keys() \
            and "Index" in spectral_parameters.keys():
        p = spectral_parameters["Index"]
        p_name = "Index"
    else:
        p = spectral_parameters[index_par_name]
        p_name = index_par_name

    # second spectral parameter
    p2 = spectral_parameters[injection_par2_name]
    if injection_par2_name == 'Cutoff':
        p2 = u.Quantity(p2).to("MeV").value

    # used normalization and scale in halo template simulation
    prefactor0 = u.Quantity(spectral_parameters[injection_norm_name]).to("s-1 cm-2 MeV-1").value
    scale0 = u.Quantity(spectral_parameters[injection_scale_name]).to("MeV").value
    outprefix = "{0:s}".format(halo_template_suffix)

    gta = GTAnalysis(config, logging={'verbosity': 3})
    gta.logger.info(
        'Fitting central source in Halo with index {0:.3f} and {1:s} {2:.3f}'.format(p, injection_par2_name, p2))

    halo_source_dict['Spatial_Filename'] = os.path.join(halo_template_files[i])
    gta.logger.info("Using Spatial File {0:s}".format(halo_source_dict['Spatial_Filename']))
    # reload base model without halo
    gta.load_xml('base')
    # set the spectral parameters of the central source
    gta.set_parameter(src_name, injection_par2_name, p2,
                      update_source=False)
    gta.set_parameter(src_name, index_par_name, -1. * p,
                      update_source=False)
    # add the halo
    try:
        gta.add_source(halo_source_name, halo_source_dict, free=False)
    except RuntimeError as e:
        gta.logger.error("Couldn't add halo source")
        gta.logger.error("Error was {0}".format(e))
        gta.logger.error("Stopping the loop and saving output")
    #    break
    gta.free_norm(halo_source_name)
    new_norm = set_halo_normalization(gta, src_name, halo_source_name,
                                      prefactor0, scale0,
                                      update_halo=True,
                                      injection_norm_name=injection_norm_name,
                                      injection_scale_name=injection_scale_name,
                                      index=p)
    logging.debug("npred halo before fit: {0}".format(gta.model_counts_spectrum(halo_source_name,
                                                                                gta.loge_bounds[0],
                                                                                gta.loge_bounds[1])))
    # perform the fit
    # treating the halo and source spectral normalizations as independent parameters
    fit_result = gta.fit(optimizer=optimizer)

    # TODO: check if fit should be performed like this:
    # gta.fit(update=False, optimizer=optimizer)
    gta.print_params(loglevel=logging.INFO)
    gta.print_model(loglevel=logging.INFO)
    # gta.update_source(halo_source_name,reoptimize=True,
    #                  optimizer={'optimizer' : optimizer})
    # gta.update_source(src_name,reoptimize=True,
    #                  optimizer={'optimizer' : optimizer})
    halo_modelname = '{0:s}-{1:03n}id'.format(outprefix, i)
    halo_modelname = halo_modelname.replace('.', 'p')
    gta.write_roi(halo_modelname,
                  make_plots=False)
    if generate_maps:
        gta.tsmap(halo_modelname, model=model_pl20,
                  # loge_bounds=loge_bounds,
                  make_plots=True)
        gta.residmap(halo_modelname, model=model3,
                     # loge_bounds=loge_bounds,
                     make_plots=True)
    # perform 2D likelihood fit over source and halo norm
    o = profile_halo_src_norms(gta, src_name, halo_source_name,
                               scale0, prefactor0,
                               index=p,
                               src_norm_name=injection_norm_name,
                               halo_norm_name='Prefactor',
                               src_scale_name=injection_scale_name,
                               reoptimize=True,
                               savestate=True,
                               sigma=3.,
                               xsteps=30,
                               ysteps=31)
    spectral_parameters.update({'Spatial_Filename': halo_template_files[i]})
    o.update(spectral_parameters)
    o.update(sim_parameters)
    # make sure to reload fit
    gta.load_xml(halo_modelname)
    # compute the SEDs
    if generate_seds:
        gta.sed(src_name,
                prefix='src_{0:s}_cov'.format(halo_modelname),
                outfile='src_{0:s}_cov'.format(outprefix),
                free_radius=free_radius_sed,
                cov_scale=cov_scale,
                optimizer={'optimizer': 'MINUIT'},
                make_plots=False)

        gta.sed(halo_source_name,
                prefix='halo_{0:s}_cov'.format(halo_modelname),
                outfile='halo_{0:s}_cov'.format(outprefix),
                free_radius=None,  # free_radius_sed, # does not work for diffuse sources
                cov_scale=cov_scale,
                optimizer={'optimizer': 'MINUIT'},
                make_plots=False)
    ul_flux = get_parameter_limits(gta.roi[halo_source_name]['flux_scan'],
                                   gta.roi[halo_source_name]['loglike_scan'])
    ul_eflux = get_parameter_limits(gta.roi[halo_source_name]['eflux_scan'],
                                    gta.roi[halo_source_name]['loglike_scan'])
    gta.roi[halo_source_name]['flux_err'] = ul_flux['err']
    gta.roi[halo_source_name]['eflux_err'] = ul_eflux['err']
    gta.logger.info('%{0:s} Injection Index: {1:6.3f} Injection {5:s}: {2:6.2f} TS: {3:6.2f} Flux: {4:8.4g}'.format(
        halo_modelname, p, p2,
        gta.roi[halo_source_name]['ts'],
        gta.roi[halo_source_name]['flux'],
        injection_par2_name
        )
    )

    halo_data = copy.deepcopy(gta.roi[halo_source_name].data)
    gta.delete_source(halo_source_name, save_template=False)
    return o, halo_data, outprefix, p_name
