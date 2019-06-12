import numpy as np
from glob import glob
from astropy.io import fits
import logging
import yaml
from astropy.table import Table
from os import path
from fermipy.utils import met_to_mjd
from scipy.stats import linregress
from LikelihoodState import LikelihoodState
import fermiAnalysis as fa 
import subprocess
import copy
from fermiAnalysis import adaptivebinning as ab

def calculate_exposure(srcconf,earray):
    """
    Calculate the Fermi LAT exposure

    :param srcconf:
    dict, source configuration
    :return:
    tuple with array-like objects for time, front, and back exposure
    """
    texp_file = path.join(srcconf['fileio']['outdir'],
                          "exposure_{0[emin]:n}-{0[emax]:n}MeV.npz".format(srcconf['selection']))

    if path.isfile(texp_file):

        npzfile = np.load(texp_file)
        texp, front, back = npzfile['texp'], npzfile['front'], npzfile['back']
        logging.info("loaded exposure from {0:s}".format(texp_file))

    else:
        # get the gta object
        config = myconf2fermipy(srcconf)
        config['fileio']['scratchdir'] = None
        config['fileio']['usescratch'] = False
        gta = GTAnalysis(config, logging={'verbosity': 3})

        front, back = [], []
        for energy in earray:
            texp, f, b = ab.comp_exposure_phi(gta, energy=energy)
            front.append(f)
            back.append(f)

        front = np.array(front)
        back = np.array(back)
        np.savez(texp_file, texp=texp, front=front,
                 back=back, earray=earray)
        logging.info("Saved exposure to {0:s}".format(texp_file))

    return texp, front, back

def excise_solar_flares_grbs(tsmin = 100.,
    sf_file = "/u/gl/omodei/SunMonitor_P8/SUN-ORB-TEST/SUN-ORB-TEST_jean.txt",
    grb_file = "/u/gl/omodei/GBMTRIGCAT-v2/DATA/LAT2CATALOG-v4-LTF/LAT2CATALOG-v4-LTF_jean.txt",
    tmax = 915148805.,outdir = "./"):
    """
    Derive a filter expression to exclude bright solar flares and GRBs

    Parameters
    ----------
    tsmin: float
        minimum ts for which solar flare or grb will be excluded (default: 100.)
    sf_file: str
        path to file with solar flare times and ts values
    grb_file: str
        path to file with solar flare times and ts values
    tmax: float
        maximum mission lifetime in MET (default: 01-01-2030)
    outidr: str
        path where output files are stored (default: ./)

    Returns
    -------
    Path to fits file containing the GTIs
    """
    sf = np.loadtxt(sf_file, usecols = (3,5,6), delimiter = ',')
    grb = np.loadtxt(grb_file, usecols = (3,5,6), delimiter = ',')
    ms = sf[:,-1] >= tsmin # cut on solar flares
    mg = grb[:,-1] >= tsmin # cut on grbs
    #filter_str = "!((START >= {start:.2f} && STOP < {stop:.2f})" \
    #             "||(START < {start:.2f} && STOP > {start:.2f} && STOP < {stop:.2f})" \
    #             "||(START < {stop:.2f} && STOP > {stop:.2f} && START >= {start:.2f}))"
    #out = ""
    #for s in sf[ms]:
        #out += filter_str.format(start = s[0],stop = s[0] + s[1])
        #out += "&&"
    #for i,g in enumerate(grb[mg]):
    #    out += filter_str.format(start = g[0],stop = g[0] + g[1])
    #    if i < np.sum(mg) - 1:
    #        out += "&&"
    #sum_excl = sf[ms][:,1].sum() + grb[mg][:,1].sum()
    #logging.info("Excluding {0:.2f} days".format(sum_excl / 3600. / 24.))

    # concetanate the two 
    btis = np.vstack([sf[ms], grb[mg]])
    # sort them
    btis = btis[np.argsort(btis[:,0])]

    # remove overlaping intervals by combining them
    bti_start = []
    bti_stop = []

    bti_start = btis[:,0]
    bti_stop = btis[:,0] + btis[:,1]

    stop,start = bti_stop[:-1] , bti_start[1:]
    m = np.where(stop > start)
    while len(m[0]):
        bti_stop = np.concatenate([np.delete(bti_stop[:-1], m[0][0]), [bti_stop[-1]]])
        bti_start = np.concatenate([[bti_start[0]],np.delete(bti_start[1:], m[0][0])])
        stop,start = bti_stop[:-1] , bti_start[1:]
        m = np.where(stop > start)

    gti_start = np.concatenate([[0.], bti_stop])
    gti_stop = np.concatenate([bti_start , [tmax]])
    gtis = np.vstack([gti_start, gti_stop])

    header ="""START D\nSTOP D"""
    np.savetxt(path.join(outdir,"gtis.txt"),gtis.T)
    with open(path.join(outdir, "header.txt"),'w') as f:
        f.write(header)
        f.close()
    subprocess.call(['ftcreate', path.join(outdir, "header.txt"),
            path.join(outdir,"gtis.txt"), path.join(outdir, 'nosolarflares_nogrbs.gti'),
            "clobber=yes","extname=GTI"])
    
    return path.join(outdir, 'nosolarflares_nogrbs.gti')

def mjd_to_met(time):
    """Convert MJD to Fermi MET"""
    return (time - 54682.65) * (86400.) + 239557414.0
    #""""Convert mission elapsed time to mean julian date."""
    #    return 54682.65 + (time - 239557414.0) / (86400.)

def fit_with_retries(gta, fit_config, target,
                        alt_spec_pars = None,
                        tsmax = 1e3):
    """Fit the model and retry if not converged"""

    try:
        o = gta.optimize() # perform an initial fit
        logging.debug(o)
        gta.print_roi()
    except RuntimeError as e:
        logging.warning("optimize failed with {0}. Trying alternative parameters".format(e))
        if not type(alt_spec_pars) == type(None):
            gta.set_source_spectrum(target,
                spectrum_type = gta.roi.get_source_by_name(target)['SpectrumType'],
                spectrum_pars= alt_spec_pars.spectral_pars)

    try:
        f = gta.fit()
    except RuntimeError as e:
        logging.warning("fit failed with {0}. Trying to alternative parameters".format(e))
        if not type(alt_spec_pars) == type(None):
            gta.set_source_spectrum(target,
                spectrum_type = gta.roi.get_source_by_name(target)['SpectrumType'],
                spectrum_pars= alt_spec_pars.spectral_pars)


    nfree0 = print_free_sources(gta)

    retries = f['config']['retries']
    tsfix = fit_config['ts_fixed']
    while not f['fit_success'] and retries > 0:
        # Fix more sources
        logging.warning("fit did not suceed, failed with status {0}, fixing sources".format(f['fit_success']))
        tsfix *= 3
        gta.free_sources(minmax_ts=[None,tsfix],free=False,
                        pars=fa.allidx + fa.allnorm)
        logging.info("retrying fit, {0:n} retries left".format(retries))
        nfree = print_free_sources(gta)

        if nfree == nfree0 and tsfix < tsmax:
            continue
        elif nfree == nfree0 and tsfix >= tsmax:
            # freeze all index parameters
            logging.info("TS fix maximum reached, fixing all indeces")
            gta.free_sources(free=False,
                                pars=fa.allidx)
            nfree0 = print_free_sources(gta)
        else:
            nfree0 = nfree
            logging.info("Retrying fit with {0:n} free parameters".format(nfree))
        try:
            o = gta.optimize() # perform an initial fit
            logging.debug(o)
            gta.print_roi()
        except RuntimeError as e:
            logging.warning("optimize failed with {0}. Trying to continue anyway".format(e))

        f = gta.fit()
        retries -= 1

    if not f['fit_success']:
        logging.warning("fit did not suceed, fixing beta and Index2")
        gta.free_source(target, pars = ['beta', 'Index2'], free = False)
        nfree = print_free_sources(gta)
        logging.warning('Now there are {0:n} free parameters'.format(nfree))
        f = gta.fit()

    return f,gta

def set_free_pars_avg(gta, fit_config, freezesupexp = False):
    """
    Freeze are thaw fit parameters for average fit
    """
    # run the fitting of the entire time and energy range
    try:
        o = gta.optimize() # perform an initial fit
        logging.debug(o)
        gta.print_roi()
    except RuntimeError as e:
        logging.warning("optimize failed with {0}. Trying to continue anyway".format(e))

    # Free all parameters of all Sources within X deg of ROI center
    #gta.free_sources(distance=fit_config['ps_dist_all'])
    # Free Normalization of all Sources within X deg of ROI center
    gta.free_sources(distance=fit_config['ps_dist_norm'],pars=fa.allnorm)
    # Free spectra parameters of all Sources within X deg of ROI center
    gta.free_sources(distance=fit_config['ps_dist_idx'],pars=fa.allidx)
    # Free all parameters of isotropic and galactic diffuse components
    gta.free_source('galdiff', pars='norm', free = fit_config['gal_norm_free'])
    gta.free_source('galdiff', pars=['index'], free = fit_config['gal_idx_free'])
    gta.free_source('isodiff', pars='norm', free = fit_config['iso_norm_free'])
    # Free sources with TS > X
    gta.free_sources(minmax_ts=[fit_config['ts_norm'],None], pars =fa.allnorm)
    # Fix sources with TS < Y
    gta.free_sources(minmax_ts=[None,fit_config['ts_fixed']],free=False,
            pars = fa.allnorm + fa.allidx)
    # Fix sources Npred < Z
    gta.free_sources(minmax_npred=[None,fit_config['npred_fixed']],free=False,
            pars = fa.allnorm + fa.allidx)
    if freezesupexp:
        logging.info("Freezing Index2 for all sources")
        gta.free_sources(pars = ['Index2'], free = False)
    gta.print_roi()
    print_free_sources(gta)
    return gta

def set_lc_bin(tmin, tmax, dt,n, ft1 = 'None'):
    """
    Set the starting and ending times of a light curve bin
    given the min and maximum times of the full time interval

    Parameters
    ----------
    tmin: float
        start time of full interval in seconds
    tmax: float
        end time of full interval in seconds
    dt: int, float or str
        length of time interval in seconds, if str, use GTI as time intervals
    n: int
        number of time interval, starts at 0

    kwargs
    ------
    ft1: str
        if type of dt is string, ft1 is path to ft1 that contains a Table with name GTI

    Returns
    -------
    tuple containing start and end time of n-th interval as well as total number of intervals
    """
    if type(dt) == float or type(dt) == int:
        nint = int(np.floor((float(tmax) - float(tmin)) / float(dt)))
        tmin_new = tmin + n * dt
        tmax_new = tmin + (n + 1) * dt
    elif type(dt) == str:
        logging.debug(ft1)
        t = Table.read(ft1, hdu = 'GTI')
        # determine minimum time
        if tmin < t['START'].min():
            min_id = 0

        elif len(np.where((tmin <= t['STOP']) & (tmin >= t['START']))[0]):
            min_id = np.where((tmin <= t['STOP']) & (tmin >= t['START']))[0][0]
        else:
            min_id = np.where(t0 < t['START'])[0]
            if not len(min_id):
                raise ValueError("Minimum time outside GTI window")
            else:
                min_id = min_id[0]

        # determine minimum time
        if tmax >= t['STOP'].max():
            max_id = t['STOP'].size 
        elif len(np.where((tmax <= t['STOP']) & (tmax >= t['START']))[0]):
            max_id = np.where((tmax <= t['STOP']) & (tmax >= t['START']))[0][0]
        else:
            max_id = np.where(tmax >= t['STOP'])[0]
            if not len(max_id):
                raise ValueError("Maximum time outside GTI window")
            else:
                max_id = max_id[-1]

        tmin_new = (t['START'][min_id:max_id + 1])[n]
        tmax_new = (t['STOP'][min_id:max_id + 1])[n]
        nint = max_id - min_id
    return tmin_new, tmax_new, int(nint)

def print_free_sources(gta):
    """Print the free source parameters"""
    logging.info('free source parameters:')
    nfree = 0
    for s in  gta.get_sources() :                                                           
        for k in s.spectral_pars.keys():
            if s.spectral_pars[k]['free']:   
                logging.info('{0:s}: {1:s}'.format(s.name, k))
                nfree += 1
    return nfree

def refit(gta, src, f, tsfix):
    """
    If fit was not successful, refit while freezing more sources

    Parameters
    ----------
    gta: gtanalysis object
    src: str, 
        name of central source
    f: dict
        output of first gta.fit
    ts_fixed: float
       initial ts below which sources are frozen 

    Returns
    -------
    tuple with gta and final fit dictionary
    """
    retries = f['config']['retries']
    while not f['fit_success'] and retries > 0:
        gta.free_source(src, pars = ['beta', 'Index2'], free = False)
        # Fix more sources
        tsfix *= 3
        gta.free_sources(minmax_ts=[None,tsfix],free=False,
                        pars=fa.allidx + fa.allnorm)
        logging.info("retrying fit")
        for s in  gta.get_sources() :                                                           
            for k in s.spectral_pars.keys():
                if s.spectral_pars[k]['free']:   
                    logging.info('{0:s}: {1:s}'.format(s.name, k))
        o = gta.optimize()
        gta.print_roi()
        f = gta.fit()
        retries -= 1
    return gta, f



ebl_st_model_list = ['kneiske', 'primack05', 'kneiske_highUV', 'stecker05',
                       'franceschini', 'finke', 'gilmore', 'stecker05_FE',
                       'salamonstecker', 'generic', 'dominguez', 'kneiskedole',
                       'kneiskedole10_CMB', 'gilmore09', 'gilmore_fiducial', 
                       'gilmore_fixed', 'scully14_lowOp', 'scully14_highOp',
                       'inoue', 'helgasonKashlinsky']

def set_src_spec_pl(gta, src, e0 = None, e0_free = False):
    """
    Set an arbitrary source spectrum to a power law
    with a fitted index (using scipy.stats.linegress)
    and a normalization that gives the same flux 
    as the original spectrum in the considered energy range

    Parameters
    ----------
    gta: `~fermipy.gtanalysis`
        fermipy analysis object

    src: string
        source name

    e0: float or None (optional)
        if float, force pivot energy to this value

    e0_free: bool
        let pivot energy free during fit 
    """
    m = gta.roi.get_source_by_name(src)
    loge = np.sort(gta.log_energies)
    logecen = 0.5 * (loge[1:] + loge[:-1])
    # get the source flux
    flux = []
    for i,le in enumerate(loge[:-1]):
        flux.append(gta.like[src].flux(10.**le,10.**loge[i+1]))
    # get the power-law index
    flux = np.array(flux)
    flux[flux <= 0.] = 1e-40 * np.ones(np.sum(flux <= 0.))
    r = linregress(logecen,np.log10(flux))
    if -1. * r.slope < 0. or -1. * r.slope > 5.:
        slope = -2.
    else:
        slope = r.slope

    if e0 is None:
        if m['SpectrumType'] == 'LogParabola':
            e0 = m['spectral_pars']['Eb']['value']
        else:
            try:
                e0 = m['spectral_pars']['Scale']['value']
            except KeyError:
                logging.warning("Did not find 'Scale' in spectral pars, using pivot energy for E0")
                e0 = gta.get_src_model(src)['pivot_energy']

    prefactor = gta.like[src].flux(10.**loge[0],10.**loge[-1])
    # compute the normalization so that integrated 
    # flux is the same
    if slope == -1.:
        prefactor /= (e0 * (loge[-1] - loge[0]))
        prefactor *= np.log(10.)
    else:
        gp1 = 1. + slope
        prefactor *=  gp1 * e0**slope
        prefactor /= (10.**(gp1 * loge[-1]) - 10.**(gp1 * loge[0]))
    # change the source spectrum
    pars = {}
    pars['Index'] =  dict(error= np.nan, 
                            free = True, 
                            max= 5.0, min= 0.0, 
                            name= 'Index', scale= -1.0,
                            value= slope * -1.)
    pars['Prefactor'] = dict(error= np.nan,
                            free= True,
                            max= prefactor * 1e5, min= prefactor * 1.0e-05,
                            name= 'Prefactor',scale= 1.0,
                            value= prefactor)
    pars['Scale'] = dict(error= np.nan,
                        free= e0_free,
                        max= e0 * 10., min= e0 / 10.,
                        name= 'Scale', scale= 1.0,
                        value= e0)

    gta.set_source_spectrum(src,spectrum_type = 'PowerLaw', spectrum_pars= pars)
    logging.info('Changed spectrum of {0:s} to power law'.format(src))
    logging.info('New source parameters are: {0}'.format(gta.roi.get_source_by_name(src).spectral_pars))
    return gta

def set_src_spec_plexpcut(gta, src, e0 = None, e0_free = False):
    """
    Set an arbitrary source spectrum to a power law
    with exponential cut off 

    Parameters
    ----------
    gta: `~fermipy.gtanalysis`
        fermipy analysis object

    src: string
        source name

    e0: float or None (optional)
        if float, force pivot energy to this value

    e0_free: bool
        let pivot energy free during fit 
    """
    m = gta.roi.get_source_by_name(src)

    if e0 is None:
        if m['SpectrumType'] == 'LogParabola':
            e0 = m['spectral_pars']['Eb']['value']
        else:
            try:
                e0 = m['spectral_pars']['Scale']['value']
            except KeyError:
                logging.warning("Did not find 'Scale' in spectral pars, using pivot energy for E0")
                e0 = gta.get_src_model(src)['pivot_energy']

    loge = np.sort(gta.log_energies)
    logecen = 0.5 * (loge[1:] + loge[:-1])
    # get the source flux
    flux = []
    for i,le in enumerate(loge[:-1]):
        flux.append(gta.like[src].flux(10.**le,10.**loge[i+1]))
    # get the power-law index
    flux = np.array(flux)
    flux[flux <= 0.] = 1e-40 * np.ones(np.sum(flux <= 0.))
    maskcen = 10.**logecen <= e0
    mask = 10.**loge <= e0

    r = linregress(logecen[maskcen],np.log10(flux[maskcen]))
    if -1. * r.slope < 0. or -1. * r.slope > 5.:
        slope = -2.
    else:
        slope = r.slope

    prefactor = gta.like[src].flux(10.**loge[mask].min(),10.**loge[mask].max())
    # compute the normalization so that integrated 
    # flux is the same
    if slope == -1.:
        prefactor /= (e0 * (loge[mask].max() - loge[mask].min()))
        prefactor *= np.log(10.)
    else:
        gp1 = 1. + slope
        prefactor *=  gp1 * e0**slope
        prefactor /= (10.**(gp1 * loge[mask].max()) - 10.**(gp1 * loge[mask].min()))

    prefactor_scale = np.floor(np.log10(prefactor))
    # change the source spectrum
    pars = {}
    pars['Index1'] =  dict(error= np.nan, 
                            free = True, 
                            max= 5.0, min= 0.0, 
                            name= 'Index', scale= -1.0,
                            value= slope * -1.)
    pars['Prefactor'] = dict(error= np.nan,
                            free= True,
                            max= prefactor * 1e5 / 10.**prefactor_scale,
                            min= prefactor * 1.0e-05 / 10.**prefactor_scale,
                            name= 'Prefactor',scale= 10.**prefactor_scale,
                            value= prefactor / 10.**prefactor_scale)
    pars['Scale'] = dict(error= np.nan,
                        free= e0_free,
                        max= e0 * 10., min= e0 / 10.,
                        name= 'Scale', scale= 1.0,
                        value= e0)
    pars['Cutoff'] = dict(error= np.nan,
                        free= False,
                        max= 1e+5, min= 1.0,
                        name= 'Cutoff', scale= 1.0,
                        value= e0 * 2.)
    pars['Index2'] =  dict(error= np.nan, 
                            free = True, 
                            max= 2.0, min= 0.0, 
                            name= 'Index2', scale= 1.0,
                            value= 1.)

    gta.set_source_spectrum(src,spectrum_type = 'PLSuperExpCutoff', spectrum_pars= pars)
    logging.info('Changed spectrum of {0:s} to PLSuperExpCutoff'.format(src))
    logging.info('New source parameters are: {0}'.format(gta.roi.get_source_by_name(src).spectral_pars))
    return gta

def add_ebl_atten(gta, src, z, eblmodel = 'dominguez', force_lp= False):
    """
    Change the source model of a source to include EBL attenuation

    Parameters
    ----------
    gta: `~fermipy.gtanalysis`
        fermipy analysis object

    src: string
        source name

    z: float
        source redshift

    kwargs
    ------
    eblmodel: string
        ebl model identifier, default : dominguez

    force_lp: bool
        if true, force the intrinsic spectrum to be a log parabola, default: False

    Returns
    -------
    modified `~fermipy.gtanalysis` object
    """
    m = gta.roi.get_source_by_name(src)
    if not eblmodel in ebl_st_model_list:
        raise ValueError("Unknown EBL model chosen. Availbale models are {0}".format(ebl_st_model_list))

    if force_lp: # bug in src maps for latter model
        new_spec = 'EblAtten::{0:s}'.format('LogParabola')
    else:
        new_spec = 'EblAtten::{0:s}'.format(m['SpectrumType'])
        new_spec_dict = m.spectral_pars
    if m['SpectrumType'] == 'PowerLaw': # only PL2 exists with EBL attenuation in glorious ST
        new_spec += '2'
        new_spec_dict['LowerLimit'] = dict(value = 10.**gta.log_energies[0])
        new_spec_dict['UpperLimit'] = dict(value = 10.**gta.log_energies[-1])
        val = gta.like[gta.get_source_name(src)].flux(
                                        10.**gta.log_energies[0],10.**gta.log_energies[-1])
        new_spec_dict['Integral'] = dict(free = True, value = val, min = val *1e-5, max = val * 1e5)
        for p in ['Scale','Prefactor']:
            new_spec_dict.pop(p,None)

    logging.info("Using EBL model {0:s} with ST id {1:n}".format(
        eblmodel, ebl_st_model_list.index(eblmodel)))
    new_spec_dict['ebl_model'] = dict(
            error = np.nan, 
            free =  False, 
            min = 0, max = len(ebl_st_model_list),
            name = 'ebl_model', scale = 1.0,
            value= ebl_st_model_list.index(eblmodel))
    new_spec_dict['redshift'] = dict(
            error = np.nan, 
            free =  False, 
            max = 2., min = 0.,
            name = 'redshift', scale = 1.0,
            value= z)
    new_spec_dict['tau_norm'] = dict(
            error = np.nan, 
            free =  False, 
            max = 2., min = 0.,
            name = 'tau_norm', scale = 1.0,
            value= 1.)
    logging.info('Changing spectrum of {0:s} to {1:s}...'.format(src, new_spec))
    logging.info('Dictionary: {0}'.format(new_spec_dict))
    gta.set_source_spectrum(src,spectrum_type=new_spec, spectrum_pars= new_spec_dict)
    logging.info('Changed spectrum of {0:s} to {1:s}'.format(src, new_spec))
    logging.info('New source parameters are: {0}'.format(gta.roi.get_source_by_name(src).spectral_pars))
    return gta

def add_columns(columns, fitsfile, hdu, createsedlc = False):
    with fits.open(fitsfile) as f:

        if hdu == 'CATALOG':
            idx = np.argmin(f[hdu].data['offset']) # index of central source
            s = slice(idx,idx+1)

            if createsedlc:
                fsed = glob(path.join(path.dirname(fi), 'lc_sed*.fits'))
                if not len(fsed):
                    logging.warning("No SED file in {0:s}".format(path.join(path.dirname(fi), 'lc_sed*.fits')))
                else:
                    tsed = Table.read(fsed[0])
                    columns['emin_sed'].append(tsed['e_min'].data)
                    columns['emax_sed'].append(tsed['e_max'].data)
                    columns['eref_sed'].append(tsed['e_ref'].data)
                    columns['norm'].append(tsed['norm'].data)
                    columns['norm_err'].append(tsed['norm_err'].data)
                    columns['norm_errp'].append(tsed['norm_errp'].data)
                    columns['norm_errn'].append(tsed['norm_errn'].data)
                    columns['norm_ul'].append(tsed['norm_ul'].data)
                    columns['ts_sed'].append(tsed['ts'].data)
                    columns['ref_dnde'].append(tsed['ref_dnde'].data)
                    columns['ref_flux'].append(tsed['ref_flux'].data)
                    columns['ref_eflux'].append(tsed['ref_eflux'].data)
                    columns['ref_npred'].append(tsed['ref_npred'].data)
                    columns['norm_scan'].append(tsed['norm_scan'].data)
                    columns['dloglike_scan'].append(tsed['dloglike_scan'].data)
                    del tsed


            c = yaml.load(f[0].header['CONFIG'])
            columns['tmin'].append(c['selection']['tmin'])
            columns['tmax'].append(c['selection']['tmax'])
            columns['emin'].append(c['selection']['emin'])
            columns['emax'].append(c['selection']['emax'])
            columns['param_names'].append(np.squeeze(f[hdu].data['param_names'][s]))
            columns['param_values'].append(np.squeeze(f[hdu].data['param_values'][s]))
            columns['param_errors'].append(np.squeeze(f[hdu].data['param_errors'][s]))
            columns['npred'].append(np.squeeze(f[hdu].data['npred'][s]))
            columns['ts'].append(np.squeeze(f[hdu].data['npred'][s]))
            columns['eflux_err'].append(np.squeeze(f[hdu].data['eflux_err'][s]))
            columns['eflux'].append(np.squeeze(f[hdu].data['eflux'][s]))
            columns['flux_err'].append(np.squeeze(f[hdu].data['flux_err'][s]))
            columns['flux'].append(np.squeeze(f[hdu].data['flux'][s]))

        elif hdu == 'LIGHTCURVE':
            s = slice(f[hdu].data.size)
            columns['tmin'] = np.concatenate([columns['tmin'], f[hdu].data['tmin'][s]])
            columns['tmax'] = np.concatenate([columns['tmax'], f[hdu].data['tmax'][s]])
            columns['emin'] = np.concatenate([columns['emin'], np.zeros_like(f[hdu].data['tmin'][s])])
            columns['emax'] = np.concatenate([columns['emax'], np.zeros_like(f[hdu].data['tmin'][s])])
            columns['npred'] = np.concatenate([columns['npred'], f[hdu].data['npred'][s]])
            columns['ts'] = np.concatenate([columns['ts'], f[hdu].data['ts'][s]])
            columns['eflux_err'] = np.concatenate([columns['eflux_err'], f[hdu].data['eflux_err'][s]])
            columns['eflux'] = np.concatenate([columns['eflux'], f[hdu].data['eflux'][s]])
            columns['flux_err'] = np.concatenate([columns['flux_err'], f[hdu].data['flux_err'][s]])
            columns['flux'] = np.concatenate([columns['flux'], f[hdu].data['flux'][s]])
            columns['param_names'] = np.vstack([columns['param_names'],f[hdu].data['param_names'][s]])
            columns['param_values'] = np.vstack([columns['param_values'],f[hdu].data['param_values'][s]])
            columns['param_errors'] = np.vstack([columns['param_errors'],f[hdu].data['param_errors'][s]])
            if 'flux_fixed' in f[hdu].columns.names:
                columns['eflux_err_fixed'] = np.concatenate([columns['eflux_err_fixed'], f[hdu].data['eflux_err_fixed'][s]])
                columns['eflux_fixed'] = np.concatenate([columns['eflux_fixed'], f[hdu].data['eflux_fixed'][s]])
                columns['flux_err_fixed'] = np.concatenate([columns['flux_err_fixed'], f[hdu].data['flux_err_fixed'][s]])
                columns['flux_fixed'] = np.concatenate([columns['flux_fixed'], f[hdu].data['flux_fixed'][s]])
                columns['ts_fixed'] = np.concatenate([columns['ts_fixed'], f[hdu].data['ts_fixed'][s]])

        for h in f:
            del h.data
        del f

    return columns

def collect_lc_results(outfiles, hdu = "CATALOG",
                                createsedlc = False):

    ff = glob(outfiles)
    ff = sorted(ff, key = lambda f: int(path.dirname(f).split('/')[-1]))
    logging.debug('collecting results from {0:n} files'.format(len(ff)))

    columns = dict(tmin = [], tmax = [],
                    emin = [], emax = [],
                    npred = [], ts = [], eflux = [],
                    eflux_err = [],flux = [],flux_err = [],
                    param_names = [], param_values = [], param_errors = [])
    if createsedlc:
        columns.update(dict(
                    emin_sed = [],
                    emax_sed = [],
                    eref_sed = [],
                    norm = [],
                    norm_err = [],
                    norm_errp = [],
                    norm_errn = [],
                    norm_ul = [],
                    ts_sed = [],
                    ref_dnde = [],
                    ref_flux = [],
                    ref_eflux = [],
                    ref_npred = [],
                    norm_scan = [],
                    dloglike_scan = []))

    with fits.open(ff[0]) as f:
        if 'flux_fixed' in f[hdu].columns.names:
            columns.update(dict(ts_fixed = [],
            flux_fixed = [],
            flux_err_fixed  = [],
            eflux_err_fixed  = [],
            eflux_fixed = []
            ))
        del f

    for i,fi in enumerate(ff):

        columns = add_columns(columns, fi, hdu = hdu, createsedlc = createsedlc)
        logging.debug('Opening {0:s}'.format(fi))

    columns['tmin'] = np.squeeze(met_to_mjd(np.array(columns['tmin'])))
    columns['tmax'] = np.squeeze(met_to_mjd(np.array(columns['tmax'])))

    t = Table(columns)
    t['tmin'].unit = 'MJD'
    t['tmax'].unit = 'MJD'
    t['emin'].unit = 'MeV'
    t['emax'].unit = 'MeV'
    t['eflux_err'].unit = 'MeV'
    t['flux'].unit = 'cm-2 s-1'
    t['flux_err'].unit = 'cm-2 s-1'
    t['eflux'].unit = 'MeV cm-2 s-1' 
    t['eflux_err'].unit = 'MeV cm-2 s-1' 

    return t

def myconf2fermipy(myconfig):
    """Convert my config files to fermipy config files"""
    config = copy.deepcopy(myconfig)
    if config['data']['ltcube'] == '':
        config['data'].pop('ltcube',None)
    for k in ['configname', 'tmp', 'log', 'fit_pars']: 
        config.pop(k,None)
    if 'adaptive' in config['lightcurve'].keys():
        config['lightcurve'].pop('adaptive',None)
    if type(config['lightcurve']['binsz']) == str:
        if config['lightcurve']['binsz'] == 'gti':
            config['lightcurve']['binsz'] = 3600. * 3.
        elif len(config['lightcurve']['binsz'].strip('gti')):
            if 'min' in config['lightcurve']['binsz'] > 0:
                config['lightcurve']['binsz'] = float(
                    config['lightcurve']['binsz'].strip('gti').strip('min')) * 60.
    return config

def read_srcprob(ft1srcprob, srcname):
    """Get the srcprob from an FT1 file or a list of files and one source"""
    if type(ft1srcprob) == str:
        ft1srcprob = [ft1srcprob]

    if not len(ft1srcprob):
        raise ValueError("No files found in {0}".format(ft1srcprob))
    for io,o in enumerate(ft1srcprob):
        t = Table.read(o, hdu = "EVENTS")
        if not io: 
            energies = t["ENERGY"].data
            times = t["TIME"].data
            conv = t["CONVERSION_TYPE"].data
            prob = t[srcname].data
        else:
            energies = np.concatenate([energies,t["ENERGY"].data])
            times = np.concatenate([times ,t["TIME"].data])
            conv = np.concatenate([conv ,t["CONVERSION_TYPE"].data])
            prob = np.concatenate([prob,t[srcname].data])
    return energies, times, conv, prob

def compute_profile2d(gta, target, prefix = '', sigma = 5., xsteps = 20, ysteps = 21):
    # save state
    saved_state = LikelihoodState(gta.like)

    s = gta.roi.get_source_by_name(target)
    if not s['SpectrumType']  == 'PowerLaw':
        logging.info("Changing spec to PL and re-fitting")
        gta = set_src_spec_pl(gta, gta.get_source_name(target))
        gta.optimize()
        gta.fit()

    loglike = -gta.like()

    # get the right index for the norm and index parameters of the target source
    name = gta.roi.get_source_by_name(target).name
    idx_norm = gta.like.par_index(name, 'Prefactor')
    idx_index = gta.like.par_index(name, 'Index')

    index = np.linspace(s.spectral_pars['Index']['value'] - sigma * s.spectral_pars['Index']['error'],
                    s.spectral_pars['Index']['value'] + sigma * s.spectral_pars['Index']['error'],
                    xsteps)

    #norm = np.logspace(np.log10(s.spectral_pars['Prefactor']['error'] * sigma/ \
                            #s.spectral_pars['Prefactor']['value']),
                        #np.log10(s.spectral_pars['Prefactor']['value'] / sigma \
                            #/ s.spectral_pars['Prefactor']['error']),
                        #ysteps)

    norm = np.logspace(np.log10(s.spectral_pars['Prefactor']['value'] / sigma / 2.), \
                        np.log10(s.spectral_pars['Prefactor']['value'] * sigma * 2.), \
                        ysteps)

    # limit range to allowed values by min max
    m = (index >= s.spectral_pars['Index']['min']) & \
        (index <= s.spectral_pars['Index']['max']) 
    if not m.sum():
        raise ValueError("No index scan points!")
    index = index[m]

    m = (norm >= s.spectral_pars['Prefactor']['min']) & \
        (norm <= s.spectral_pars['Prefactor']['max']) 
    if not m.sum():
        raise ValueError("No Prefactor scan points!")
    norm= norm[m]

    # do the profiling
    output = {'loglike' : []}
    logging.info("Starting profiling over {0:n} Index and {1:n}" \
                    " Prefactor values".format(index.size, norm.size))
    for i,x in enumerate(index):
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
            gta.like[idx_index] = x
            gta.like.freeze(idx_index)
        except RuntimeError:
            logging.warning(
                "Caught failure on set for %s::%s: %.2f" % (name, parName, x))
        o = gta.profile_norm(target, xvals=norm, reoptimize = True)
        output['loglike'].append(o['loglike'])
        gta.like.thaw(idx_index)

    output['loglike'] = np.array(output['loglike'])
    output['dloglike'] = output['loglike'] - loglike
    output['Index'] = index
    output['Prefactor'] = norm
    logging.info("Finished")
    saved_state.restore()
    np.save(path.join(gta.workdir, prefix + '_{0:s}_profile2d'.format(target.replace(' ','').lower())),
            output)
    return output

def get_best_fit_covar(gta, target, prefix = ''):
    """
    Return best fit values and covariance matrix for best-fit parameters 
    for a given source as a dictionary
    """
    s = gta.roi.get_source_by_name(target)
    name = gta.roi.get_source_by_name(target).name
    idx = []
    o = {'value':{}, 'error': {}}
    covar = np.array(gta.like.covariance)
    for k in list(s.spectral_pars.keys()):

        if s.spectral_pars[k]['free']:
            findpos = np.where(np.round(np.sqrt(np.diag(covar)), 4) == \
                    np.round(s.spectral_pars[k]['error'],4))

            if len(findpos):
                idx.append(findpos[0][0])

        for j in ['value', 'error']:
            o[j][k] = s.spectral_pars[k][j]
            if j == 'value':
                o[j][k] *= s.spectral_pars[k]['scale']
    if len(idx):
        o['covar'] = covar[np.min(idx):np.max(idx)+1, np.min(idx): np.max(idx)+1]
    np.save(path.join(gta.workdir, prefix + '_{0:s}_best-fit-pars'.format(target.replace(' ','').lower())), 
            o)
    return o
