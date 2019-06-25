import matplotlib
matplotlib.use('Agg')
import fermipy
from fermipy.gtanalysis import GTAnalysis
from fermipy.skymap import Map
from astropy.table import Table
from os import path, environ
import logging
from haloanalysis.batchfarm import utils,lsf
import fermiAnalysis as fa
from fermiAnalysis.utils import * 
import argparse
import yaml
import os
import copy
import numpy as np
from glob import glob
from fermiAnalysis.prepare import PreparePointSource 

def rebin(c,b,minc = 10):
    """
    Rebin the light curve with at least minc counts per bin
    
    Parameters
    c: list with counts
    b: list with bin bounds
    minc: min counts per bin
    """
    csum = 0
    newt = [b[0]]
    argt = [0]
    for i,ci in enumerate(c[:-1]):
        csum += ci
        cremain = c[i+1:].sum()

        if csum >= minc:
            if cremain >= minc:
                newt.append(b[i+1])
                argt.append(i+1)
                csum = 0

    newt.append(b[-1])
    argt.append(b.size)

    return newt

def rm_trailing_zeros(counts):
    """Remove bins with trailing zeros"""
    mcounts_pre = np.cumsum(counts) > 4.
    mcounts_post = np.cumsum(counts[::-1])[::-1] > 4.
    if not np.sum(mcounts_post & mcounts_pre):
        logging.error("*** After removing trailing bins with less than 4 counts, no bins left! Exiting")
        raise Exception
    return mcounts_post, mcounts_pre

def calc_counts(t,bins):
    """
    Use photons in table t and bins to calculate number of counts 
    per bin
    """
    counts = []
    for i,ti in enumerate(bins[:-1]):
        m = (t["TIME"] >= ti) & (t["TIME"] < bins[i+1])
        counts.append(np.sum(m))
    return np.array(counts)

def main():
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the lc analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('-i', required=False, default = 0, 
                        help='Set local or scratch calculation', type=int)
    parser.add_argument('--forcepl', default = 0,  
                        help='Force the target source to have power-law shape',
                        type=int)
    parser.add_argument('--createsed', default = 0,  
                        help='Create SED from best fit model',
                        type=int)
    parser.add_argument('--adaptive', default = 0,  
                        help='Use adaptive binning for minute scale light curves',
                        type=int)
    parser.add_argument('--srcprob', default = 0,  
                        help='Calculate the source probability for the photons,' \
                            ' only works when no sub orbit time scales are used',
                        type=int)
    parser.add_argument('--mincounts', default = 2,  
                        help='Minimum number of counts within LC bin to run analysis',
                        type=int)
    parser.add_argument('--simulate', default = None,  
                        help='None or full path to yaml file which contains src name' \
                        'and spec to be simulated',
                        )
    parser.add_argument('--randomize', default = 1,  
                        help='If you simulate, use Poisson realization. If false, use Asimov data set',
                        type=int)
    args = parser.parse_args()

    utils.init_logging('DEBUG')
    config = yaml.load(open(args.conf))
    tmpdir, job_id = lsf.init_lsf()
    if not job_id:
        job_id = args.i
    logging.info('tmpdir: {0:s}, job_id: {1:n}'.format(tmpdir,job_id))
    os.chdir(tmpdir)    # go to tmp directory
    logging.info('Entering directory {0:s}'.format(tmpdir))
    logging.info('PWD is {0:s}'.format(os.environ["PWD"]))

    # copy the ft1,ft2 and ltcube files 
    #for k in ['evfile','scfile','ltcube']:
# don't stage them, done automatically by fermipy if needed
#        config[k] = utils.copy2scratch(config[k], tmpdir)
    # set the scratch directories
    logging.debug(config['data'])
    config['fileio']['scratchdir'] = tmpdir

    # set the log file
    logdir = copy.deepcopy(config['fileio']['logfile'])
    config['fileio']['logfile'] = path.join(tmpdir,'fermipy.log')
    # debugging: all files will be saved (default is False)
    #config['fileio']['savefits'] = True

    # if simulating an orbit, save fits files
    if args.simulate is not None:
        config['fileio']['savefits'] = True

    # copy all fits files already present in outdir
    # run the analysis
    lc_config = copy.deepcopy(config['lightcurve'])
    fit_config = copy.deepcopy(config['fit_pars'])

    # remove parameters from config file not accepted by fermipy
    for k in ['configname', 'tmp', 'log', 'fit_pars']: 
        config.pop(k,None)
    if 'adaptive' in config['lightcurve'].keys():
        config['lightcurve'].pop('adaptive', None)

    # set the correct time bin
    config['selection']['tmin'],config['selection']['tmax'],nj = set_lc_bin(
                                                    config['selection']['tmin'],
                                                    config['selection']['tmax'],
                                                    config['lightcurve']['binsz'], 
                                                    job_id - 1 if job_id > 0 else 0, 
                                                    ft1 = config['data']['evfile'])
    logging.debug('setting light curve bin' + \
        '{0:n}, between {1[tmin]:.0f} and {1[tmax]:.0f}'.format(job_id, config['selection']))
    if args.adaptive:
        config['fileio']['outdir'] = utils.mkdir(path.join(config['fileio']['outdir'],'adaptive{0:.0f}/'.format(
                                                    lc_config['adaptive'])))
    config['fileio']['outdir'] = utils.mkdir(path.join(config['fileio']['outdir'],
                        '{0:05n}/'.format(job_id if job_id > 0 else 1)))

    logging.info('Starting with fermipy analysis')
    logging.info('using fermipy version {0:s}'.format(fermipy.__version__))
    logging.info('located at {0:s}'.format(fermipy.__file__))

    if config['data']['ltcube'] == '':
        config['data'].pop('ltcube',None)

    compute_sub_gti_lc = False
    if type(config['lightcurve']['binsz']) == str:
        if len(config['lightcurve']['binsz'].strip('gti')):
            compute_sub_gti_lc = True
            if config['lightcurve']['binsz'].find('min') > 0:
                config['lightcurve']['binsz'] = float(
                    config['lightcurve']['binsz'].strip('gti').strip('min')) * 60.
                logging.info("set time bin length to {0:.2f}s".format(
                                        config['lightcurve']['binsz']))
        else:
            config['lightcurve']['binsz'] = 3. * 3600.
    try:
        gta = GTAnalysis(config,logging={'verbosity' : 3})
    except Exception as e:
        logging.error("{0}".format(e))
        config['selection']['target'] = None
        gta = GTAnalysis(config,logging={'verbosity' : 3})
        sep = gta.roi.sources[0]['offset'] 
        logging.warning("Source closets to ROI center is {0:.3f} degree away".format(sep))
        if sep < 0.1:
            config['selection']['target'] = gta.roi.sources[0]['name']
            gta.config['selection']['target'] = config['selection']['target']
            logging.info("Set target to {0:s}".format(config['selection']['target']))

    # stage the full time array analysis results to the tmp dir
    # do not copy png images
    files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 or fn.find('.npy') > 0]
    files += [config['data']['evfile']]
    utils.copy2scratch(files, gta.workdir)

    # check before the analysis start if there are any events in the master file
    # in the specified time range
    logging.info('Checking for events in initial ft1 file')
    t = Table.read(path.join(gta.workdir,
                    path.basename(config['data']['evfile'])),
                    hdu = 'EVENTS')
    logging.info("times in base ft1: {0} {1} {2}".format(t["TIME"].max(), t["TIME"].min(), t["TIME"].max() - t["TIME"].min()))
    m = (t["TIME"] >= config['selection']['tmin']) & (t["TIME"] <= config['selection']['tmax'])
    if np.sum(m) < args.mincounts + 1:
        logging.error("*** Only {0:n} events between tmin and tmax! Exiting".format(np.sum(m)))
        assert np.sum(m) > args.mincounts
    else:
        logging.info("{0:n} events between tmin and tmax".format(np.sum(m)))

    # check how many bins are in each potential light curve bin
    if compute_sub_gti_lc:
# select time of first and last 
# photon instead of GTI time
        m = (t["TIME"] >= config['selection']['tmin']) & \
             (t["TIME"] <= config['selection']['tmax'])

        tmin = t["TIME"][m].min() - 1.
        tmax = t["TIME"][m].max() + 1.
        logging.info("There will be up to {0:n} time bins".format(np.ceil(
            (tmax - tmin) / \
            config['lightcurve']['binsz'])))

        bins = np.arange(tmin,tmax,config['lightcurve']['binsz'])
        bins = np.concatenate([bins, [config['selection']['tmax']]])
        counts = calc_counts(t,bins)
        # remove the starting times of the bins with zero counts
        # and rebin the data
        logging.info("Counts before rebinning: {0}".format(counts))
        mincounts = 10.
        mc = counts < mincounts
        if np.sum(mc):
            # remove trailing zeros
            if np.any(counts == 0.):
                mcounts_post, mcounts_pre = rm_trailing_zeros(counts)
                counts = counts[mcounts_post & mcounts_pre]
                bins = np.concatenate([bins[:-1][mcounts_post & mcounts_pre],
                    [bins[1:][mcounts_post & mcounts_pre].max()]])
            bins = rebin(counts, bins)
            logging.info("Bin lengths after rebinning: {0}".format(np.diff(bins)))
            logging.info("Bin times after rebinning: {0}".format(bins))
            counts = calc_counts(t, bins)
            logging.info("Counts after rebinning: {0}".format(counts))
        else:
            logging.info("Regular time binning will be used")
        bins = list(bins)

    logging.info('Running fermipy setup')
    try:
        gta.setup()
    except (RuntimeError,IndexError) as e:
        logging.error('Caught Runtime/Index Error while initializing analysis object')
        logging.error('Printing error:')
        logging.error(e)
        if e.message.find("File not found") >= 0 and e.message.find('srcmap') >= 0:
            logging.error("*** Srcmap calculation failed ***")

        logging.info("Checking if there are events in ft1 file")
        ft1 = path.join(gta.workdir,'ft1_00.fits')
        f = glob(ft1)
        if not len(f):
            logging.error("*** no ft1 file found at location {0:s}".format(ft1))
            raise
        t = Table.read(f[0], hdu = 'EVENTS')
        if not len(t):
            logging.error("*** The ft1 file contains no events!! ***".format(len(t)))
        else:
            logging.info("The ft1 file contains {0:n} event(s)".format(len(t)))
        raise


    logging.info('Loading the fit for the average spectrum')
    gta.load_roi('avgspec') # reload the average spectral fit
    logging.info('Running fermipy optimize and fit')

    if args.forcepl:
        gta = set_src_spec_pl(gta, gta.get_source_name(config['selection']['target']))
# to do add EBL absorption at some stage ...
#        gta = add_ebl_atten(gta, gta.get_source_name(config['selection']['target']), fit_config['z'])

    if compute_sub_gti_lc:
        if args.adaptive:
            # do import only here since root must be compiled
            from fermiAnalysis import adaptivebinning as ab
             # compute the exposure
            energy = 1000.
            texp, front, back = ab.comp_exposure_phi(gta, energy = 1000.)
            # compute the bins
            result = ab.time_bins(gta, texp,
                                    0.5 * (front + back), 
                                    #critval = 20., # bins with ~20% unc
                                    critval = lc_config['adaptive'],
                                    Epivot = None, # compute on the fly
            #                        tstart = config['selection']['tmin'],
            #                        tstop = config['selection']['tmax']
                                    )

            # cut the bins to this GTI
            mask = result['tstop'] > config['selection']['tmin']
            mask = mask & (result['tstart'] < config['selection']['tmax'])

            # try again with catalog values
            if not np.sum(mask):
                logging.error("Adaptive bins outside time window, trying catalog values for flux")
                result = ab.time_bins(gta, texp,
                                    0.5 * (front + back), 
                                    critval = lc_config['adaptive'], # bins with ~20% unc
                                    Epivot = None, # compute on the fly
                                    forcecatalog = True,
            #                        tstart = config['selection']['tmin'],
            #                        tstop = config['selection']['tmax']
                                    )

                # cut the bins to this GTI
                mask = result['tstop'] > config['selection']['tmin']
                mask = mask & (result['tstart'] < config['selection']['tmax'])
                if not np.sum(mask):
                    logging.error("Adaptive bins do not cover selected time interval!")
                    logging.error("Using original bins")

                else:
                    bins = np.concatenate((result['tstart'][mask], [result['tstop'][mask][-1]]))
                    bins[0] = np.max([config['selection']['tmin'], bins[0]])
                    bins[-1] = np.min([config['selection']['tmax'], bins[-1]])
                    bins = list(bins)

                    # removing trailing zeros
                    counts = calc_counts(t,bins)
                    mcounts_post, mcounts_pre = rm_trailing_zeros(counts)
                    logging.info("count masks: {0} {1}, bins: {2}, counts: {3}".format(mcounts_post, mcounts_pre, bins, counts))
                    counts = counts[mcounts_post & mcounts_pre]
                    bins = np.concatenate([np.array(bins)[:-1][mcounts_post & mcounts_pre],
                        [np.array(bins)[1:][mcounts_post & mcounts_pre].max()]])
                    logging.info("Using bins {0}, total n={1:n} bins".format(bins, len(bins)-1))
                    logging.info("bins widths : {0}".format(np.diff(bins)))
                    logging.info("counts per bin: {0} ".format(calc_counts(t,bins)))
                    bins = list(bins)

# TODO: test that this is working also with GTIs that have little or no counts

        lc = gta.lightcurve(config['selection']['target'],
                                binsz = config['lightcurve']['binsz'],
                                free_background = config['lightcurve']['free_background'],
                                free_params = config['lightcurve']['free_params'],
                                free_radius = config['lightcurve']['free_radius'],
                                make_plots = True,
                                multithread = True,
                                nthread = 4,
                                #multithread = False,
                                #nthread = 1,
                                save_bin_data = True,
                                shape_ts_threshold = 16.,
                                use_scaled_srcmap = True,
                                use_local_ltcube = True,
                                write_fits = True, 
                                write_npy = True,
                                time_bins = bins, 
                                outdir = '{0:.0f}s'.format(config['lightcurve']['binsz'])
                                )
    else:
    # run the fitting of the entire time and energy range
        try:
            o = gta.optimize() # perform an initial fit
            logging.debug(o)
        except RuntimeError as e:
            logging.error("Error in optimize: {0}".format(e))
            logging.info("Trying to continue ...")

        # Freeze all parameters
        gta.free_sources(free = False, pars = fa.allidx + fa.allnorm)
        # Free parameters of central source
        gta.free_source(config['selection']['target'], pars = config['lightcurve']['free_params'])
        # Free Normalization of all Sources within X deg of ROI center
        gta.free_sources(distance=config['lightcurve']['free_radius'],pars=fa.allnorm)
        # Fix sources with TS < Y
        gta.free_sources(minmax_ts=[None,fit_config['ts_fixed']],free=False,
                pars = fa.allnorm + fa.allidx)
        # Fix sources Npred < Z
        gta.free_sources(minmax_npred=[None,fit_config['npred_fixed']],free=False,
                pars = fa.allnorm + fa.allidx)

        # Free all parameters of isotropic and galactic diffuse components
        if config['lightcurve']['free_background']:
            gta.free_source('galdiff', pars=fa.allnorm, free = True)
#        gta.free_source('galdiff', pars=['index'], free = True)
            gta.free_source('isodiff', pars=fa.allnorm, free = True)

        print_free_sources(gta)
        f = gta.fit()
        gta, f = refit(gta, config['selection']['target'],f, fit_config['ts_fixed'])
        gta.print_roi()
        gta.write_roi('lc')

        if args.createsed:
            gta.load_roi('lc') # reload the average spectral fit
            logging.info('Running sed for {0[target]:s}'.format(config['selection']))
            sed = gta.sed(config['selection']['target'],
                        prefix = 'lc_sed',
                        free_radius = None if config['sed']['free_radius'] == 0. \
                            else config['sed']['free_radius'],
                        free_background= config['sed']['free_background'],
                        free_pars = fa.allnorm,
                        make_plots = config['sed']['make_plots'],
                        cov_scale = config['sed']['cov_scale'],
                        use_local_index = config['sed']['use_local_index'],
                        bin_index = config['sed']['bin_index']
                        )

    # debugging: calculate sed and resid maps for each light curve bin
    #logging.info('Running sed for {0[target]:s}'.format(config['selection']))
    #sed = gta.sed(config['selection']['target'], prefix = 'lc')
    #model = {'Scale': 1000., 'Index' : fit_config['new_src_pl_index'], 'SpatialModel' : 'PointSource'}
    #resid_maps = gta.residmap('lc',model=model, make_plots=True, write_fits = True, write_npy = True)

        if args.srcprob:
            logging.info("Running srcprob with srcmdl {0:s}".format('lc'))
            gta.compute_srcprob(xmlfile = 'lc', overwrite = True)

        # simulate a source 
        if args.simulate is not None:
            with open(args.simulate) as f:
                simsource = np.load(f, allow_pickle = True).flat[0]

            # set the source to the simulation value
            gta.set_source_spectrum(simsource['target'],
                spectrum_type = simsource['spectrum_type'],
                spectrum_pars = simsource['spectrum_pars'][job_id - 1])

            logging.info("changed spectral parameters to {0}".format(
                gta.roi.get_source_by_name(simsource['target']).spectral_pars))

            # simulate the ROI
            gta.simulate_roi(randomize = bool(args.randomize))
            # freeze spectral shapes
            #gta.free_sources(free = False, pars = fa.allidx + fa.allnorm)
            # Free parameters of central source
            gta.free_source(config['selection']['target'],
                pars = config['lightcurve']['free_params'])
            # Free Normalization of all Sources within X deg of ROI center
            #gta.free_sources(distance=0.5,pars=fa.allnorm)
            # Fix spectral shape parameters of central source
            #gta.free_source(config['selection']['target'],
            #    pars = fa.allidx,
            #    free = False)
                
            # fit the simulation
            f = gta.fit()
            gta, f = refit(gta, config['selection']['target'],f, fit_config['ts_fixed'])
            gta.print_roi()
            gta.write_roi('lc_simulate_{0:s}'.format(simsource['suffix']))
    return gta 

if __name__ == '__main__':
    gta = main()
