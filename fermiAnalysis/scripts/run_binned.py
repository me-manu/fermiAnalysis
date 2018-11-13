import matplotlib
matplotlib.use('Agg')
from fermipy.gtanalysis import GTAnalysis
#from fermipy.skymap import Map
from gammapy.maps import WcsNDMap
from os import path
import logging
from haloanalysis.batchfarm import utils,lsf
import fermiAnalysis as fa
import argparse
import yaml
import os
import copy
import numpy as np
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
from fermiAnalysis import setup
from fermiAnalysis.utils import set_free_pars_avg,fit_with_retries
from fermiAnalysis.utils import set_src_spec_pl, set_src_spec_plexpcut
from fermiAnalysis.prepare import PreparePointSource 

def main():
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('--overwrite', required=False, default = 0, 
                        help='Overwrite existing files', type=int)
    parser.add_argument('--state', default = ['setup'],
                        choices = ['setup','avgspec','avgspec_ebl','lcbin', 'lcmonthly'],
                            help='Analysis state')
    parser.add_argument('-i', required=False, default = 0, 
                        help='Set local or scratch calculation', type=int)
    parser.add_argument('--specin', required=False, default = -2.1, 
                        help='Spectral index used for gtexposure in lieu of target in xml file',
                        type=float)
    parser.add_argument('--addnewsrcs', default = 0,  
                        help='Search for and add new sources and create residual map',
                        type=int)
    parser.add_argument('--reloadfit', default = 0,  
                        help='Reload ROI from avgspec xml file',
                        type=int)
    parser.add_argument('--relocalize', default = 0,  
                        help='Relocalize central source',
                        type=int)
    parser.add_argument('--createsed', default = 0,  
                        help='Create SED from best fit model',
                        type=int)
    parser.add_argument('--forcespec', default = 0,  
                        help='Recompute model parameters',
                        type=int)
    parser.add_argument('--freezesupexp', default = 0,  
                        help='freeze super exponential index parameters',
                        type=int)
    parser.add_argument('--restorecatspec', default = 0,  
                        help='Restore intitial catalog spectrum',
                        type=int)
    parser.add_argument('--sethardexpcutoff', default = 0,  
                        help='Manually change parameters of PL with SuperExpCutoff',
                        type=int)
    parser.add_argument('--pivotE_free', default = 0,  
                        help='let the pivot energy free during fit if spectrum is changed',
                        type=int)
    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i = args.i, logging_level = "INFO")
    logging.info('Running fermipy setup')

    if args.reloadfit:
        files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 or fn.find('.npy') > 0]
        if len(files):
            utils.copy2scratch(files, gta.workdir)
        else:
            logging.error("No files found in {0:s}".format(fit_config['avgspec']))
            args.reloadfit = False

    if args.state == 'lcbin':
        pps = PreparePointSource(config,logging={'verbosity' : 3})
        pps.setup()
        pps._bin_data_lc(overwrite = args.overwrite)
        pps._compute_lc_exp(overwrite = args.overwrite,
                        specin = args.specin)
        return pps

    elif args.state == 'setup':
        gta.setup()

    elif args.state.find('avgspec') >= 0:
        gta.setup()
        if not type(config['selection']['target']) == str:
            # target name not given
            # take closest source to ROI center if separation 
            # is less then 0.1 degree
            logging.warning("Target name is {0}".format(config['selection']['target']))
            sep = gta.roi.sources[0]['offset'] 
            logging.warning("Source closets to ROI center is {0:.3f} degree away".format(sep))
            if sep < 0.1:
                config['selection']['target'] = gta.roi.sources[0]['name']
                logging.info("Set target to {0:s}".format(config['selection']['target']))
            else: # add source at center of ROI
                csrc = SkyCoord(ra = config['selection']['ra'],
                        dec = config['selection']['dec'], frame = 'fk5', unit = 'degree')
                if csrc.dec.value < 0.:
                    sign = '-'
                else:
                    sign = '+'
                newname = 'j{0:02.0f}{1:02.0f}{2:s}{3:02.0f}{4:02.0f}'.format(
                                        csrc.ra.hms.h, csrc.ra.hms.m, sign, 
                                        np.abs(csrc.dec.dms.d), np.abs(csrc.dec.dms.m))
                gta.add_source(newname,{
                                'ra' : config['selection']['ra'], 'dec' : config['selection']['dec'],
                                'SpectrumType' : 'PowerLaw', 'Index' : fit_config['new_src_pl_index'],
                                'Scale' : fit_config['pivotE'] if 'pivotE' in fit_config.keys() else 1000.,
                                'Prefactor' : 1e-11,
                                'SpatialModel' : 'PointSource' })
                config['selection']['target'] = newname
                logging.info("Set target to {0:s}".format(config['selection']['target']))

        if args.reloadfit:
            # save old spectrum
            spec_cat = gta.roi.get_source_by_name(config['selection']['target'])
            try:
                gta.load_roi(args.state) # reload the average spectral fit
            except:
                logging.error("Could not reload fit. Continuing anyway.")
        logging.info('Running fermipy optimize and fit')
        # gives "failed to create spline" in get_parameter_limits function

        if 'source_spec' in fit_config.keys():
            m = gta.roi.get_source_by_name(config['selection']['target'])
            if not m['SpectrumType'] == fit_config['source_spec'] or args.forcespec:
                if fit_config['source_spec'] == 'PowerLaw':
                    gta = set_src_spec_pl(gta, gta.get_source_name(config['selection']['target']), 
                                        fit_config['pivotE'] if 'pivotE' in fit_config.keys() else None,
                                        e0_free = args.pivotE_free)
                elif fit_config['source_spec'] == 'PLSuperExpCutoff':
                    gta = set_src_spec_plexpcut(gta, gta.get_source_name(config['selection']['target']),
                                        fit_config['pivotE'] if 'pivotE' in fit_config.keys() else None,
                                        e0_free = args.pivotE_free)
                #elif fit_config['source_spec'] == 'LogParabola':
                    #gta = set_src_spec_lp(gta, gta.get_source_name(config['selection']['target']))
                else:
                    logging.warning("Spectrum {0:s} not supported, spectrum not changed".format(fit_config['source_spec']))

        # restore spectrum from catalog
        if args.restorecatspec:
            #for k in ['alpha','Index','Index1']:
                #if k in spec_cat.spectral_pars.keys():
                    #spec_cat.spectral_pars[k]['value'] -= 0.5
            #spec_cat.spectral_pars['Eb']['free'] = True
            #spec_cat.spectral_pars['Eb']['value'] = 2000.
            #spec_cat.spectral_pars['alpha']['value'] = 1.8
            #print spec_cat.spectral_pars
            gta.set_source_spectrum(config['selection']['target'],
                spectrum_type = spec_cat['SpectrumType'],
                spectrum_pars= spec_cat.spectral_pars)
            logging.info("restored catalog spectrum")

        # for some sources modeled with PL with super exponential cutoff I 
        # have to do this to get a nice SED, but not for 3C454.3!
        if gta.roi.get_source_by_name(config['selection']['target'])['SpectrumType'] ==\
            'PLSuperExpCutoff' and args.sethardexpcutoff :
            pars = {}
            old_spec_pars = copy.deepcopy(gta.roi.get_source_by_name(config['selection']['target']))
            for k in ['Prefactor','Scale','Index1','Index2','Cutoff']:
                pars[k] = old_spec_pars.spectral_pars[k]
            if config['selection']['target'] == '3C454.3':
                # values from Romoli et al 2017
                pars['Prefactor']['value'] = 4.7
                pars['Index1']['value'] = 1.87
                pars['Index2']['value'] = 0.4
                pars['Cutoff']['value'] = 1100.
                pars['Cutoff']['scale'] = 1.
                pars['Cutoff']['min'] = 100.
                pars['Cutoff']['max'] = 10000.
            else:
                pars['Index1']['value'] = 1.8
                pars['Index2']['value'] = 1.
                pars['Cutoff']['value'] = 5e4
            gta.set_source_spectrum(config['selection']['target'],
            #    spectrum_type = 'PLSuperExpCutoff2',
                spectrum_type = 'PLSuperExpCutoff',
                spectrum_pars= pars)
            logging.info("changed spectral parameters to {0}".format(
                gta.roi.get_source_by_name(config['selection']['target']).spectral_pars))
        else:
            old_spec_pars = None

        if args.state.find('ebl') >= 0:
            gta = fa.utils.add_ebl_atten(gta,config['selection']['target'],fit_config['z'])

        gta = set_free_pars_avg(gta, fit_config, freezesupexp = args.freezesupexp)
        f,gta = fit_with_retries(gta, fit_config, config['selection']['target'], 
                                    alt_spec_pars = old_spec_pars)

        logging.debug(f)

        #relocalize central source and refit
        if args.relocalize and type(config['selection']['target']) == str:
            loc = gta.localize(config['selection']['target'], make_plots=True,
                                free_background = fit_config['reloc_bkg'],
                                free_radius = fit_config['reloc_rad'],
                                update=True)

            logging.info('new position is {0[pos_offset]:.3f} degrees from old position'.format(loc))
            logging.info('Pos uncertainty is {0[pos_r68]:.3f} (68%); {0[pos_r95]:.3f} (95%) degrees'.format(loc))
            logging.info('Refitting with new source position ...')
            logging.info('free source parameters:')
            for s in  gta.get_sources() :                                                           
                for k in s.spectral_pars.keys():
                    if s.spectral_pars[k]['free']:   
                        logging.info('{0:s}: {1:s}'.format(s.name, k))
            f = gta.fit()
            #f,gta = fit_with_retries(gta, fit_config, config['selection']['target'])

        gta.print_roi()
        gta.write_roi(args.state)

    elif args.state == 'lcmonthly':
        gta.setup()
        gta.load_roi('avgspec') # reload the average spectral fit
        logging.info('Running the 30-day bin' + \
            'light curve for {0[target]:s}'.format(config['selection']))
        lc = gta.lightcurve(config['selection']['target'],
                                binsz = 30. * 24. * 60. * 60.)

    model = {'Scale': 1000., 'Index' : fit_config['new_src_pl_index'], 'SpatialModel' : 'PointSource'}
    if args.addnewsrcs:
        gta.load_roi(args.state) # reload the average spectral fit

        max_sqrt_ts = 1000.
        irun = 0
        # define the test source
        #model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}

        # run ts map and add new sources with sqrt(ts) > 5
        # reoptimize iteratively for each new source
        # this is only done for outer RoI
        while max_sqrt_ts >= fit_config['max_sqrt_ts']:
            # run ts and residual maps
            ts_maps = gta.tsmap(args.state,model=model, 
                write_fits = True, write_npy = True, make_plots = True)
            # get the skydirs
            #coords = ts_maps['sqrt_ts'].get_pixel_skydirs()
            coords = ts_maps['sqrt_ts'].geom.get_coord()
            if ts_maps['sqrt_ts'].geom.coordsys == 'CEL':
                frame = 'fk5'
            elif ts_maps['sqrt_ts'].geom.coordsys == 'GAL':
                frame = 'galactic'
            c = SkyCoord(coords[0], coords[1], unit = 'deg', 
                frame = frame)

            #sqrt_ts = ts_maps['sqrt_ts'].get_map_values(coords.ra, coords.dec) # these are all nans. workaround: load fits file
            #sqrt_ts_map = Map.create_from_fits(path.join(gta.workdir,ts_maps['file']),hdu = 'SQRT_TS_MAP')
            #n_map = Map.create_from_fits(path.join(gta.workdir,ts_maps['file']),hdu = 'N_MAP')

            sqrt_ts_map = WcsNDMap.read(path.join(gta.workdir,ts_maps['file']),hdu = 'SQRT_TS_MAP')
            n_map = WcsNDMap.read(path.join(gta.workdir,ts_maps['file']),hdu = 'N_MAP')

            sqrt_ts = sqrt_ts_map.data
            amplitudes = n_map.data
            
            # get the angular separation from RoI center
            sep = gta.roi.skydir.separation(c)
            # mask nan values and pixels close to central source
            m = np.isfinite(sqrt_ts) & (sep.value > fit_config['new_src_search_rad'])
            if not np.sum(m):
                logging.warning('No pixels that are finite at distance > {0[new_src_search_rad]:.2f}'.format(fit_config))
                raise RuntimeError

            # get max ts value 
            max_sqrt_ts = np.max(sqrt_ts[m])

            if max_sqrt_ts < fit_config['max_sqrt_ts']: break

            # get the coords of max ts
            idx = np.argmax(sqrt_ts[m])
            logging.info('Found new source with sqrt(ts) = {0:.2f} at ra,dec = {1:.2f}, {2:.2f}'.format(
                            sqrt_ts[m][idx], c.ra[m][idx].value, c.dec[m][idx].value))
            # add a new source
            csrc = SkyCoord(ra = c.ra[m][idx], dec = c.dec[m][idx], frame = 'fk5')
            if csrc.dec.value < 0.:
                sign = '-'
            else:
                sign = '+'
            newname = 'j{0:02.0f}{1:02.0f}{2:s}{3:02.0f}{4:02.0f}'.format(
                                        csrc.ra.hms.h, csrc.ra.hms.m, sign, 
                                        np.abs(csrc.dec.dms.d), np.abs(csrc.dec.dms.m))
            gta.add_source(newname,{
                                'ra' : c.ra[m][idx].value, 'dec' : c.dec[m][idx].value,
                                'SpectrumType' : 'PowerLaw', 'Index' : fit_config['new_src_pl_index'],
                                'Scale' : 1000, 'Prefactor' : amplitudes[m][idx],
                                'SpatialModel' : 'PointSource' })
            logging.debug('Amplitude of source: {0}'.format(amplitudes[m][idx]))

            gta.free_source(newname, pars=['norm','index'], free = True)
            f = gta.fit()
            gta.print_roi()
            irun += 1

        # if new sources where added, save output
        if irun > 0:
            gta.print_roi()
#            gta = reset_diff_filenames(gta)
            # refit the model with new sources present
            gta = set_free_pars_avg(gta, fit_config)
            f,gta = fit_with_retries(gta, fit_config, config['selection']['target'])
            gta.write_roi(args.state)

    else:
        ts_maps = gta.tsmap(args.state,model=model, 
            write_fits = True, write_npy = True, make_plots = True)
    try:
        resid_maps = gta.residmap(args.state,model=model, make_plots=True, write_fits = True, write_npy = True)
    except:
        logging.error("Residual map computation and plotting failed")


    if args.createsed:
        gta.load_roi(args.state) # reload the average spectral fit
        logging.info('Running sed for {0[target]:s}'.format(config['selection']))
        sed = gta.sed(config['selection']['target'],
                        prefix = args.state,
                        #outfile = 'sed.fits',
                        #free_radius =  #sed_config['free_radius'],
                        #free_background= #sed_config['free_background'],
                        free_pars = fa.allnorm,
                        #make_plots = sed_config['make_plots'],
                        #cov_scale = sed_config['cov_scale'],
                        #use_local_index = sed_config['use_local_index'],
                        #use_local_index = True, # sed_config['use_local_index'],
                        #bin_index = sed_config['bin_index']
                        )

    # run the analysis for the full flare durations
    #if args.state == 'fullflare-avg':
    #if args.state.find('-avg') >= 0:
        ## stage the full time array analysis results to the tmp dir
        ## do not copy png images
        #files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 or fn.find('.npy') > 0]
        #utils.copy2scratch(files, gta.workdir)
#
        #if args.state == 'fullflare-avg':
            #gta.load_roi('avgspec') # reload the average spectral fit
        #elif args.state == 'gtiflare-avg':
            #gta.load_roi('fullflare-avg') # reload the average spectral fit
#
        #o = gta.optimize() # perform an initial fit
        #logging.debug(o)
        #gta.print_roi()
#
        ## Free all parameters of all Sources within X deg of ROI center
        ##gta.free_sources(distance=fit_config['ps_dist_all'])
        ## Free Normalization of all Sources within X deg of ROI center
        #gta.free_sources(distance=fit_config['ps_dist_norm_fullflare'],pars=fa.allnorm)
        ## Free spectra parameters of all Sources within X deg of ROI center
        #gta.free_sources(distance=fit_config['ps_dist_idx_fullflare'],pars=fa.allidx)
        ## Free all parameters of isotropic and galactic diffuse components
        #gta.free_source('galdiff', pars= fa.allnorm, free = fit_config['gal_norm_free_fullflare'])
        #gta.free_source('galdiff', pars= fa.allidx, free = fit_config['gal_idx_free_fullflare'])
        #gta.free_source('isodiff', pars= fa.allnorm, free = fit_config['iso_norm_free_fullflare'])
        ## Free sources with TS > X
        #gta.free_sources(minmax_ts=[fit_config['ts_norm'],None], pars = fa.allnorm)
        ## Fix sources with TS < Y
        #gta.free_sources(minmax_ts=[None,fit_config['ts_fixed']],free=False,
                            #pars=fa.allidx + fa.allnorm)
        ## Fix indeces for sources with TS < Z
        #gta.free_sources(minmax_ts=[None,fit_config['ts_fixed_idx']],free=False,
                #pars = fa.allidx)
        ## Fix sources Npred < Z
        #gta.free_sources(minmax_npred=[None,fit_config['npred_fixed']],free=False,
                            #pars=fa.allidx + fa.allnorm)
#
        ## gives "failed to create spline" in get_parameter_limits function
        ##gta = fa.utils.add_ebl_atten(gta,config['selection']['target'],fit_config['z'])
        #logging.info('free source parameters:')
        #for s in  gta.get_sources() :                                                           
            #for k in s.spectral_pars.keys():
                #if s.spectral_pars[k]['free']:   
                    #logging.info('{0:s}: {1:s}'.format(s.name, k))
        #f = gta.fit()
#
        #retries = f['config']['retries']
        #tsfix = fit_config['ts_fixed']
        #while not f['fit_success'] and retries > 0:
            #gta.free_source(config['selection']['target'], pars = ['beta', 'Index2'], free = False)
            ## Fix more sources
            #tsfix *= 3
            #gta.free_sources(minmax_ts=[None,tsfix],free=False,
                            #pars=fa.allidx + fa.allnorm)
            #logging.info("retrying fit")
            #for s in  gta.get_sources() :                                                           
                #for k in s.spectral_pars.keys():
                    #if s.spectral_pars[k]['free']:   
                        #logging.info('{0:s}: {1:s}'.format(s.name, k))
            #o = gta.optimize()
            #gta.print_roi()
            #f = gta.fit()
            #retries -= 1
#
        #gta.write_roi(args.state)
#
        #logging.info('Running sed for {0[target]:s}'.format(config['selection']))
        #sed = gta.sed(config['selection']['target'], prefix = args.state)
        #model = {'Scale': 1000., 'Index' : fit_config['new_src_pl_index'], 'SpatialModel' : 'PointSource'}
        #resid_maps = gta.residmap(args.state,model=model, make_plots=True, write_fits = True, write_npy = True)
    return gta

if __name__ == '__main__':
    gta = main()
