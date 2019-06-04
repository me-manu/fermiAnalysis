import matplotlib
matplotlib.use('Agg')
from fermipy.gtanalysis import GTAnalysis
from fermipy.skymap import Map
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
from fermiAnalysis.prepare import PreparePointSource 

def main():
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('--overwrite', required=False, default = 0, 
                        help='Overwrite existing files', type=int)
    parser.add_argument('--state', default = ['setup'],
                           choices = ['setup','lcbinexp','avgspec','tsresid','avgsed','lcbin','srcprob'],
                               help='Analysis state')
    parser.add_argument('-i', required=False, default = 0, 
                        help='Set local or scratch calculation', type=int)
    parser.add_argument('--dt', required=False, default = 0., 
                        help='time interval for light curve binning, if none, use binning of config file', type=float)
    parser.add_argument('--specin', required=False, default = -2.1, 
                        help='Spectral index used for gtexposure in lieu of target in xml file',
                           type=float)
    parser.add_argument('--targetname', required=False, default = "", 
                        help='Source name used for binned light curve exposure calculation')
    parser.add_argument('--srcmdl', required=False, default = "srcmdl", 
                        help='srcmdl name for gtsrcprob / gtexposure calculation' )
    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i = args.i, logging_level = "INFO")
    logging.info('Running fermipy setup')
    pps = PreparePointSource(config,logging={'verbosity' : 3})

    if args.state == 'setup':
         pps.setup()
         return pps

    if args.state == 'lcbinexp':
         pps._bin_data_lc(overwrite = args.overwrite, dtime = args.dt)
         if 'target' in config['selection'].keys():
             target = pps.roi.sources[0].name # assume central source
             files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 ]
             utils.copy2scratch(files, gta.workdir)
         else:
             target = args.targetname

         logging.info("Running lc bin with target {0:s} and specin {1:.2f}".format(target,args.specin))
         pps._compute_lc_exp(overwrite = args.overwrite,
                           specin = args.specin, 
                           target = target,
                           srcmdl = args.srcmdl
                           )
         return pps

    if args.state == 'srcprob':
         # copy srcmodels of average fit
         files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 ]
         utils.copy2scratch(files, gta.workdir)

         logging.info("Running diffrsp with srcmdl {0:s}".format(args.srcmdl))
         pps._compute_diffrsp(args.srcmdl,
                           overwrite = args.overwrite)
         logging.info("Running srcprob with srcmdl {0:s}".format(args.srcmdl))
         pps._compute_srcprob(args.srcmdl,
                           overwrite = args.overwrite)
         return pps

    elif args.state == 'avgspec':
         gta.setup()
         logging.info('Running fermipy optimize and fit')
         # run the fitting of the entire time and energy range
         o = gta.optimize() # perform an initial fit
         logging.debug(o)
         gta.print_roi()

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

         # gives "failed to create spline" in get_parameter_limits function
         #gta = fa.utils.add_ebl_atten(gta,config['selection']['target'],fit_config['z'])
         f = gta.fit()
         #logging.debug(f)
         #relocalize central source and refit
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

         gta.print_roi()
         gta.write_roi(args.state)

    if args.state == 'tsresid':
         gta.setup()
         gta.load_roi('avgspec') # reload the average spectral fit

         max_sqrt_ts = 1000.
         irun = 0
         # define the test source
         model = {'Scale': 1000., 'Index' : fit_config['new_src_pl_index'], 'SpatialModel' : 'PointSource'}
         #model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}

         # run ts map and add new sources with sqrt(ts) > 5
         # reoptimize iteratively for each new source
         # this is only done for outer RoI
         while max_sqrt_ts >= fit_config['max_sqrt_ts']:
             # run ts and residual maps
             ts_maps = gta.tsmap('avgspec',model=model, 
                  write_fits = True, write_npy = True, make_plots = True)
             # get the skydirs
             coords = ts_maps['sqrt_ts'].get_pixel_skydirs()
             #sqrt_ts = ts_maps['sqrt_ts'].get_map_values(coords.ra, coords.dec) # these are all nans. workaround: load fits file
             sqrt_ts_map = Map.create_from_fits(path.join(gta.workdir,ts_maps['file']),hdu = 'SQRT_TS_MAP')
             n_map = Map.create_from_fits(path.join(gta.workdir,ts_maps['file']),hdu = 'N_MAP')

             sqrt_ts = sqrt_ts_map.get_map_values(coords.ra, coords.dec)
             amplitudes = n_map.get_map_values(coords.ra, coords.dec)
             
             # get the angular separation from RoI center
             sep = gta.roi.skydir.separation(coords)
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
                               sqrt_ts[m][idx], coords.ra[m][idx].value, coords.dec[m][idx].value))
             # add a new source
             c = SkyCoord(ra = coords.ra[m][idx], dec = coords.dec[m][idx], frame = 'icrs')
             if c.dec.value < 0.:
                  sign = '-'
             else:
                  sign = '+'
             newname = 'j{0:02.0f}{1:02.0f}{2:s}{3:02.0f}{4:02.0f}'.format(c.ra.hms.h, c.ra.hms.m, sign, 
                                                                                              np.abs(c.dec.dms.d), np.abs(c.dec.dms.m))
             gta.add_source(newname,{
                                    'ra' : coords.ra[m][idx].value, 'dec' : coords.dec[m][idx].value,
                                    'SpectrumType' : 'PowerLaw', 'Index' : fit_config['new_src_pl_index'],
                                    'Scale' : 1000, 'Prefactor' : amplitudes[m][idx],
                                    'SpatialModel' : 'PointSource' })
             logging.debug('Amplitude of source: {0}'.format(amplitudes[m][idx]))

             gta.free_source(newname, pars=['norm','index'], free = True)
             f = gta.fit()
             gta.print_roi()
             irun += 1

         resid_maps = gta.residmap('avgspec',model=model, make_plots=True, write_fits = True, write_npy = True)

         # if new sources where added, save output
         if irun > 0:
             gta.print_roi()
#             gta = reset_diff_filenames(gta)
             gta.write_roi('avgspec')

    if args.state == 'avgsed':
         gta.setup()
         gta.load_roi('avgspec') # reload the average spectral fit
         logging.info('Running sed for {0[target]:s}'.format(config['selection']))
         sed = gta.sed(config['selection']['target'],
                           #outfile = 'sed.fits',
                           #free_radius = sed_config['free_radius'],
                           #free_background= sed_config['free_background'],
                           #make_plots = sed_config['make_plots'],
                           #cov_scale = sed_config['cov_scale'],
                           #use_local_index = sed_config['use_local_index'],
                           #bin_index = sed_config['bin_index']
                           )
    if args.state == 'lcmonthly':
         gta.setup()
         gta.load_roi('avgspec') # reload the average spectral fit
         logging.info('Running the 30-day bin light curve for {0[target]:s}'.format(config['selection']))
         lc = gta.lightcurve(config['selection']['target'],
                                    binsz = 30. * 24. * 60. * 60.)

    # run the analysis for the full flare durations
    #if args.state == 'fullflare-avg':
    if args.state.find('-avg') >= 0:
         # stage the full time array analysis results to the tmp dir
         # do not copy png images
         files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 or fn.find('.npy') > 0]
         utils.copy2scratch(files, gta.workdir)

         if args.state == 'fullflare-avg':
             gta.load_roi('avgspec') # reload the average spectral fit
         elif args.state == 'gtiflare-avg':
             gta.load_roi('fullflare-avg') # reload the average spectral fit

         o = gta.optimize() # perform an initial fit
         logging.debug(o)
         gta.print_roi()

         # Free all parameters of all Sources within X deg of ROI center
         #gta.free_sources(distance=fit_config['ps_dist_all'])
         # Free Normalization of all Sources within X deg of ROI center
         gta.free_sources(distance=fit_config['ps_dist_norm_fullflare'],pars=fa.allnorm)
         # Free spectra parameters of all Sources within X deg of ROI center
         gta.free_sources(distance=fit_config['ps_dist_idx_fullflare'],pars=fa.allidx)
         # Free all parameters of isotropic and galactic diffuse components
         gta.free_source('galdiff', pars= fa.allnorm, free = fit_config['gal_norm_free_fullflare'])
         gta.free_source('galdiff', pars= fa.allidx, free = fit_config['gal_idx_free_fullflare'])
         gta.free_source('isodiff', pars= fa.allnorm, free = fit_config['iso_norm_free_fullflare'])
         # Free sources with TS > X
         gta.free_sources(minmax_ts=[fit_config['ts_norm'],None], pars = fa.allnorm)
         # Fix sources with TS < Y
         gta.free_sources(minmax_ts=[None,fit_config['ts_fixed']],free=False,
                               pars=fa.allidx + fa.allnorm)
         # Fix indeces for sources with TS < Z
         gta.free_sources(minmax_ts=[None,fit_config['ts_fixed_idx']],free=False,
                  pars = fa.allidx)
         # Fix sources Npred < Z
         gta.free_sources(minmax_npred=[None,fit_config['npred_fixed']],free=False,
                               pars=fa.allidx + fa.allnorm)

         # gives "failed to create spline" in get_parameter_limits function
         #gta = fa.utils.add_ebl_atten(gta,config['selection']['target'],fit_config['z'])
         logging.info('free source parameters:')
         for s in  gta.get_sources() :                                                           
             for k in s.spectral_pars.keys():
                  if s.spectral_pars[k]['free']:   
                      logging.info('{0:s}: {1:s}'.format(s.name, k))
         f = gta.fit()

         retries = f['config']['retries']
         tsfix = fit_config['ts_fixed']
         while not f['fit_success'] and retries > 0:
             gta.free_source(config['selection']['target'], pars = ['beta', 'Index2'], free = False)
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

         gta.write_roi(args.state)

         logging.info('Running sed for {0[target]:s}'.format(config['selection']))
         sed = gta.sed(config['selection']['target'], prefix = args.state)
         model = {'Scale': 1000., 'Index' : fit_config['new_src_pl_index'], 'SpatialModel' : 'PointSource'}
         resid_maps = gta.residmap(args.state,model=model, make_plots=True, write_fits = True, write_npy = True)
    return f,gta

if __name__ == '__main__':
    gta = main()
