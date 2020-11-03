#from fermipy.skymap import Map
try:
    from fermipy.utils import init_matplotlib_backend
    init_matplotlib_backend()
except:
    pass
import argparse
import os
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
from fermiAnalysis import setup
from fermiAnalysis.batchfarm import utils,lsf
from fermiAnalysis.ext_funcs import fit_region, fit_halo, fit_halo_scan

def main():
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('--state', default = ['avgspec'],
                        choices = ['avgspec','avgspec_ebl'],
                            help='Analysis state')
    parser.add_argument('-i', required=False, default = 0, 
                        help='Set local or scratch calculation', type=int)
    parser.add_argument('--create_ts_maps', required=False, default = 1, 
                        help='Generate TS maps', type=int)
    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i = args.i, logging_level = "INFO")
    gta.logger.info('Running fermipy setup')

    init_matplotlib_backend()

    #files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 or fn.find('.npy') > 0]
    #if len(files):
    #    utils.copy2scratch(files, gta.workdir)
    #else:
    #    gta.logger.error("No files found in {0:s}".format(fit_config['avgspec']))

    gta.setup()

    gta.logger.info('reloading {0:s}'.format(args.state))
    gta.load_roi(args.state) # reload the average spectral fit

    modelname = "{0:s}".format(args.state)

    # change the outdir
    # not the greatest that I'm not using the API here, 
    # but no other possibility
    gta._outdir = os.path.join(gta.outdir, 'extension_' + modelname + '/')
    if not os.path.exists(gta.outdir):
        os.makedirs(gta.outdir)
    gta.logger.info("Set new outdir: {0:s}".format(gta.outdir))

    free_radius_sed = fit_config.get('extension', dict(free_radius_sed=1.)).pop('free_radius_sed', 1.)
    force_ps = fit_config.get('extension', dict(force_ps=False)).pop('force_ps', False)
    distance_free_norm = fit_config.get('extension', dict(distance_free_norm=1.5)).pop('distance_free_norm', 1.5)
    distance_free_shape = fit_config.get('extension', dict(distance_free_shape=1.)).pop('distance_free_shape', 1.)
    halo_fit = fit_config.get('extension', dict(halo_fit=False)).pop('halo_fit', False)
    halo_scan = fit_config.get('extension', dict(halo_scan=False)).pop('halo_scan', False)
    fit_halo_kwargs = fit_config.get('extension', dict(fit_halo_kwargs={})).pop('fit_halo_kwargs', {})
    scan_halo_kwargs = fit_config.get('extension', dict(scan_halo_kwargs={})).pop('scan_halo_kwargs', {})


    fit_region(gta, modelname, gta.config['selection']['target'],
               loge_bounds=None, 
               skip_opt=list(fit_config.get('fix_sources', {}).keys()),
               shape_ts_threshold=9.0,
               force_ps=force_ps,
               create_maps=args.create_ts_maps, 
               create_sed=False,
               free_radius_sed=free_radius_sed,
               distance_free_norm=distance_free_norm,
               distance_free_shape=distance_free_shape,
               **fit_config.get('extension', {})
               )

    if halo_fit:
        fit_halo(gta, modelname, gta.config['selection']['target'],
                 **fit_halo_kwargs
                )

    if halo_scan:
        fit_halo_scan(gta, modelname, gta.config['selection']['target'],
                 **fit_halo_kwargs
                )

    return gta, fit_config

if __name__ == '__main__':
    res = main()
