try:
    from fermipy.utils import init_matplotlib_backend
    init_matplotlib_backend()
except:
    pass
import logging
import argparse
import numpy as np
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
from fermiAnalysis import setup
from fermipy.gtanalysis import GTAnalysis
from fermiAnalysis.batchfarm import utils,lsf

if __name__ == "__main__":
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run gta setup and return the gta object for testing in ipython"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('--overwrite', required=False, default = 0, 
                        help='Overwrite existing files', type=int)
    parser.add_argument('--model')
    parser.add_argument('-i', required=False, default = 0, 
                        help='Set local or scratch calculation', type=int)
    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i = args.i, logging_level = "INFO")
    init_matplotlib_backend()

    files = [fn for fn in glob(fit_config['avgspec']) if fn.find('.xml') > 0 or fn.find('.npy') > 0]
    if len(files):
        utils.copy2scratch(files, gta.workdir)
    else:
        logging.error("No files found in {0:s}".format(fit_config['avgspec']))
        args.reloadfit = False

    if args.model is None:
        try:
            logging.info('Running fermipy setup')
            gta.setup()
        except RuntimeError as e:
            logging.error("setup ended with runtime error:\n{0}.".format(e))

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
    else:
        gta.load_roi(args.model)
