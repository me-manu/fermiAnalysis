from haloanalysis.batchfarm import utils,lsf
from os import path
import os
import yaml
import copy
import logging 
from fermipy.gtanalysis import GTAnalysis
from fermiAnalysis.utils import excise_solar_flares_grbs

def init_gta(configfile, i = 1, logging_level = "INFO", tsmin = 100):
    """
    Initialize the fermipy analysis

    Add filter to config expression excising the brightest GRBs and 
    solar flares
    """
    utils.init_logging(logging_level)
    config = yaml.load(open(configfile))
    tmpdir, job_id = lsf.init_lsf()
    if not job_id:
        job_id = i 
        tmpdir = os.environ["PWD"]
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
    logging.debug(config['selection'])
    config['fileio']['scratchdir'] = tmpdir

    # set the log file
    logdir = copy.deepcopy(config['fileio']['logfile'])
    config['fileio']['logfile'] = path.join(tmpdir,'fermipy.log')

    # copy the background files if we are on cluster
    #if not args.state == 'setup' and job_id:
        #for k in ['galdiff','isodiff']:
            #config['model'][k] = utils.copy2scratch(path.expandvars(config['model'][k]), tmpdir)

    # copy all fits files already present in outdir
    # run the analysis
    fit_config = copy.deepcopy(config['fit_pars'])

    # create a file with GTIs excluding solar flares and GRBs

    if 'gtiexclude' in fit_config.keys():
    # create fits file with gtis to be excluded
    # uses ftcreate which causes fermipy to crash with segmentation fault
        if not path.isfile(fit_config['gtiexclude']):
            fit_config['gtiexclude'] = excise_solar_flares_grbs(tsmin = tsmin, outdir = tmpdir)
        config['selection']['filter'] += " && gtifilter('{0:s}', START)".format(fit_config['gtiexclude'])
        #config['selection']['filter'] += " && gtifilter('{0:s}', STOP)".format(fit_config['gtiexclude'])


    # remove parameters from config file not accepted by fermipy
    for k in ['configname', 'tmp', 'log', 'fit_pars']: 
        config.pop(k,None)

    if config['data']['ltcube'] == '':
        config['data'].pop('ltcube',None)

    logging.info('Starting with fermipy analysis')

    if type(config['lightcurve']['binsz']) == str:
        config['lightcurve']['binsz'] = 3. * 3600.

    try:
        gta = GTAnalysis(config,logging={'verbosity' : 3})
    # exception could be caused by unknown source, 
    # remove target name and try again
    except Exception as e:
        logging.warning('Cought Exception {0}'.format(e))
        logging.warning('Trying to remove target and working with coordinates')
        if 'target' in config['selection'] and \
            ('ra' in config['selection'] or 'glon' in config['selection']):
            config['selection']['target'] = None
            gta = GTAnalysis(config,logging={'verbosity' : 3})
        else:
            raise Exception("No coordinates specified in config file selection")
    return gta, config, fit_config, job_id
