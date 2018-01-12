from haloanalysis.batchfarm import utils,lsf
from os import path
import os
import yaml
import copy
import logging 
from fermipy.gtanalysis import GTAnalysis

def init_gta(configfile, i = 1, logging_level = "INFO"):
    """
    Initialize the fermipy analysis
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
#	config[k] = utils.copy2scratch(config[k], tmpdir)
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

    # remove parameters from config file not accepted by fermipy
    for k in ['configname', 'tmp', 'log', 'fit_pars']: 
	config.pop(k,None)

    if config['data']['ltcube'] == '':
	config['data'].pop('ltcube',None)

    logging.info('Starting with fermipy analysis')

    if type(config['lightcurve']['binsz']) == str:
	config['lightcurve']['binsz'] = 3. * 3600.
    gta = GTAnalysis(config,logging={'verbosity' : 3})
    return gta, config, fit_config, job_id
