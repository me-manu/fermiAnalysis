"""
Helper functions of analysis preparation
"""
# --- Imports -------------- #
import logging
import time
from math import sqrt,floor
from os.path import *
# -------------------------- #

def set_params(gtobject,infile,outfile,config,params,**kwargs):
    """
    Set parameters of a gt analysis object

    Parameters
    ----------
    gtobject:	a GtApps object
    infile:	string, full path to input file
    outfile:	sting, full path to output file
    config:	dictionary containing parameters for gtobject
    params:	list of parameters that should be set

    Returns
    -------
    GtApps object with new set of parameters
    """
    if not exists(infile.lstrip('@')):
	logging.error('*** infile {0:s} does not exist! Returning -2'.format(infile))
	return -2

    for k in params:
	if k in ['infile','evfile','expcube','specfile']:
	    gtobject[k] = infile
	elif k == 'outfile':
	    gtobject[k] = outfile
	else:
	    if k in ['emin','emax','tmin','tmax']:
		gtobject[k] = config['selection'][k]
	    else:
		gtobject[k] = config[k]

    # --- Set hidden parameters
    for k in kwargs.keys():
	gtobject[k] = kwargs[k]

    return gtobject

def set_t_binning(gtobject, config):
    """
    Set the time binning for gtbin. S/N binning not implemented
    """
    if isinstance(config['lightcurve']['binsz'], str): 
	gtobject['tbinalg']	= "FILE"
	gtobject['tbinfile']	= config['lightcurve']['binsz']
    else:
	gtobject['tstart']	= config['tmin']
	gtobject['tstop']	= config['tmax']
	gtobject['tbinalg']	= "LIN"
	gtobject['dtime']	= config['lightcurve']['binsz']
    return gtobject

def set_e_binning(gtobject, config):
    """
    Set the energy binning for gtexpcube2 or gtbin
    """
    if isinstance(config['ebins'], int):
	gtobject['ebinalg']	= "LOG"
	gtobject['enumbins']	= config['ebins']
    elif isinstance(self.params['bins'], str):
	gtobject['ebinalg']	= "FILE"
	gtobject['ebinfile']	= config['ebins']
    return gtobject


def set_spatial_e_binning(gtobject,config,kwargs):
    """
    Set the spatial and energy binning for gtbin and gtexposure2

    Parameters
    ----------
    gtobject:	GtApps.GtApps object, either of gtbin or gtexpcube2
    config:	dictionary with config parameters
    kwargs:	dictionary with additional kwargs

    Returns
    -------
    the updated gtobject
    """

    if kwargs['coordsys'] == 'CEL':
	gtobject['xref'] = config['ra']
	gtobject['yref'] = config['dec']

    elif kwargs['coordsys'] == 'GAL':
	gtobject['xref'] = config['glon']
	gtobject['yref'] = config['glat']

    gtobject['binsz'] = kwargs['binsz']

    try:
	gtobject['nxpix'] = int(kwargs['nxpix'])
	gtobject['nypix'] = int(kwargs['nypix'])
    except KeyError:
	if gtobject.app == 'gtexpcube2':
	    r = config['exprad']
	elif gtobject.app == 'gtbin':
	    r = config['rad']
	gtobject['nxpix'] = int(floor(sqrt(2.) * r / 10.) * 10 / kwargs['binsz'])
	gtobject['nypix'] = int(floor(sqrt(2.) * r / 10.) * 10 / kwargs['binsz'])

    # --- Set the binning
    gtobject = set_e_binning(gtobject, config)

    return gtobject

def run_gtobject(gtobject, tries = 5):
    """
    Wrapper for running a gtobject command

    Parameters
    ----------
    gtobject:	a GtApp.GtApp object

    kwargs
    ------
    tries:	int, if outfile has size zero, retry [tries] number of times to create outfile (default: 5)

    Returns
    -------
    Integer:
	0: everything went fine
	-1: encountered a runtime error
	-3: outfile has size 0
    """
    logging.info("\nRunning {0} at {1} ...\n".format(gtobject.app,time.strftime("%d %b %Y %H:%M", time.gmtime())))
    logging.info("Parameters: {0}\n".format(gtobject.pars()))
    size,nTry = 0,0
    while( (not size) and nTry < tries ): 
	try:
	    gtobject.run()
	except RuntimeError as e:
	    logging.error("***Runtime error in {0}: {1}".format(gtobject.app,e.args[0]))
	    return -1
	nTry += 1
	size = getsize(gtobject['outfile'])

    logging.info("\nDone with {0} at {1} ...\n".format(gtobject.app,time.strftime("%d %b %Y %H:%M", time.gmtime())))
    if size:
	return 0
    else:
	return -3
