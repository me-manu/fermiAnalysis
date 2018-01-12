"""
utility functions 
"""

# --- Imports -------------- #
import sys
from glob import glob
from subprocess import call
from os.path import *
from os import environ
from math import ceil
from numpy import isscalar,array,where,sqrt,where,log10,linspace,cumsum,argmin,abs,empty,pi
import numpy as np
from fermiAnalysis.defaults import AnalysisDefaults as ad
from time import sleep
from numpy import string_
import logging
import cPickle
import gzip
from scipy.interpolate import UnivariateSpline as USpline
try:
    from skymaps import SkyDir
except ImportError as e:
    logging.error('{0}'.format(e))
# -------------------------- #

def save(object, filename, protocol = -1):
    """
    Save an object to a compressed disk file.
    Works well with huge objects.

    Parameters
    ----------
    object:	some python object
    filename:	str, full path to file

    kwargs
    ------
    protocol:	int, see Pickler documentation
    """
    file = gzip.GzipFile(filename, 'wb')
    cPickle.dump(object, file, protocol)
    file.close()

def load(filename):
    """
    Loads a compressed object from disk

    Parameters
    ----------
    filename:	str, full path to file

    Returns
    -------
    The loaded object.
    """
    file = gzip.GzipFile(filename, 'rb')
    object = cPickle.load(file)
    file.close()

    return object

def set_njobs(njobs, missing, irf = "None", clobber = False):
    """
    determine number of jobs

    Parameters
    ----------
    njobs:	int, total number of jobs
    missing:	list, list with job numbers

    kwargs
    ------
    irf:	str, irf name
    clobber: bool, if true, overwrite existing files

    Returns
    -------
    list or int with job ids to be submitted to lsf cluster
    """
    if len(missing) < njobs:
	if not len(missing):
	    logging.debug('all files present for irf {0:s}.'.format(irf))
	    if clobber:
		logging.warning('Clobber is set to yes, will overwrite files')
	    else:
		njobs = 0
	else:
	    njobs = missing
	    logging.info('there are {0:n} files missing'.format(len(missing)))

    logging.debug('missing: {0}'.format(missing))
    return njobs

def determine_simdir(simflag):
    """
    determine the directory of the simulations

    Parameters
    ----------
    simflag:	list or str, idenfitiers for simulations, e.g. cen_src or ['cen_src','ps','bkg']

    Returns:
    --------
    str, 	str, path to simulations
    """
    if type(simflag) == list:
	load_dir = ''
	for i,s in enumerate(simflag):
	    if not i:
		load_dir = s
	    else:
		load_dir += '+' + s
    elif type(simflag) == str:
	load_dir = simflag

    return load_dir


def sannity(kwargs):
    """
    Check options for their sannity. Raise ValueError if an error occurs.
    """
    for k in kwargs.keys():
	if k == 'coordsys':
	    if not (kwargs['coordsys'] == 'CEL' or kwargs['coordsys'] == 'GAL'):
		raise ValueError('Coordsys must be either or CEL or GAL, at the moment: {0:s}'.format(kwargs['coordsys']))
	if k == 'proj':
	    if not (kwargs['proj'] == 'AIT' or kwargs['proj'] == 'CAR'):
		raise ValueError('Projection must be either or AIT or CEL, at the moment: {0:s}'.format(kwargs['proj']))
	if k == 'ra' or k == 'dec' or k == 'rad' or k == 'exprad':
	    if not (isscalar(kwargs[k])):
		raise ValueError('{0:s} must be a scalar!'.format(k))
	if k == 'irfs':
	    if not type(kwargs[k]) == list:
		raise ValueError('irfs keyword must be a list of irfs!')
    return

def mkdir(directory):
    """
    Create a new directory

    Paramters
    ---------
    directory: string of new directory

    Returns
    -------
    directory: string of new directory
    """
    if not exists(directory):
	call(['mkdir','-p',directory])
	sleep(1.)
    return directory

def rm(path, dry = False):
    """
    Remove a file or a directory

    Parameters
    ----------
    path:	str, full path to file or directory that should be deleted

    kwargs
    ------
    dry:	bool, if true, do not remove but pring messages

    Returns
    -------
    int, if zero, deletion successfull.
    """

    if type(path) == str or type(path) == string_:
	if not exists(path):
	    logging.error("Path {0:s} does not exist".format(path))
	    return -1
	if isdir(path):
	    logging.warning("*** Deleting {0:s}".format(path))
	    sleep(0.5)
	    if not dry:
		call(['rm','-r',path])
	elif isfile(path):
	    logging.warning("*** Deleting {0:s}".format(path))
	    sleep(0.5)
	    if not dry:
		call(['rm',path])
	else:
	    logging.error("Path {0:s} does not point to a file or directory")
	    return -2

    elif type(path) == list:
	logging.warning("*** Deleting {0:s}".format(path))
	sleep(0.5)
	if isdir(path[0]):
	    if not dry:
		call(['rm','-r'] + path)
	elif isfile(path[0]):
	    if not dry:
		call(['rm'] + path)
    return 0

def update_diff(dict1,dict2):
    """
    Update dict1 with values of dict2 without adding additional keys from dict2 to dict1

    Parameters
    ----------
    dict1:	old dictionary
    dict2:	dictionary with updated values and possibly additional keys

    Returns
    -------
    updated dict1 and dict2 with updated items removed
    """
    diff_keys1	= set(dict2) - set(dict1)
    diff_keys2	= set(dict1) - set(dict2)
    dict1.update(dict2)

    for k in diff_key1s: del dict1[k]	# remove unwanted / additional keywords from dict1
    for k in diff_key2s: del dict2[k]	# remove unwanted / additional keywords from dict1

    return dict1,dict2

def set_defaults(kwargs,defaults):
    """
    Set default kwargs

    Parameters
    ----------
    kwargs:	dictionary with kwargs
    defaults:	dictionary with default kwargs

    Returns:
    --------
    dictionary with kwargs set to defaults if not given in input
    """
    for k in defaults.keys():
	kwargs.setdefault(k,defaults[k])

    try:
	defaults.keys().index('key')
	if (kwargs['key'] == 'SED' or kwargs['key'] == 'edisp'):
	    raise ValueError('*** key keyword has to be either edisp or SED not {0:s}'.format(kwargs['key']))
    except ValueError:
	pass

    try:
	defaults.keys().index('clobber')
	kwargs['clobber'] = clobber2Bool(kwargs['clobber'])
    except ValueError:
	pass

    return kwargs

def determine_tbins(dt,config):
    """
    Determine number of time bins 

    Parameters
    ----------
    config: dictionary containing the keys tmin and tmax in seconds
    dt:	float, time binning in seconds

    Returns
    -------
    integer with number of time bins (ceil)
    """
    return int(ceil((config['tmax'] - config['tmin']) / float(dt)))

def init_logging(level, color = False):
    """
    Setup logger.

    Parameters 
    ----------
    level:	string, level of logging: DEBUG,INFO,WARNING,ERROR. (default: INFO).

    kwargs
    ------
    color:	bool, if true, enable colored output for bash output

    Notes
    -----
    for color see
	stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
	https://wiki.archlinux.org/index.php/Color_Bash_Prompt
    """
    for handler in logging.root.handlers[:]:
	logging.root.removeHandler(handler)
    if level.upper() == 'INFO':
	level = logging.INFO
    elif level.upper() == 'DEBUG':
	level = logging.DEBUG
    elif level.upper() == 'WARNING':
	level = logging.WARNING
    elif level.upper() == 'ERROR':
	level = logging.ERROR


    if color:
	logging.basicConfig(level=level,stream = sys.stderr, format='\033[0;36m%(filename)10s:\033[0;35m%(lineno)4s\033[0;0m --- %(levelname)7s: %(message)s')
	logging.addLevelName( logging.DEBUG, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
	logging.addLevelName( logging.INFO, "\033[1;36m%s\033[1;0m" % logging.getLevelName(logging.INFO))
	logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
	logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))
    else:
	logging.basicConfig(level=level,stream = sys.stderr, format='%(filename)10s:%(lineno)4s --- %(levelname)7s: %(message)s')

    return

def set_tbins(nbin,config):
    """
    Change tmin and tmax to current bin (assuming linear binning)

    Parameters
    ----------
    nbin:	number of bin, must be > 0, if 0, set to one.
    config:	dict with config parameters - needs to contain tmin, tmax, and dt keywords

    Return
    ------
    tuple with new tmin and tmax for current bin
    """
    if not nbin:
	logging.warning('** time bin no. is 0, set to 1.')
	nbin = 1
    tmin = config['tmin'] + (nbin - 1) * config['dt']
    tmax = config['tmin'] + nbin * config['dt']
    if tmax > config['tmax']:
	tmax = config['tmax']
    return tmin,tmax

def set_ebins(nbin,config):
    """
    Change emin and emax to current bin (assuming log binning)

    Parameters
    ----------
    nbin:	number of bin, must be > 0, if 0, set to one.
    config:	dict with config parameters - needs to contain emin, emax, and ebins keywords

    Return
    ------
    tuple with new emin and emax for current bin
    """
    if not nbin:
	logging.warning('** energy bin no. is 0, set to 1.')
	nbin = 1
    ebin = 10.**linspace(log10(config['emin']),log10(config['emax']),config['ebins'] + 1)
    emin = ebin[nbin - 1]
    emax = ebin[nbin]
    if emax > config['emax']:
	tmax = config['emax']
    return emin,emax

from scipy.interpolate import interp1d
def calc_Eresol(edisp,etrue,ereco,axis = 0, idx = -1, conf = 0.68):
    """
    Calculate energy resolution from energy dispersion matrix

    Parameters
    ----------
    edisp: nxm dim np.array, energy dispersion matrix 
    etrue: n dim np.array, true energy values (bin centers)
    ereco: m dim np.array, reconstructed energy values (bin centers)

    kwargs
    ------
    axis: int, axis of true energy (default: 0)
    conf: float, interval around median for energy resolution (default: 0.68)
    idx: int, if >=0 : exit at this energy and return

    Return
    ------
    m dim np.array with energy resolution dE / E
    """
    # form the cumulative distripution along the true energy axis
    invert_axis = int(np.invert(bool(axis)))
# new version 
    idmax = np.argmax(edisp, axis = axis)
    edmax = np.max(edisp, axis = axis)

# old version 
    cdf = cumsum(edisp, axis = axis)
    dE = np.empty(edisp.shape[int(np.invert(bool(axis)))])

    cen = 0.5

    for j in range(edisp.shape[int(np.invert(bool(axis)))]):	# loop over ereco


	if not axis:
	    p = interp1d(cdf[:,j],etrue)
	    spline = USpline(np.log(etrue), edisp[:,j], s = 5e-5, ext = 0)
	else:
	    p = interp1d(cdf[j,:],etrue)
	    spline = USpline(np.log(etrue), edisp[j,:], s = 5e-5, ext = 0)

#	cen = spline.integral(np.log(etrue[0]), np.log(etrue[idmax[j]])) / spline.integral(np.log(etrue[0]), np.log(etrue[-1])) 
#	if cen <= 0. or cen >= 1. : cen = 0.5


	try:
	    Elo = p(cen - conf / 2.)
	except ValueError:
	    Elo = etrue[0]
	try:
	    Eup = p(cen + conf / 2.)
	except ValueError:
	    Eup = etrue[-1]
	#dE[j] = (Eup - Elo) / ereco[j]

	if idx >=0 and j == idx: 
	    return p(cen), Elo, Eup, etrue, p, ereco[j]

	dE[j] = (Eup - Elo) / p(cen)

    # compute the energies +/- 32 % around median and the width
    #ielo = argmin(abs(edisp - (np.diag(edisp[idmax,:]) - conf / 2.)), axis = axis)
    #ieup = argmin(abs(edisp - (np.diag(edisp[idmax,:]) + conf / 2.)), axis = axis)
    #ielo = argmin(abs(edisp - (np.diag(edisp[idmax,:]) - conf / 2.)), axis = axis)
    #ieup = argmin(abs(edisp - (np.diag(edisp[ielo,:]) + conf)), axis = axis)
    #print idmax
    #print ielo
    #print ieup
    #print etrue[ieup],etrue[idmax],etrue[ielo]
    #return (etrue[ieup] - etrue[ielo]) / ereco
    return dE


def missing_files(fname,fn, look = 'ltcube', missing = True, minimum = 0, nbins = 0):
    """
    Check for missing files.

    Parameters
    ----------
    fname:	string, file names (with wildcards), last 4 spaces need to be file numbers, e.g. file0001.file
    fn:		int, number of files that should be present

    kwargs
    ------
    minimum:	int, minimum job id that will be resubmitted (default: 0)
    missing:	bool, if True, look for missing files, if false, look for present files.
    look:	str, specify which files you are looking for. Possibilities:
    			- ltcube
			- mc
			- SED
			- FullSED
			(default: ltcube)
    nbins:	int, if look == FullSED, this number gives the number of bin files that should be present

    Returns
    -------
    list with missing file numbers
    """
    files = glob(fname)
    if look == 'ltcube' or look == 'SED':
	idxs  = array(map( lambda f: int(basename(f).split('.')[0][-4:]), files))
    elif look == 'mc':
	try:
	    idxs  = array(map( lambda f: int(f.split('/')[-2]), files))
	except ValueError:
	    idxs  = array(map( lambda f: int(f.split('/')[-3]), files))	# srcmap with additional folder in logL_prefix
    elif look == 'FullSED':
	idxs  = array(map( lambda f: int(f.split('/')[-3]), files))
    else:
	raise ValueError('*** Unknown option {0:s} for look kwargs'.format(look))

    miss = []
    logging.debug('number of files that should be there: {0:n}'.format(fn))
    if not fn == len(files):
	for idxi in range(1,fn + 1):
	    if missing:
		if look == 'FullSED' :
		    if not len(glob(fname.replace('*','{0:05n}'.format(idxi),1))) == nbins:
			miss.append(idxi)
		else:
		    if not len(where(idxi == idxs)[0]):
			miss.append(idxi)
	    else:
		if look == 'FullSED':
		    if len(glob(fname.replace('*','{0:05n}'.format(idxi),1))) == nbins:
			miss.append(idxi)
		else:
		    if len(where(idxi == idxs)[0]):
			miss.append(idxi)
    if fn == len(files) and not missing:
	miss = range(1,fn + 1)

    logging.debug('missing: {0}{1}, min: {2}'.format(missing,miss,minimum))

    if minimum:
	miss = np.array(miss)
	miss = miss[miss > minimum]
	logging.debug('2nd missing: {0}{1}, min: {2}'.format(missing,miss,minimum))
	return list(miss)

    else:
	return miss

def clobber2Bool(clobber):
    """
    Convert clobber string (yes / no) to boolean (1 / 0).
    """
    if type(clobber) == str:
	if clobber.lower() == 'no':	return False
	elif clobber.lower() == 'yes':	return True
    elif type(clobber) == bool:
	return clobber
    elif type(clobber) == int:
	if clobber == 0:	return False
	else: 			return True
    else:
	raise ValueError('Clobber must be either string, bool, or int, not {0} of type {1}'.format(clobber,type(clobber)))
    return

def filesExists(*files):
    """
    test if file(s) exists.

    Parameters
    ----------
    *files: list of file names

    Returns
    -------
    bool, 0 if one file is missing, 1 if all files are present.
    """
    res = 0
    for f in files:
	if exists(f):
	    logging.info('Found file {0}'.format(f))
	    res += 1
    if res == len(files): 
	return True
    else: 
	return False


def binCenter(binBounds):
    """
    Calculate central energies from an array of bin bounds for logarithmic spaced bins

    Parameters
    ----------
    binBounds:	n-dim np.array with bin bounds

    Returns
    -------
    (n-1)-dim np.array with bin centers
    """
    return sqrt(binBounds[:-1] * binBounds[1:])

def binWidth(binBounds):
    """
    Calculate bin width from an array of bin bounds for logarithmic spaced bins

    Parameters
    ----------
    binBounds:	n-dim np.array with bin bounds

    Returns
    -------
    (n-1)-dim np.array with bin centers
    """
    return binBounds[1:] / binBounds[:-1]

def binBounds(binCenter):
    """
    calculate bin boundaries for array assuming that array values lie at logarithmic bin center

    Parameters
    ----------
    binCenter:	n-dim array with logarithmic center values

    Returns
    -------
    (n+1) dim array with bin boundaries
    """
    bin_bounds = np.empty(binCenter.shape[0] + 1)
    for i,x in enumerate(binCenter[:-1]):
	bin_bounds[i + 1] = np.sqrt(x * binCenter[i + 1])
    bin_bounds[0]	= binCenter[0] **2. / bin_bounds[1]
    bin_bounds[-1]	= binCenter[-1]**2. / bin_bounds[-2]
    return bin_bounds


def angdist(ra,dec,ps,ds = None):
    """
    Calculate angular distance of point sources to point in the sky

    Parameters
    ----------
    ra:		float, R.A. of ROI (in degree)
    dec:	float, DEC. of ROI (in degree)
    ps:		list, list of point sources as generated with pointlike

    kwargs
    ------
    ds:		list, list of point sources as generated with pointlike, if None, don't use it and 
    		don't return it

    Returns
    -------
    list that contains the angular distance in degrees to point sources
    """

    skydir	= SkyDir(ra,dec, SkyDir.EQUATORIAL)
    diffPs	= empty(len(ps))

    for i,p in enumerate(ps):
	diffPs[i] = skydir.difference(p.skydir) * 180. / pi	# angular separation in degree
    if not ds == None:
	diffDs	= empty(len(ds))
	for i,d in enumerate(ds):
	    diffDs[i] = skydir.difference(d.skydir) * 180. / pi	# angular separation in degree
	return diffPs,diffDs
    else:
	return diffPs

def angsep(r1,r2,d1,d2):
    """Calculate angulate seperation in degrees for two (RA,DEC) tuples in degrees

    Parameters:
    -----------
    r1: float or n-dim array, RA of 1st source in degrees
    r2: float or n-dim array, DEC of 1st source in degrees
    d1: float or n-dim array, RA of 2nd source in degrees
    d2: float or n-dim array, DEC of 2nd source in degrees

    Returns
    -------
    angular seperation in degrees
    """

    if np.isscalar(r1):
	r1 = np.array([r1])
    if np.isscalar(r2):
	r2 = np.array([r2])
    if np.isscalar(d1):
	d1 = np.array([d1])
    if np.isscalar(d2):
	d2 = np.array([d2])
    th1 = (90.0*np.ones(d1.shape)-d1)*pi/180.
    th2 = (90.0*np.ones(d2.shape)-d2)*pi/180.  
    ph1 = pi/180. * r1  
    ph2 = pi/180. * r2  
    tt1,tt2 = np.meshgrid(th1,th2)
    pp1,pp2 = np.meshgrid(ph1,ph2)

    cth = np.sin(tt1)*np.sin(tt2)*\
	(np.cos(pp1)*np.cos(pp2)+np.sin(pp1)*np.sin(pp2))+np.cos(tt1)*np.cos(tt2)
    return np.arccos(cth)*180./pi * ( cth <= 1. )

def split_simflag(option, opt, value, parser):
    """
    Split command line argument simflag 
    """
    if value == None:
	return
    if len(value.split(',')) > 1:
	setattr(parser.values, option.dest, value.split(','))
    else:
	setattr(parser.values, option.dest, value)
    return

def check4events(eventfile, dataid = 1):
    """
    Check if FT1 event file contains any events

    Parameters
    ----------
    eventfile:	str, path to ft1 event file

    kwargs
    ------
    dataid:	int, index in fits file that contains the data (default = 1)

    Returns
    -------
    int, number of events in ft1 file
    """
    hdu = pyfits.open(eventfile)
    return hdu[dataid].data.field(0).shape[0]

def printOneLine(string):
    """
    Print an output on the same line

    Parameter
    ---------
    string:	str,the string to be printed
    """
    sys.stdout.write(string)
    sys.stdout.flush()
    lenStr = len(string)
    sys.stdout.write("\b"*lenStr)
    return

# To do:
#def make_ebin_fits_file():
