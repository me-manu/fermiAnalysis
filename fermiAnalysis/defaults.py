from os.path import join
from os import environ

class AnalysisDefaults(object):
    """
    Class that stores analysis defaults in dictionaries.
    """
    lsf_defaults = {
	'queue'	: 'time',
	'time'	: '04:00',
	'jname'	: 'lsf',
	'sleep'	: 10.,
	'lsb_steps'	: 1,
	'concurrent'	: 0,
	'dependency'	: None
    }

    gtselect_defaults = {
	'ra'	: None,
	'dec'	: None,
	'tmin'	: 239557417,	# fermi mission time of first photon detected in week 9 minus 1 in seconds
	'tmax'	: 428859819,	# 72 months after tmin, 2014-04-04 15:43:36 UTC
	'zmax'	: 90.,
	'zmin'	: 0.,
	'rad'   : None,
	'exprad': None,
	'evclass': None,
	'evtype': "INDEF",
	'phasemin': 0.,
	'phasemax': 1.
    }
    gtselect_hidden_defaults = {
	'chatter'	: 2,		# verbosity, 0,2,4
	'clobber'	: True,		# boolean! 
	'convtype'	: -1,		# 0: front, 1: back, -1: both
	#'evclsmin'	: "INDEF",
	#'evclsmax'	: "INDEF",
    }
    gtmktime_defaults = {
	'roicut': 	'yes',
	#'roicut': 	'no',
	#'filter': 	"DATA_QUAL>0 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52",
	'filter': 	"DATA_QUAL>0 && LAT_CONFIG==1",
    }
    gtmktime_hidden_defaults = {
	'chatter'	: 2,		# verbosity, 0,2,4
	'clobber'	: True		# boolean! 
    }
    gtbin_defaults = {
	'axisrot'	: 0.,
	'coordsys'	: 'CEL',
	'proj'		: 'AIT',
	'binsz'		: 0.2,
	'chatter'	: 2,		# verbosity, 0,2,4
	'clobber'	: True		# boolean! 
    }
    gtltcube_defaults = {
	'binsz'		: 1.,
	'phibins'	: 0,
	'dcostheta'	: 0.025,
	'chatter'	: 2,		# verbosity, 0,2,4
	'clobber'	: True		# boolean! 
    }
    gtpsf_defaults = {
	'nenergies': 20,
	'thetamax': 60., 
	'ntheta': 300
    }
    gtrspgen_defaults = {
	'respalg': 'PS',
	'thetacut': 60., 
	'dcostheta': 0.5
    }
    gtexpcube2_defaults = {
#	'cmap'		: 'none', # do not use output of ccube here in order to accomodate contribution outside the ROI
#				 # see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/binned_likelihood_tutorial.html#computeExposure
	'axisrot'	: 0.,
	'coordsys'	: 'CEL',
	'proj'		: 'AIT',
	'binsz'		: 0.2,		# should be the same as in gtbin_defaults
	'chatter'	: 2,		# verbosity, 0,2,4
	'clobber'	: True		# boolean! 
    }
    gtsrcmaps_defaults = {
	'convol'	: True,		# convolve with psf
	'ptsrc'		: True,		# compute src maps for point sources
	'psfcorr'	: True,		# apply psf correction
	'chatter'	: 2,		# verbosity, 0,2,4
	'clobber'	: True		# boolean! 
    }

def set_default_lc(func):
    """
    Decorator for the Lightcurve class
    """
    def init(*args,**kwargs):
        kwargs.setdefault('dt',None)
	for d in [AnalysisDefaults.lsf_defaults,
		AnalysisDefaults.makeDefaults,
		AnalysisDefaults.gtrspgen_defaults
	    ]:
	    for k in d.keys():
		kwargs.setdefault(k,d[k])
        return func(*args, **kwargs)
    return init
