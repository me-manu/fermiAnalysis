# --- Imports ----------------------------------------------------- #
from specFeat.base.analysis import PointAnalysis
from specFeat.utils.tools import *
from specFeat.utils.lsf import *
from specFeat.base.defaults import AnalysisDefaults as ad, set_default_lc
from specFeat.base.prepare import PreparePointSource as pps
from os.path import *
import pyfits
# ----------------------------------------------------------------- #

class Lightcurve(PointAnalysis):
    """
    Class for generating light curve analysis files
    """
    @set_default_lc
    def __init__(self, *args, **kwargs):
	"""
	Initialize the class
	"""
	kwargs
	super(Lightcurve,self).__init__(*args, **kwargs)
	self.gti = None
	for k in ['SED','edisp']:
	    kwargs['key'] = k
	    self.__setup_t_outfiles(**kwargs)		# setup the outfile names
	return

    def __setup_t_outfiles(self, **kwargs):
	"""
	set the file names for the lightcurve files
	"""
	key = kwargs['key']
	self.lcdirs = {key : {}}
	for d in self.config['irfs']:
	    if self.gti == None:
		self.lcdirs[key][d]			= mkdir(join(self.savedirs[key][d], 'lc/full/'))
	    else:
		self.lcdirs[key][d]			= mkdir(join(self.savedirs[key][d], 
								'lc/gti{0:05n}/'.format(self.gti)))
	    if kwargs['dt'] == None:
		self.lcdirs[key][d]			= mkdir(join(self.lcdirs[key][d], 
								'dt_full/'))
	    else:
		self.lcdirs[key][d]			= mkdir(join(self.lcdirs[key][d], 
								'dt_{0[dt]:.0f}/'.format(kwargs)))
	    self.outfiles[key][d]['ft1lc']	= join(self.lcdirs[key][d],'ft1lc.fits')
	    self.outfiles[key][d]['rsp']	= join(self.lcdirs[key][d],'rsp.fits')
	    self.outfiles[key][d]['psf']	= join(self.lcdirs[key][d],'psf.fits')
	return

    def set_gti(self,fits, gti, **kwargs):
	"""
	Restrict analysis to one GTI

	Parameters
	----------
	fits:	string, fits file whose last extension contains the GTIs
	gti:	int, the GTI that should be used
	"""
	f = pyfits.open(fits)
	self.config['tmin'] = float(f[-1].data.field(0)[gti])
	self.config['tmax'] = float(f[-1].data.field(1)[gti])
	self.__dict__.update()
	self.gti = gti
	logging.info('Set time range to {0[tmin]:.2f}-{0[tmax]:.2f}'.format(self.config))
	for k in ['SED','edisp']:
	    kwargs['key'] = k
	    self.__setup_t_outfiles(**kwargs)		# setup the outfile names
	return

    @set_default_lc
    def makeFT1LC(self, **kwargs):
	"""
	Make a light curve from an ft1 eventfile by sending it
	it to the lsf cluster.

	kwargs
	------
	dt:		float, length of time interval in seconds (default: None, use full time)

	The remaining kwargs are the standard lsf cluster keywords.
	"""
	# --- determine the number of time bins
	if kwargs['dt'] == None:
	    kwargs['dt'] = self.config['tmax'] - self.config['tmin']

	# --- Loop over the IRFs
	for i,irf in enumerate(self.config['irfs']):
	    tmpconf,tmpdir,logdir = self.makeTmpConf(
				    'ft1lc{0:s}'.format(kwargs['key']),# a name for the logdir
				    # the savedir directory
				    dirname(self.outfiles[kwargs['key']][irf]['ft1lc']),
				    irf,	# the irf name, anyway replaced by CALDB
				    int(np.sum(self.config['evtype']).astype(int)) \
					if self.stVersion.find('10') == 0 else "INDEF",	# the event type
				    # the input file
				    self.config['ft1'] if self.config['select'] == None else self.config['select'],
				    self.config['ft2'],		# the spacecraft file
				    key = kwargs['key']
				    )

	    # --- extra keyword required for lc preparation
	    tmpconf['dtime']		= kwargs['dt']
	    tmpconf['tstart']		= float(self.config['tmin'])
	    tmpconf['tstop']		= float(self.config['tmax'])
	    logging.info("time window: {0[tstart]:.1f} - {0[tstop]:.1f}, dt: {0[dtime]:.1f}".format(tmpconf))

	    if not len(glob(self.outfiles[kwargs['key']][irf]['ft1lc'])) or kwargs['clobber'] == 'yes':

		kwargs['jname'] = 'ft1lc{0:s}'.format(kwargs['key'])	# job name
		submit_lsf(join(self.basedir,pps.fname),	# the script that will be called on cluster
			tmpdir,				# tmp dir to store bash script
			logdir,				# the log directory
			tmpconf,			# the config dict
			'-r ft1lc',			# additional option string
			#10,# testing
			1,				# the number of jobs to be submitted
			**kwargs
			)
	    else:
		logging.info('ft1 lc file for irf {0:s} present and clobber is set to False.'.format(irf))
	return 0

    @set_default_lc
    def makeRsp(self, **kwargs):
	"""
	Make the diffuse response and pha1 file for the observation with possibility to divide the time interval 
	and sending each interval to the lsf cluster.

	kwargs
	------
	dt:		float, length of time interval in seconds (default: None, use full time)

	The remaining kwargs are the standard lsf cluster keywords.
	"""
	# --- determine the number of time bins
	if kwargs['dt'] == None:
	    kwargs['dt'] = self.config['tmax'] - self.config['tmin']
	    optStr = '-r rspgen'
	else:
	    optStr = '-r rspgen_lc'
	nt = determine_tbins(kwargs['dt'],self.config)

	# --- Loop over the IRFs
	for i,irf in enumerate(self.config['irfs']):
	    tmpconf,tmpdir,logdir = self.makeTmpConf(
				    'rspgen{0:s}'.format(kwargs['key']),# a name for the logdir
				    # the savedir directory
				    dirname(self.outfiles[kwargs['key']][irf]['rsp']),
				    irf,	# the irf name, anyway replaced by CALDB
				    int(np.sum(self.config['evtype']).astype(int)) \
					if self.stVersion.find('10') == 0 else "INDEF",	# the event type
				    # the input file
				    self.config['ft1'] if self.config['select'] == None else self.config['select'],
				    self.config['ft2'],		# the spacecraft file
				    key = kwargs['key']
				    )

	    # --- extra keyword required for Rsp Gen
	    tmpconf['dt']		= kwargs['dt']
	    for p in ['respalg','thetacut','dcostheta']:
		tmpconf[p] = kwargs[p]
	    tmpconf['ebinalg']		= "LOG"
	    tmpconf['emin'] = self.config[kwargs['key']]['emin'] * 3. / 5.
	    tmpconf['emax'] = self.config[kwargs['key']]['emax'] * 2.
	    tmpconf['enumbins'] = self.config[kwargs['key']]['ebins'] * 4
	    tmpconf['irfs'] = irf # set the irf explicitely

	    # --- check for missing files
	    rspFiles	= join(tmpconf['savedir'],'rsp*.fits')
	    if len(glob(rspFiles)):
		missing	= missing_files(rspFiles,nt,look = 'ltcube')
	    else:
		missing = list(np.ones(nt))

	    if len(missing) < nt:
		njobs = missing
		logging.info('there are {0:n} files missing in {1:s}'.format(len(missing),rspFiles))
	    else:
		njobs = nt

	    if len(missing) or kwargs['clobber'] == 'yes':
		if not len(missing): njobs = nt

		kwargs['jname'] = 'rg{0:s}'.format(kwargs['key'])	# job name
		submit_lsf(join(self.basedir,pps.fname),	# the script that will be called on cluster
			tmpdir,				# tmp dir to store bash script
			logdir,				# the log directory
			tmpconf,			# the config dict
			optStr,				# additional option string
			njobs,				# the number of jobs to be submitted
			**kwargs
			)
	    else:
		logging.info('All response files present for irf {0:s} and clobber is set to False.'.format(irf))
	return

    @set_default_lc
    def makePSF(self, **kwargs):
	"""
	Make the psf for the observation with possibility to divide the time interval 
	and sending each interval to the lsf cluster.

	kwargs
	------
	dt:		float, length of time interval in seconds (default: None, use full time)

	The remaining kwargs are the standard lsf cluster keywords.
	"""
	# --- determine the number of time bins
	if kwargs['dt'] == None:
	    kwargs['dt'] = self.config['tmax'] - self.config['tmin']
	nt = determine_tbins(kwargs['dt'],self.config)

	# --- Loop over the IRFs
	for i,irf in enumerate(self.config['irfs']):
	    tmpconf,tmpdir,logdir = self.makeTmpConf(
				    'psf{0:s}'.format(kwargs['key']),# a name for the logdir
				    # the savedir directory
				    dirname(self.outfiles[kwargs['key']][self.config['irfs'][0]]['psf']),
				    irf,	# the irf name, anyway replaced by CALDB
				    int(np.sum(self.config['evtype']).astype(int)) \
					if self.stVersion.find('10') == 0 else "INDEF",	# the event type
				    # the input file
				    self.config['ft1'] if self.config['select'] == None else self.config['select'],
				    self.config['ft2'],		# the spacecraft file
				    key = kwargs['key']
				    )

	    # --- extra keyword required for LTCube production
	    tmpconf['dt']		= kwargs['dt']
	    tmpconf['irfs']		= irf	# needs to be set explicitedly for psf

	    # --- check for missing files
	    psfFiles	= join(tmpconf['savedir'],'psf*.fits')
	    if len(glob(psfFiles)):
		missing	= missing_files(psfFiles,nt,look = 'ltcube')
	    else:
		missing = list(np.ones(nt))

	    if len(missing) < nt:
		njobs = missing
		logging.info('there are {0:n} files missing in {1:s}'.format(len(missing),psfFiles))
	    else:
		njobs = nt

	    if len(missing) or kwargs['clobber'].lower() == 'yes':
		if not len(missing): njobs = nt

		kwargs['jname'] = 'psf{0:s}'.format(kwargs['key'])	# job name
		submit_lsf(join(self.basedir,pps.fname),	# the script that will be called on cluster
			tmpdir,				# tmp dir to store bash script
			logdir,				# the log directory
			tmpconf,			# the config dict
			'-r psf',			# additional option string
			#10,# testing
			njobs,				# the number of jobs to be submitted
			**kwargs
			)
	    else:
		logging.info('All psf files present for irf {0:s} and clobber is set to False.'.format(irf))
	return
