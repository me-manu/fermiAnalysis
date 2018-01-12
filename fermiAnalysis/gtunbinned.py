"""
Functions to run prepratory steps for Fermi analysis - wrapper to gt_apps and GtApps
"""

# --- Imports ------------------------------------------- #
import logging
try:
    import gt_apps as gt
except ImportError as e:
    logging.error(e)
import sys
import yaml
from os import chdir,environ, path
import os
from subprocess import call
from optparse import OptionParser
from glob import glob
from fermiAnalysis.defaults import AnalysisDefaults as ad
from fermiAnalysis.tools import *
from fermiAnalysis.gtprep import *
from haloanalysis.batchfarm.utils import copy2scratch
from haloanalysis.batchfarm.lsf import init_lsf
from numpy import ndarray,savetxt
from fermipy.gtanalysis import GTAnalysis
from fermipy.ltcube import LTCube
from fermipy.gtanalysis import run_gtapp
import fermipy.irfs as irfs
import subprocess
from fermipy.skymap import Map, HpxMap
# ------------------------------------------------------- #

# === the class ============================================================ #
class GTUnbinned(GTAnalysis):
    """ 
    Class that contains functions for the preparation of data products 
    of the Fermi point source analysis in the unbinned case
    """
    fname = __module__.replace('.','/') + '.py'

    def __init__(self, *args, **kwargs):
    	super(GTUnbinned,self).__init__(*args, **kwargs)
	return

    def setup(self, phibins = 0, force_create_srcmap=False, overwrite=False, **kwargs):
	"""
	Same as GTBinnedAnalysis.setup function, which does not create 
	srcmap, running the setup for each component

	{options} 

	overwrite : bool

	Run all pre-processing steps even if the output file of
	that step is present in the working directory.

	phibins : int
	number of phibins for ltcube (default: 0; averaged ltcube)
	"""

	for c in self.components:
	    loglevel = kwargs.get('loglevel', c.loglevel)
												    
	    c.logger.log(loglevel, 'Running setup for component %s',c.name)
	    #use_external_srcmap = c.config['gtlike']['use_external_srcmap']
														    
	    # Run data selection
	    #if not use_external_srcmap:
		#c._select_data(overwrite=overwrite, **kwargs)
	    # Create LT Cube
	    if c._ext_ltcube is not None:
		c.logger.log(loglevel, 'Using external LT cube.')
	    else:
		c._create_ltcube(phibins = phibins, overwrite=overwrite, **kwargs)

	    c.logger.debug('Loading LT Cube %s', c.files['ltcube'])
	    c._ltc = LTCube.create(c.files['ltcube'])
	    
	    # Extract tmin, tmax from LT cube
	    c._tmin = c._ltc.tstart
	    c._tmax = c._ltc.tstop
	    c.logger.debug('Creating PSF model')
	    c._psf = irfs.PSFModel.create(c.roi.skydir, c._ltc,
		c.config['gtlike']['irfs'],
		c.config['selection']['evtype'],
		c.energies)

	    # we don't need this for unbinned analysis
	    # Bin data and create exposure cube
	    #if not use_external_srcmap:
		#c._bin_data(overwrite=overwrite, **kwargs)
		#c._create_expcube(overwrite=overwrite, **kwargs)
		#c._bexp = Map.create_from_fits(c.files['bexpmap'])

	    # compute exposure 
	    c._create_expmap(overwrite=overwrite, **kwargs)

	    # Make spatial maps for extended sources
	    for s in c.roi.sources:
		if s.diffuse:
		    continue
		if not s.extended:
		    continue
		c.make_template(s, c.config['file_suffix'])
	    # Write ROI XML
	    c.roi.write_xml(c.files['srcmdl'])

	    # compute diffuse response 
	    c._diffrsp(c.files['srcmdl'], overwrite=overwrite, **kwargs)

	    # we don't need this for unbinned
	    # Create source maps file
	    #if force_create_srcmap:
	#	if not use_external_srcmap:
	#	    c._create_srcmaps(overwrite=overwrite)
	    if not c.config['data']['cacheft1'] and os.path.isfile(c.files['ft1']):
		c.logger.debug('Deleting FT1 file.')
		os.remove(c.files['ft1'])

	    c.logger.log(loglevel, 'Finished setup for component %s',
		c.name)
	return


    def _bin_data_lc(self,overwrite=False, dtime = 0., **kwargs):
	"""
	Run gtbin for a light curve counts map 
	"""
	if dtime > 0.:
	    self.config['lightcurve']['binsz'] = dtime
	logging.info("Binning for LC: {0}".format(self.config['lightcurve']['binsz']))

	loglevel = kwargs.get('loglevel', self.loglevel)

	for i,c in enumerate(self.components):
	    self.components[i]._files['lcmap'] = path.join(self.workdir, 
			    'lcmap{0[file_suffix]:s}.fits'.format(c.config))

	    kw = dict(algorithm='lc',
		evfile=c.files['ft1'],
		outfile=c.files['lcmap'],
		scfile=c.data_files['scfile'],
		emin = c.config['selection']['emin'],
		emax = c.config['selection']['emax'],
		tstart = c.config['selection']['tmin'],
		tstop = c.config['selection']['tmax'],
		chatter=self.config['logging']['chatter'])

	    if isinstance(self.config['lightcurve']['binsz'], str): 
		kw['tbinalg']= 'FILE'
		kw['tbinfile'] = self.config['lightcurve']['binsz'],
	    else:
		kw['tbinalg']= 'LIN'
		kw['dtime'] = self.config['lightcurve']['binsz']

	    if not os.path.isfile(c.files['lcmap']) or overwrite:
		run_gtapp('gtbin', self.logger, kw, loglevel=loglevel)
	    else:
		self.logger.debug('Skipping gtbin.')
	return 

    def _compute_lc_exp(self,srcmdl = 'none', target = "", specin = -2.1, overwrite=False, **kwargs):
	"""
	Run gtexposure to compute the exposure for a light curve
	"""


	loglevel = kwargs.get('loglevel', self.loglevel)

	for i,c in enumerate(self.components):

	    kw = dict(infile=c.files['lcmap'],
		scfile=c.data_files['scfile'],
		irfs = c.config['gtlike']['irfs'],
		srcmdl = path.join(self.workdir,'{1:s}_{0:02n}.xml'.format(i,srcmdl)) \
			    if not target == "" else "none",
		target = target, 
		specin = specin,
		emin = c.config['selection']['emin'],
		emax = c.config['selection']['emax'],
		enumbins = self.enumbins
		)

	    print kw['target'], kw['srcmdl']
	    run_gtapp('gtexposure', self.logger, kw, loglevel=loglevel)
	return 

    def _compute_diffrsp(self,srcmdl, overwrite=False, **kwargs):
	"""
	Run gtsrcprob on an ft1 file
	"""


	loglevel = kwargs.get('loglevel', self.loglevel)

	for i,c in enumerate(self.components):

	    kw = dict(evfile=c.files['ft1'],
		scfile=c.data_files['scfile'],
		irfs = c.config['gtlike']['irfs'],
		evtype = c.config['selection']['evtype'],
		srcmdl = path.join(self.workdir,'{1:s}_{0:02n}.xml'.format(i,srcmdl))
		)

	    logging.info("Using srcmdl {0:s}".format(kw['srcmdl']))
	    run_gtapp('gtdiffrsp', self.logger, kw, loglevel=loglevel)
	return 
    def _compute_srcprob(self,srcmdl, overwrite=False, **kwargs):
	"""
	Run gtsrcprob on an ft1 file
	"""


	loglevel = kwargs.get('loglevel', self.loglevel)

	for i,c in enumerate(self.components):
	    self.components[i]._files['ft1srcprob'] = path.join(self.workdir, 
			    'ft1srcprob{0[file_suffix]:s}.fits'.format(c.config))

	    kw = dict(evfile=c.files['ft1'],
		scfile=c.data_files['scfile'],
		outfile=c.files['ft1srcprob'],
		irfs = c.config['gtlike']['irfs'],
		srcmdl = path.join(self.workdir,'{1:s}_{0:02n}.xml'.format(i,srcmdl))
		)

	    logging.info("Using srcmdl {0:s}".format(kw['srcmdl']))
	    cmd = ['time','-p','gtsrcprob']
	    for k,v in kw.items():
		cmd.append('{0:s}={1:s}'.format(k,v))
	    cmd.append("clobber='yes'")
	    logging.info(' '.join(cmd))
	    os.system(' '.join(cmd))
	    #child = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
	    #while True:
	#	out = child.stderr.read(1)
	#	if out == '' and child.poll() != None:
	#	    break
	#	if out != '':
	#	    logging.info(out)
	#	    sys.stdout.write(out)
	#	    sys.stdout.flush()
	    #run_gtapp('gtsrcprob', self.logger, kw, loglevel=loglevel)
	return 

