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
from fermiAnalysis.batchfarm.utils import copy2scratch
from fermiAnalysis.batchfarm.lsf import init_lsf
from numpy import ndarray,savetxt
from fermipy.gtanalysis import GTAnalysis
from fermipy.ltcube import LTCube
from fermipy.gtanalysis import run_gtapp
import fermipy.irfs as irfs
import subprocess
from fermipy.skymap import Map, HpxMap
# ------------------------------------------------------- #

# === the class ============================================================ #
class PreparePointSource(GTAnalysis):
    """ 
    Class that contains functions for the preparation of data products 
    of the Fermi point source analysis
    """
    fname = __module__.replace('.','/') + '.py'

    def __init__(self, *args, **kwargs):
        super(PreparePointSource,self).__init__(*args, **kwargs)
        return

    def setup(self, force_create_srcmap=False, overwrite=False, **kwargs):
        """
        Same as GTBinnedAnalysis.setup function, which does not create 
        srcmap, running the setup for each component

        Parameters
        ----------
        overwrite : bool

        Run all pre-processing steps even if the output file of
        that step is present in the working directory.
        """

        for c in self.components:
            loglevel = kwargs.get('loglevel', c.loglevel)
                                                                                                    
            c.logger.log(loglevel, 'Running setup for component %s',c.name)
            use_external_srcmap = c.config['gtlike']['use_external_srcmap']
                                                                                                                    
            # Run data selection
            if not use_external_srcmap:
                c._select_data(overwrite=overwrite, **kwargs)
            # Create LT Cube
            if c._ext_ltcube is not None:
                c.logger.log(loglevel, 'Using external LT cube.')
            else:
                c._create_ltcube(overwrite=overwrite, **kwargs)

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
            # Bin data and create exposure cube
            if not use_external_srcmap:
                c._bin_data(overwrite=overwrite, **kwargs)
                c._create_expcube(overwrite=overwrite, **kwargs)
                c._bexp = Map.create_from_fits(c.files['bexpmap'])

            # Make spatial maps for extended sources
            for s in c.roi.sources:
                if s.diffuse:
                    continue
                if not s.extended:
                    continue
                c.make_template(s, c.config['file_suffix'])
            # Write ROI XML
            c.roi.write_xml(c.files['srcmdl'])

            # Create source maps file
            if force_create_srcmap:
                if not use_external_srcmap:
                    c._create_srcmaps(overwrite=overwrite)
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

    def _compute_lc_exp(self,srcmdl = 'none',
                        target = "", specin = -2.1, overwrite=False, **kwargs):
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
        #        out = child.stderr.read(1)
        #        if out == '' and child.poll() != None:
        #            break
        #        if out != '':
        #            logging.info(out)
        #            sys.stdout.write(out)
        #            sys.stdout.flush()
            #run_gtapp('gtsrcprob', self.logger, kw, loglevel=loglevel)
        return 

    def _compute_psf(self,thetamax=30., ntheta=300, overwrite=False, **kwargs):
        """
        Run gtpsf
        """
        ra = kwargs.pop('ra',c.config['selection']['ra'])
        dec = kwargs.pop('dec',c.config['selection']['dec'])

        loglevel = kwargs.get('loglevel', self.loglevel)

        for i,c in enumerate(self.components):
            self.components[i]._files['psf'] = path.join(self.workdir, 
                            'psf{0[file_suffix]:s}.fits'.format(c.config))

            kw = dict(expcube=c.files['ltcube'],
                outfile=c.files['psf'],
                irfs = c.config['gtlike']['irfs'],
                evtype = c.config['selection']['evtype'],
                ra = ra,
                dec = dec,
                emin =  c.config['selection']['emin'],
                emax =  c.config['selection']['emax'],
                nenergies = c.enumbins * 2,
                thetamax = thetamax,
                ntheta = ntheta,
                )

            cmd = ['time','-p','gtpsf']
            for k,v in kw.items():
                cmd.append('{0:s}={1:s}'.format(k,v))
            cmd.append("clobber='yes'")
            logging.info(' '.join(cmd))
            os.system(' '.join(cmd))
            #child = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            #while True:
        #        out = child.stderr.read(1)
        #        if out == '' and child.poll() != None:
        #            break
        #        if out != '':
        #            logging.info(out)
        #            sys.stdout.write(out)
        #            sys.stdout.flush()
            #run_gtapp('gtsrcprob', self.logger, kw, loglevel=loglevel)
        return 

    def _compute_edisp(self,thetamax=30., ntheta=300, overwrite=False, **kwargs):
        """
        Run gtedisp
        """
        ra = kwargs.pop('ra',c.config['selection']['ra'])
        dec = kwargs.pop('dec',c.config['selection']['dec'])

        loglevel = kwargs.get('loglevel', self.loglevel)

        for i,c in enumerate(self.components):
            self.components[i]._files['psf'] = path.join(self.workdir, 
                            'psf{0[file_suffix]:s}.fits'.format(c.config))

            kw = dict(expcube=c.files['ltcube'],
                outfile=c.files['psf'],
                irfs = c.config['gtlike']['irfs'],
                evtype = c.config['selection']['evtype'],
                ra = ra,
                dec = dec,
                emin =  c.config['selection']['emin'],
                emax =  c.config['selection']['emax'],
                nenergies = c.enumbins * 2,
                thetamax = thetamax,
                ntheta = ntheta,
                )

            cmd = ['time','-p','gtpsf']
            for k,v in kw.items():
                cmd.append('{0:s}={1:s}'.format(k,v))
            cmd.append("clobber='yes'")
            logging.info(' '.join(cmd))
            os.system(' '.join(cmd))
            #child = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            #while True:
        #        out = child.stderr.read(1)
        #        if out == '' and child.poll() != None:
        #            break
        #        if out != '':
        #            logging.info(out)
        #            sys.stdout.write(out)
        #            sys.stdout.flush()
            #run_gtapp('gtsrcprob', self.logger, kw, loglevel=loglevel)
        return 

    def run_rspgen(self,infile,outfile,config, **kwargs):
        """
        Run gtrspgen to create the response file

        Parameters
        ----------
        infile:        string, full path to input file (output of gtbin with algorithm PHA1)
        outfile:        string, full path to output file
        config:        dictionary with main analysis parameters:
                    'scfile','ebins','rad'
        kwargs
        ------
        Correspond to parameters which are usually not changed. See Notes and specFeat.base.AnalysisDefaults. 

        Returns
        -------
        boolean with error codes:
            0: everything went ok
            -1: Runtime error

        Notes
        -----
        For extensive help see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtbin.txt
        """
        # ---- set defaults:
        kwargs = set_defaults(kwargs,ad.gtrspgen_defaults)

        # ---- Set parameters:
        params        = ['specfile','outfile','scfile','emin','emax','ebinalg','enumbins','irfs']
        gtrspgen= set_params(gt.rspgen,infile,outfile,config,params, **kwargs)

        if gtrspgen == -2:
            return gtrspgen

        return run_gtobject(gtrspgen)

    def run_ltcube(self,infile,outfile,config, **kwargs):
        """
        Run gtltcube for the liftime exposure cube

        Parameters
        ----------
        infile:        string, full path to input file (output of gtmaketime)
        outfile:        string, full path to output file
        config:        dictionary with main analysis parameters:
                    'scfile'

        kwargs
        ------
        Correspond to parameters which are usually not changed. See Notes and specFeat.base.AnalysisDefaults. 

        Returns
        -------
        boolean with error codes:
            0: everything went ok
            -1: Runtime error
            -2: Input file not found

        Notes
        -----
        For extensive help see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtltcube.txt
        """
        # ---- set defaults:
        kwargs = set_defaults(kwargs,ad.gtltcube_defaults)

        # ---- Set parameters:
        params        = ['evfile','outfile','scfile']
        gtltcube        = set_params(gt.expCube,infile,outfile,config,params, **kwargs)
        if gtltcube == -2:
            return gtltcube

        return run_gtobject(gtltcube)

    def run_psf(self,ltcube,outfile,config, **kwargs):
        """
        Run gtpsf to calculate the exposure averaged psf

        Parameters
        ----------
        ltcube:        string, full path to input file (output of gtltcube)
        outfile:        string, full path to output file
        config:        dictionary with main analysis parameters:
                    'expcube','outfile','irfs','evtype','ra','dec','emin','emax'
        kwargs
        ------
        Correspond to parameters which are usually not changed. See Notes and specFeat.base.AnalysisDefaults. 

        Returns
        -------
        boolean with error codes:
            0: everything went ok
            -1: Runtime error
            -2: Input file not found

        Notes
        -----
        For extensive help see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtpsf.txt
        """
        # ---- set defaults:
        kwargs = set_defaults(kwargs,ad.gtpsf_defaults)

        # ---- Set parameters:
        params        = ['expcube','outfile','irfs','evtype','ra','dec','emin','emax']
        gtpsf         = set_params(gt.GtApp('gtpsf'),ltcube,outfile,config,params, **kwargs)
        if gtpsf == -2:
            return gtpsf

        return run_gtobject(gtpsf)

    def run_expcube2(self,infile,outfile,config, **kwargs):
        """
        Run expcube.

        Parameters
        ----------
        infile:        string, full path to input file (output of gtltcube)
        outfile:        string, full path to output file
        config:        dictionary with main analysis parameters:
                    'emin','emax','ebins','exprad','irfs'
        kwargs
        ------
        Correspond to parameters which are usually not changed. See Notes and specFeat.base.AnalysisDefaults. 

        Returns
        -------
        boolean with error codes:
            0: everything went ok
            -1: Runtime error
            -2: Input file not found

        Notes
        -----
        For extensive help see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtexpcube2.txt
        """
        # ---- set defaults:
        kwargs = set_defaults(kwargs,ad.gtexpcube2_defaults)

        # ---- Set parameters:
        params        = ['infile','outfile','emin','emax','irfs','cmap']
        gtexpcube2        = set_params(gt.gtexpcube2,infile,outfile,config,params, **kwargs)
        if gtexpcube2 == -2:
            return gtexpcube2
        gtexpcube2        = set_spatial_e_binning(gtexpcube2,config,kwargs)

        return run_gtobject(gtexpcube2)

    def run_srcmaps(self,ltcube,cmap,expcube2,srcmdl,outfile,config, tries = 5, **kwargs):
        """
        Run gtsrcmaps.

        Parameters
        ----------
        ltcube:        string, full path to ltcube file (output of gtltcube)
        cmap:        string, full path to ccube file (output of gtbin)
        expcube2:        string, full path to expcube2 file (output of gtexpcube2)
        srcmdl:        string, full path to xml source model file
        outfile:        string, full path to output file
        irf:        string, instrumental response function
        config:        dictionary with main analysis parameters:
                    'scfile','irfs'
        kwargs
        ------
        tries:        if hdu list of srcmap is less than 3, try recomputing srcmap [tries] number of times
        Remaining kwargs correspond to parameters which are usually not changed. See Notes and specFeat.base.AnalysisDefaults. 


        Returns
        -------
        boolean with error codes:
            0: everything went ok
            -1: Runtime error
            -2: Input file not found
            -3: If hdu number in srcmap is less than 3

        Notes
        -----
        For extensive help see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtexpcube2.txt
        """
        # ---- set defaults:
        kwargs = set_defaults(kwargs,ad.gtsrcmaps_defaults)
        for k in kwargs.keys():
            gt.srcMaps[k] = kwargs[k]

        # ---- Set parameters:
        gt.srcMaps['scfile']        = config['scfile']
        gt.srcMaps['outfile']        = outfile
        gt.srcMaps['irfs']        = config['irfs']
        gt.srcMaps['expcube']        = ltcube
        gt.srcMaps['cmap']        = cmap
        gt.srcMaps['bexpmap']        = expcube2
        gt.srcMaps['srcmdl']        = srcmdl

        for f in ['expcube','cmap','bexpmap','srcmdl']:
            if not exists(gt.srcMaps[f]):
                logging.error('*** {0:s} not found in {1:s}, returning -2'.format(f,gt.srcMaps[f]))
                return -2
        ver,nTry = -3,0
        while(ver < 0 and nTry < tries):
            stat        = run_gtobject(gt.srcMaps)
            ver                = verify_srcmap(outfile)
            stat        = stat if ver == 0 else ver
            if ver < 0:
                call(['rm',outfile])
            nTry        += 1

        return stat

    def run_ltsum(self,infile1,outfile, infile2 = None, tmpdir = None, **kwargs):
        """
        Run expcube.

        Parameters
        ----------
        infile1:        string, full path to input file(s) (output of gtltcube) - if wildcard and more than one file
                        is found, these ltcubes are summed
        outfile:        string, full path to output file

        kwargs
        ------
        infile2:        second ltcube. Only required if infile1 is only one ltcube file.
        Correspond to parameters which are usually not changed. See Notes and specFeat.base.AnalysisDefaults. 
        tmpdir:                str, temporary directory for all_cubes.txt, if not given, dirname of infile1 is used

        Returns
        -------
        boolean with error codes:
            0: everything went ok
            -1: Runtime error
            -2: Input file not found

        Notes
        -----
        For extensive help see http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtltsum.txt
        """
        # ---- set defaults:
        kwargs = set_defaults(kwargs,ad.gtexpcube2_defaults)

        files = glob(infile1)
        if len(files) == 1:
            gt.addCubes['infile1']        = infile1
            if infiles2 == None:
                logging.error('*** infile2 cannot be one if infile1 is only one ltcube, returning -2')
                return -2
            else:
                gt.addCubes['infile2']        = infile2
        elif len(files) > 1:
            if tmpdir == None:
                listname        = join(dirname(infile1),'all_cubes.txt')
            else: 
                listname        = join(tmpdir,'all_cubes.txt')
            f                = open(listname,'w')
            for fi in files:
                f.write('{0:s}\n'.format(fi))
            f.close()
            gt.addCubes['infile1']        = '@' + listname
        else:
            logging.error('*** infile1 {0:s} not found, returning -2'.format(infile1))
            return -2

        gt.addCubes['outfile']        = outfile

        return run_gtobject(gt.addCubes)
# ========================================================================== #
# === the script =========================================================== #
# ========================================================================== #

if __name__ == "__main__":
    parser=OptionParser()
    parser.add_option("-c","--config",dest="c",help="config yaml file to use",action="store")
    parser.add_option("-r","--run_analysis",dest="r",help="specify what should be prepared",action="store", 
            choices = ['ft1','ltcube','obsExpFull','obsExpFullMC','srcMapFull',\
            'srcMapFullMC','SED','FullSED','FullSEDMC','test','select','psf','ft1lc','rspgen','rspgen_lc'])
    (opt, args) = parser.parse_args()

    # set the logger
    init_logging('INFO')

    # load config
    config = yaml.load(open(opt.c))

    # check if we are on the lsf cluster, if so, tmpdir on cluster scratch is made, 
    # pfiles are set to it, and job_id is determined and returned
    tmpdir, job_id = init_lsf()
    chdir(tmpdir)    # go to tmp directory
    logging.info('Entering directory {0:s}'.format(tmpdir))

    # start analysis
    logging.info('starting {0:s} preparation ...'.format(opt.r))
    pps = PreparePointSource()

    # some test
    if opt.r == 'test':        
        logging.info('files present: {0}'.format(glob('../*')))
        logging.info('sys path: {0:s}'.format(sys.path))
        logging.info('job id: {0:n}'.format(job_id))
        sys.exit(0)

    # extract FT1 file from Jean Ballet's files
    if opt.r == 'ft1':        
        outfile                                = {'ft1': join(tmpdir,basename(config['ft1']))}
        ops                                = ['select']
        localFilesKeys                        = ['infile']

    # run gtselect
    if opt.r == 'select':        
        outfile                                = {'select': join(tmpdir,basename(config['select']))}
        ops                                = ['select']
        localFilesKeys                        = ['infile']

    # LT cube calculation
    if opt.r == 'ltcube':        
        config['tmin'],config['tmax']        = set_tbins(job_id,config)
        outfile                                = {'ltcube': join(tmpdir,'ltcube{0:04n}.fits'.format(job_id))}
        ops                                = ['select','mktime','ltcube']
        localFilesKeys                        = ['infile','scfile']

    if opt.r == 'ft1lc':        
        outfile                                = {'ft1lc': join(tmpdir,'ft1lc.fits')}
        ops                                = ['select','mktime','lc']
        localFilesKeys                        = ['infile','scfile']

    if opt.r.find('rspgen') >= 0:        
        if opt.r == ('rspgen_lc'):
            config['tmin'],config['tmax']        = set_tbins(job_id,config)
            outfile                                = {
                                                    'rsp': join(tmpdir,'rsp{0:04}.fits').format(job_id),
                                                    'pha1': join(tmpdir,'pha1{0:04}.fits').format(job_id)
                                                    }
        else:
            outfile                                = {
                                                    'rsp': join(tmpdir,'rsp.fits'),
                                                    'pha1': join(tmpdir,'pha1.fits')
                                                    }

        ops                                = ['select','mktime','pha1','rspgen']
        localFilesKeys                        = ['infile','scfile']

    # PSF calculation + number of counts
    if opt.r == 'psf':        
        config['tmin'],config['tmax']        = set_tbins(job_id,config)
        outfile                                = {
                                            'ltcube': join(tmpdir,'ltcube{0:04n}.fits'.format(job_id)),
                                            'psf': join(tmpdir,'psf{0:04n}.fits'.format(job_id))
                                        }
        ops                                = ['select','mktime','ltcube','psf']
        localFilesKeys                        = ['infile','scfile']


    # SED calculation
    if opt.r == 'SED':        
        config['emin'],config['emax']        = set_ebins(job_id,config)
        config['ebins']                        = 1
        outfile        = {        'ccube' :        join(tmpdir,'ccube{0:04n}.fits'.format(job_id)),
                        'expcube':        join(tmpdir,'expcube{0:04n}.fits'.format(job_id)),
                        'cmap':                join(tmpdir,'cmap{0:04n}.fits'.format(job_id))
                    }
        ops                                = ['select','mktime','ccube','expcube2','cmap']
        localFilesKeys                        = ['infile','scfile','ltcube','xmlmodel']
        for k in ['cmap', 'ccube','expcube']:
            config[k]        = outfile[k]
    # Full SED calculation
    elif opt.r.find('FullSED') >= 0:        
        ops                                = ['select','mktime','ccube','expcube2','cmap']
        localFilesKeys                        = ['infile','scfile','ltcube','xmlmodel']
        emin,emax,ebins                        = config['emin'],config['emax'],config['ebins']        # save the original values
        if opt.r.find('MC') >= 0:
            config['infile']                = config['infile'].replace('*','{0:05n}'.format(job_id))
            config['savedir']                = mkdir(join(config['savedir'],'{0:05n}/SED/'.format(job_id)))

    # CCube and ExpCube calculation for full E and T range
    elif opt.r == 'obsExpFull':        
        outfile        = {        'ccube' :        join(tmpdir,'ccube.fits'),
                        'expcube':        join(tmpdir,'expcube.fits'),
                        'mktime':        join(tmpdir,'mktime.fits'),
                        'cmap':                join(tmpdir,'cmap.fits')
                    }
        ops                = ['select','mktime','ccube','cmap','expcube2']
        localFilesKeys        = ['infile','scfile','ltcube']
        config['cmap']        = outfile['cmap']

    # Source Map calculation for full E and T range
    elif opt.r == 'srcMapFull':        
        outfile                = {'srcmap' : join(tmpdir,'srcmap.fits')}
        ops                = ['srcmap']
        localFilesKeys        = ['scfile','ltcube','ccube','expcube','xmlmodel']

    # CCube and ExpCube calculation for full E and T range for MC data
    elif opt.r == 'obsExpFullMC':        
        config,localFilesKeys        = determine_infilesMC(config, tmpdir, job_id)        
        config['savedir']        = mkdir(join(config['savedir'],'{0:05n}'.format(job_id)))
        outfile        = {        'ccube' :        join(tmpdir,'ccube.fits'),
                        'expcube':        join(tmpdir,'expcube.fits'),
                        'mktime':        join(tmpdir,'mktime.fits'),
                        'cmap':                join(tmpdir,'cmap.fits')
                    }
        ops                = ['select','mktime','ccube','expcube2','cmap']
        config['cmap']        = outfile['cmap']


    # Source Map calculation for full E and T range for MC data
    elif opt.r == 'srcMapFullMC':        
        config['savedir']        = mkdir(join(config['savedir'],'{0:05n}'.format(job_id)))
        config['ccube']                = join(config['savedir'],'ccube.fits')
        config['expcube']        = join(config['savedir'],'expcube.fits')
        if not config['xmlkey'] == '':
            config['savedir']        = mkdir(join(config['savedir'],config['xmlkey']))

        outfile                = {'srcmap' : join(tmpdir,'srcmap.fits')}
        ops                = ['srcmap']
        localFilesKeys        = ['scfile','ltcube','ccube','expcube','xmlmodel']

    # if on cluster, copy infiles to cluster scratch
    if job_id:
        for k in localFilesKeys:
            if config[k].find('@') >= 0:
                config[k] = '@' + copy2scratch(config[k].lstrip('@'),tmpdir)
            else:
                if opt.r.find('FullSED') >= 0 and k == 'infile':        
                    copy2scratch(config[k],join(tmpdir,'ft1.fits'))        # different name for original ft1 file in SED production required
                    config[k] = join(tmpdir,'ft1.fits')
                else:
                    config[k] = copy2scratch(config[k],tmpdir)


    nsteps = config['ebins'] if opt.r.find('FullSED') >= 0 else 1
    for i in range(nsteps):
        if nsteps > 1:
            config['emin'],config['emax'],config['ebins']        = emin,emax,ebins        # restore original values
            config['emin'],config['emax']        = set_ebins(i+ 1,config)
            config['ebins']                        = 1
            outfile        = { 'ccube' :        join(tmpdir,'ccube{0:04n}.fits'.format(i + 1)),
                            'expcube':        join(tmpdir,'expcube{0:04n}.fits'.format(i + 1)),
                            'cmap':        join(tmpdir,'cmap{0:04n}.fits'.format(i + 1))
                    }
            for k in ['cmap','ccube','expcube']:
                config[k]        = outfile[k]
    # --- GT Select ------------ #
        status = {}
        outSelect = outfile['ft1'] if opt.r == 'ft1' else join(tmpdir,'select.fits')
        outSelect = outfile['select'] if opt.r == 'select' else outSelect

        try:
            ops.index('select')

            ### For pulsar analysis: loop over phasemin / phasemax lists if given ###
            ### Does not work because of this error in gtselect
            ### Caught St13runtime_error at the top level: DSS keywords in /scratch/mmeyer.991784/991784.2MrxLH/select00001.fits do not match those 
            ### in /scratch/mmeyer.991784/991784.2MrxLH/select00000.fits
#            if type(config['phasemin']) == list:
#                if not type(config['phasemax']) == list or not len(config['phasemin']) == len(config['phasemax']):
#                    logging.error('Incompatible phasemin and phasemax keywords: {0[phasemin]}, {0[phasemax]}'.format(config))
#                    sys.exit(101)
#
#                phasemax = config['phasemax']
#                phasemin = config['phasemin']

#                out = []
#                for i,p in enumerate(phasemin):
#                    config['phasemax'] = phasemax[i]
#                    config['phasemin'] = p
#                    out.append(join(tmpdir,'select{0:05n}.fits'.format(i)))
#                    pps.run_select(config['infile'],out[-1],config)

#                config['phasemin'] = 0.
#                config['phasemax'] = 1.

#                config['infile'] = join(tmpdir,'combined.txt')
#                savetxt(config['infile'],np.array(out),fmt = '%s')
#                config['infile'] = '@' + config['infile']
                
            logging.info('using infile {0:s}'.format(config['infile']))
            status['select'] = pps.run_select(config['infile'],outSelect,config)

            if not check4events(outSelect) and not status['select']:
                logging.warning('*** No events in select.fits -> setting status to -3! Check ROI cuts!')
        except ValueError:
            status['select'] = 0

    # --- GT Maketime ---------- #
        try:
            ops.index('mktime')
            outMktime = join(tmpdir,'mktime.fits')
            if not status['select']:
                status['mktime'] = pps.run_mktime(
                        outSelect,
                        outMktime,
                        config
                        )
                # check if there are any events in the mktime.fits output
                if not check4events(outMktime) and not status['mktime']:
                    logging.warning('*** No events in mktime.fits! Check ROI cuts!')
            else:
                sys.exit(status['select'])
        except ValueError:
            status['mktime'] = 0


        # --- gt bin ccube --------- #
        try:
            ops.index('ccube')
            if not status['mktime']:
                status['ccube'] = pps.run_bin(outMktime,outfile['ccube'],config, algorithm = 'CCUBE')
            else:
                sys.exit(status['mktime'])
        except ValueError:
            status['ccube'] = 0

        # --- gt bin ccube --------- #
        try:
            ops.index('cmap')
            if not status['mktime']:
                status['cmap'] = pps.run_bin(outMktime,outfile['cmap'],config, algorithm = 'CMAP')
            else:
                sys.exit(status['mktime'])
        except ValueError:
            status['cmap'] = 0

        # --- gt bin lc ------------ #
        try:
            ops.index('lc')
            if not status['mktime']:
                status['lc'] = pps.run_bin(outMktime,outfile['ft1lc'],config, algorithm = 'LC')
            else:
                sys.exit(status['mktime'])
        except ValueError:
            status['lc'] = 0

        # --- gt bin pha1 ---------- #
        try:
            ops.index('pha1')
            if not status['mktime']:
                status['pha1'] = pps.run_bin(outMktime,outfile['pha1'],config, algorithm = 'PHA1')
            else:
                sys.exit(status['mktime'])
        except ValueError:
            status['pha1'] = 0

        # --- gt rsp gen ---------- #
        try:
            ops.index('rspgen')
            if not status['pha1']:
                status['rspgen'] = pps.run_rspgen(outfile['pha1'],outfile['rsp'], config)
            else:
                sys.exit(status['pha1'])
        except ValueError:
            status['rspgen'] = 0

        # --- GT LTcube ------------ #
        try:
            ops.index('ltcube')
            if not status['mktime']:
                status['ltcube'] = pps.run_ltcube(outMktime,outfile['ltcube'],config)
            else:
                sys.exit(status['mktime'])
        except ValueError:
            status['ltcube'] = 0

        # --- PSF ------------------ #
        try:
            ops.index('psf')
            if not status['ltcube']:
                logging.info('starting psf calc')
                status['psf'] = pps.run_psf(outfile['ltcube'],outfile['psf'],config)
                outfile.pop('ltcube',None) #remove the ltcube
            else:
                sys.exit(status['mktime'])
        except ValueError as e:
            logging.info(e)
            status['ltcube'] = 0

        # --- GT ExpCube2 ---------- #
        try:
            ops.index('expcube2')
            if not status['ltcube']:
                status['expcube'] = pps.run_expcube2(config['ltcube'],outfile['expcube'],config)
            else:
                sys.exit(status['ltcube'])
        except ValueError:
            status['expcube'] = 0

        # --- GT SrcMap ------------ #
        try:
            ops.index('srcmap')
            if not (status['ltcube'] or status['ccube'] or status['expcube']):
                status['srcmap'] = pps.run_srcmaps(config['ltcube'],config['ccube'],config['expcube'],config['xmlmodel'],outfile['srcmap'],config)
        except ValueError:
            status['srcmap'] = 0

        # -------------------------- #
        # copy the final outfile
        copy2scratch(outfile,config['savedir'])
        
    sys.exit(status)
