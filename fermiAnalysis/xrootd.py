# --- Imports ------------------------------------------- #
import logging
import sys
import yaml
from os import path
import os 
from subprocess import call
import argparse
from glob import glob
from fermiAnalysis.tools import *
from haloanalysis.batchfarm import lsf,utils
import shlex
import copy
from time import sleep
from fermipy.gtanalysis import GTAnalysis
import fermiAnalysis
# ------------------------------------------------------- #

# === the class ============================================================ #
class XrootdSelect(object):
    """ 
    Class that contains functions for gtselect using xrootd
    of the Fermi point source analysis
    """
    fname = __module__.replace('.','/') + '.py'
    #xrdft1list = path.join(path.dirname(fermiAnalysis.__file__), 'xrootd-ft1.txt')
    xrdft1list = path.join(path.dirname(fermiAnalysis.__file__), 'xrootd-ft1-new.txt')

    def __init__(self, ft1 = True):
        if ft1:
            self.flist = XrootdSelect.xrdft1list
        return

    def create_shortlist(self, tmin, tmax, outdir):
        """
        Create a short file list that contains 
        only the files with events within tmin and tmax

        Parameters
        ----------
        tmin: float
            minimum observation time in MET

        tmax: float
            maximum observation time in MET

        outdir: str
            path to output file

        Returns
        -------
        path to new output file
        """
        with open(self.flist) as f: 
            li = f.readlines()
        files = []

        for l in li:
            tstart = int(path.basename(l).split('_')[-2].strip('r'))
            if tstart >= tmin and tstart < tmax:
                files.append(l)

        if not len(files):
            raise Exception("No files selected!")
        else:
            logging.info("Shortlist contains {0:n} files".format(len(files)))

        with open(path.join(outdir,'tmpfilelist.txt'), 'w') as f:
            for fi in files:
                f.write(fi)

        return path.join(outdir,'tmpfilelist.txt')


    def make_cmd_str(self,**conf):
        """
        Make a string to run gtselect with chosen parameters

        kwargs        
        ------
        conf:        dict with parameters

        Returns
        -------
        str with options for gtselect
        """
        conf.setdefault("evclass", 128)
        conf.setdefault("ra", "INDEF")
        conf.setdefault("dec","INDEF")
        conf.setdefault("rad", 180.)
        conf.setdefault("emin", 100.)
        conf.setdefault("emax", 300000.)
        conf.setdefault("tmin", "INDEF")
        conf.setdefault("tmax", "INDEF")
        conf.setdefault("zmax", 180.)
        conf.setdefault("evtype", 3)
        conf.setdefault("phasemin", 0.)
        conf.setdefault("phasemax", 1.)
        params = ['evclass','ra','dec','rad','emin','emax','tmin','tmax','infile','outfile','zmax','evtype','phasemin','phasemax']

        cmdStr = ''
        for k in params:
            if k == 'phasemin' and conf[k] == None or conf[k] == 0.:
                #conf[k] = 0.
                continue
            elif k == 'phasemax' and conf[k] == None or conf[k] == 1.:
                #conf[k] = 1.
                continue
            cmdStr +='{0:s}={1} '.format(k,conf[k])
        return cmdStr
