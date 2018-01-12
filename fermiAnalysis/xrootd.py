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
    xrdft1list = path.join(path.dirname(fermiAnalysis.__file__), 'xrootd-ft1.txt')

    def __init__(self, ft1 = True):
        if ft1:
            self.flist = XrootdSelect.xrdft1list
        return

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
