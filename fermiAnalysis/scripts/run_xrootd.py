from os import path
import logging
import argparse
import shlex
import fermiAnalysis as fa
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fermiAnalysis import setup
from fermiAnalysis import xrootd
from fermiAnalysis.utils import set_free_pars_avg,fit_with_retries
from fermiAnalysis.prepare import PreparePointSource 
from fermiAnalysis.batchfarm import utils,lsf

if __name__ == '__main__':
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('-i', required=False, default = 0, type = int)
    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i = args.i, logging_level = "INFO")

    # start analysis
    xrd = xrootd.XrootdSelect()
    gta.config['selection']['infile'] = '@' + xrd.create_shortlist(gta.config['selection']['tmin'], 
                                                                    gta.config['selection']['tmax'],
                                                                    gta.workdir)
    gta.config['selection']['rad'] = gta.config['binning']['roiwidth']
    gta.config['selection']['outfile'] = gta.config['data']['evfile']

    command = xrd.make_cmd_str(**gta.config['selection'])
    command = """/afs/slac/g/glast/applications/xrootd/PROD/bin/xrdprel -g gtselect {0:s}""".format(command)

    logging.info('Executing command: {0:s}'.format(command))
    utils.sleep(1.)

    call(shlex.split(command))
