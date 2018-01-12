#!/usr/bin/env python

# ---- Imports --------------------------------- #
from optparse import OptionParser
from specFeat.base.lightcurve import Lightcurve
from specFeat.base.defaults import parser_defaults
from os.path import *
from sys import exit
from specFeat.utils.psf import PSF
from glob import glob
import logging
import numpy as np
# ---------------------------------------------- #

"""submit script to cluster to make ft1 lc file"""

if __name__ == "__main__":
    parser=OptionParser()
    parser = parser_defaults(parser)

    parser.add_option(
	"-k","--key",dest="k",
	    help="key for analysis, either edisp or SED, use edisp for finer binning and larger energy window (default: edisp)",
	    action="store", default = 'edisp',choices = ['edisp','SED']
	    )
    parser.add_option(
	"-T","--dT",dest="T",
	    help="Time interval. If 0., use full intervall", type = 'float', default = 0.
	    )
    parser.add_option("-S","--Step",dest="S",help="analysis step",
	action="store", default = 'ft1lc', choices = ['ft1lc','psf','rspgen'])
    parser.add_option("-F","--file",dest="F",help="fits file with GTI extensions",action="store")
    parser.add_option("-n","--nGTI",dest="n",help="number of GTI you want to use",action="store", type = 'int', default = 0)
    parser.add_option("-P","--PSFcut",dest="P",
	help="if -1, change the radius of the analysis to the 68% PSF containment. If > 0. change to this radius in deg.",
	action="store", type = 'float', default = 0)

    (opt, args) = parser.parse_args()

    try:
	if exists(opt.c):
	    lc = Lightcurve(opt.c, color = True, dt = opt.T if opt.T > 0. else None)
	else:
	    raise IOError('File {0} not found or not ready for opening.'.format(opt.c))
	    exit(-2)
    except TypeError as e:
	print "You have to provide a config yaml file with the -c option\n{0}".format(e.args)
	exit(-1)

    if not opt.F == None:
	lc.set_gti(opt.F, opt.n, dt = opt.T if opt.T > 0. else None)

# make a cut on the psf radius
    if opt.P:
	logging.info('Cutting on PSF containment radius ...')

	if opt.P == -1.:
	    r68irf = []
	    for i,irf in enumerate(lc.config['irfs']):
		pfiles = glob(join(dirname(lc.outfiles['edisp'][irf]['psf']),'psf*.fits'))
		if len(pfiles) > 1:
		    pfiles = sorted(pfiles, key = lambda f: int(basename(f).split('.')[0][-4:]))

	    # get the 68 containment radius and exposure for each file
		r68,avgexp = [],[]
		for j,f in enumerate(pfiles):
		    p = PSF(f)
		    p.get_containment(0.68, slice = slice(1))
			    
		    avgexp.append(p.get_avg_exposure())
		    r68.append( p.r[0][0] )
		    del p

		r68,avgexp = np.array(r68),np.array(avgexp)
		r68irf.append(r68.mean())
		logging.info('For Irf {0:s}: 68% containment radius is {1:.2f}'.format(irf,r68irf[i]))

	lc.config['rad'] = opt.P if opt.P > 0. else r68irf[0]
	logging.info('Radius is now {0[rad]:.2f}'.format(lc.config))
	lc.__dict__.update()
	logging.info('Setting the RoI size to {0[rad]:.2f}'.format(lc.config))

# generate a LC from FT1 file
    if opt.S == 'ft1lc':
	lc.makeFT1LC(clobber = opt.C,
		    queue = opt.q,
		    time = opt.W,
		    dry = opt.d,
		    concurrent	= opt.Co,
		    key = opt.k,
		    dt = opt.T if opt.T > 0. else None
	)
    elif opt.S == 'psf':
	lc.makePSF(clobber = opt.C,
		    queue = opt.q,
		    time = opt.W,
		    dry = opt.d,
		    concurrent	= opt.Co,
		    key = opt.k,
		    dt = opt.T if opt.T > 0. else None
	)
    elif opt.S == 'rspgen':
	lc.makeRsp(clobber = opt.C,
		    queue = opt.q,
		    time = opt.W,
		    dry = opt.d,
		    concurrent	= opt.Co,
		    key = opt.k,
		    dt = opt.T if opt.T > 0. else None
	)

    exit(0)
