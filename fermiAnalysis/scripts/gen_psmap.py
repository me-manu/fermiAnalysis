import argparse
import os
import logging
from fermipy.gtanalysis import GTAnalysis
from fermiAnalysis import setup
from fermiAnalysis.plotting import tsmap_plot_nice
from fermiAnalysis.tsps import generate_psmap, PSFReader

if __name__ == "__main__":
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('--overwrite', 
                        help='Overwrite existing files', action="store_true")
    parser.add_argument('--make-plots', 
                        help='generate psmap and ts map plots', action="store_true")
    parser.add_argument('-s', '--state', default = ['avgspec'],
                        choices = ['setup','avgspec','avgspec_ebl','lcbin', 'lcmonthly'],
                        help='Analysis state')
    parser.add_argument('-i', required=False, default = 0, 
                    help='Set local or scratch calculation', type=int)
    parser.add_argument('--prob_epsilon', default = 1e-7, 
                    help='epsilon parameter for gtpsmap', type=float)

    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i=args.i, logging_level="INFO")
    logging.info('Reloading roi')
    gta.load_roi(args.state)

    psmaps = generate_psmap(gta, args.state, prob_epsilon=args.prob_epsilon, overwrite=args.overwrite, make_plots=args.make_plots)

    if args.make_plots:
        # make ps map plots
        psf_file = os.path.join(gta.config['fileio']['workdir'], "psf_01.fits")
        psf_file2 = os.path.join(gta.config['fileio']['workdir'], "psf_00.fits")
        psf = PSFReader(psf_file)
        psf2 = PSFReader(psf_file2)
        tsmap_plot_nice(gta, args.state, psf, psf2=psf2)
