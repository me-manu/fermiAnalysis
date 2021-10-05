from fermipy.utils import init_matplotlib_backend
import argparse
import os
from fermiAnalysis import setup
from fermiAnalysis.ext_funcs import fit_igmf_halo_scan

def main():
    usage = "usage: %(prog)s -c config.yaml"
    description = "Run the analysis"
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('-c', '--conf', required=True)
    parser.add_argument('--halo-template-dir', required=True, help='Directory with halo template fits files')
    parser.add_argument('--halo-template-suffix', required=True, help='suffix for halo template fits fits')
    parser.add_argument('--file-suffix', help='additional suffix for output files', default='')
    parser.add_argument('--state', default = ['avgspec'],
                        choices = ['avgspec','avgspec_ebl'],
                        help='Analysis state')
    parser.add_argument('-i', required=False, default = 0, 
                        help='Set local or scratch calculation', type=int)
    parser.add_argument('--overwrite', action="store_true", help='overwrite existing single files')
    parser.add_argument('--generate-seds', action="store_true", help='Generate SEDs during analysis')
    parser.add_argument('--generate-maps', action="store_true", help='Generate TS and residual maps during analysis')

    args = parser.parse_args()

    gta, config, fit_config, job_id  = setup.init_gta(args.conf, i = args.i, logging_level = "INFO")
    gta.logger.info('Running fermipy setup')

    init_matplotlib_backend()

    gta.logger.info('reloading {0:s}'.format(args.state))
    gta.load_roi(args.state) # reload the average spectral fit

    modelname = "{0:s}_{1:s}{2:s}".format(args.state,
                                     '_'.join([k for k in args.halo_template_dir.split('/')[-4:] if not 'spec' in k]), args.file_suffix)

    gta.logger.info("Using modelname: {0:s}".format(modelname))
    # change the outdir
    # not the greatest that I'm not using the API here, 
    # but no other possibility
    gta._outdir = os.path.join(gta.outdir, 'igmf_' + modelname + '/')
    if not os.path.exists(gta.outdir):
        os.makedirs(gta.outdir)
    gta.logger.info("Set new outdir: {0:s}".format(gta.outdir))

    gta.logger.info("reloaded ROI had log likelihood value: {0:.2f}".format(-gta.like()))
    halo_profile_tied = fit_igmf_halo_scan(gta, modelname,
                                           config['selection']['target'],
                                           args.halo_template_dir,
                                           model_idx=job_id,
                                           halo_template_suffix=args.halo_template_suffix,
                                           injection_spectrum='PLSuperExpCutoff',
                                           injection_par2_name='Cutoff',
                                           injection_norm_name='Prefactor',
                                           injection_scale_name='Scale',
                                           index_par_name='Index',
                                           free_bkgs=True,
                                           generate_maps=args.generate_maps,
                                           generate_seds=args.generate_seds,
                                           distance_free_norm=3.,  # at 1e-14 G, above 2. deg about 10% of cascade photons are beyond 2 deg at 1GeV
                                           z=fit_config['z'], 
                                           ebl_model_name='dominguez',
                                           optimizer='MINUIT')

    return gta, halo_profile_tied

if __name__ == '__main__':
    res = main()
