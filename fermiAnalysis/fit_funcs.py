import os
import glob
import sys
import copy
import itertools
import logging

import numpy as np
from .utils import stack_files
from astropy.table import Column
from fermipy.utils import get_parameter_limits

    
def fit_region(gta,modelname,src_name,loge_bounds=None, **kwargs):

    skip_opt = kwargs.get('skip_opt',[])
    
    gta.logger.info('Starting Region Fit %s'%(modelname))
    lnl0 = -gta.like()    
    gta.logger.info('%s Model Likelihood: %f'%(modelname,lnl0))
    gta.print_params()
    
    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0],loge_bounds[1])
    
    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model_pl20 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model_pl27 = { 'SpatialModel' : 'PointSource', 'Index' : 2.7 }
    model3 = { 'SpatialModel' : 'Gaussian', 'Index' : 2.0, 'SpatialWidth' : 0.1 }
    model4 = { 'SpatialModel' : 'RadialDisk', 'Index' : 2.0,
               'SpatialWidth' : 0.1 * 0.8246211251235321 }
    
    gta.optimize(skip=skip_opt, shape_ts_threshold=9.0)

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    skydir = gta.roi[src_name].skydir
    
    gta.free_sources(False)
    gta.free_sources(skydir=skydir,distance=1.5, pars='norm')
    gta.free_sources(skydir=skydir,distance=1.0, pars='shape', exclude=diff_sources)
    gta.free_source(src_name)
    gta.fit()
    gta.update_source(src_name,reoptimize=True)
    gta.write_roi(modelname + '_roi', make_plots=True)

    gta.print_roi()
    gta.print_params()

    lnl1 = -gta.like()
    gta.logger.info('%s Model Likelihood: %f'%(modelname,lnl1))
    gta.logger.info('%s Model Likelihood Delta: %f'%(modelname,lnl1-lnl0))
    
    # TS Maps
    maps_model_pl20 = gta.tsmap(modelname, model=model_pl20,
                            loge_bounds=loge_bounds, make_plots=True)
    gta.tsmap(modelname, model=model_pl27,
              loge_bounds=loge_bounds, make_plots=True)
    maps_model_pl20_nosource = gta.tsmap('%s_nosource'%modelname,
                                     model=model_pl20, exclude=[src_name],
                                     loge_bounds=loge_bounds, make_plots=True)
    maps_model_pl27_nosource = gta.tsmap('%s_nosource'%modelname,
                                     model=model_pl27, exclude=[src_name],
                                     loge_bounds=loge_bounds, make_plots=True)
    #maps_model4_nosource = gta.tsmap('%s_nosource'%modelname,
    #                                 model=model4, exclude=[src_name],
    #                                 loge_bounds=loge_bounds, make_plots=True)    
    gta.residmap(modelname, model=model3,
                 loge_bounds=loge_bounds, make_plots=True)
                                                                        
    # SED Analysis
    gta.sed(src_name, outfile=modelname + '_sed_fixed',
            prefix=modelname + '_fixed',
            make_plots=True)
    
    gta.sed(src_name, outfile=modelname + '_sed',
            prefix=modelname,
            free_radius=1.0, make_plots=True)    

    gta.sed(src_name,outfile=modelname + '_sed_bin4',
            prefix=modelname + '_bin4', loge_bins=gta.log_energies[::2],
            free_radius=1.0, make_plots=True)
    
    psf_syst_scale = np.array([0.05,0.05,0.2])
    psf_fnlo = ([3.0,4.0,5.5],list(-1.0*psf_syst_scale))
    psf_fnhi = ([3.0,4.0,5.5],list(1.0*psf_syst_scale))

    # -----------------------------------------------------------------
    # Gaussian Analysis
    # -----------------------------------------------------------------
    kw = dict(spatial_model='RadialGaussian',
              free_radius=1.0, make_tsmap=False)
    
    gta.extension(src_name, outfile=modelname + '_ext_gauss_ext',
                  prefix=modelname + '_gauss',
                  fit_position=True, free_background=True,
                  make_plots=True, update=True, **kw)

    gta.extension(src_name, outfile=modelname + '_ext_gauss_ext_psflo',
                  prefix=modelname + '_gauss_psflo',
                  psf_scale_fn=psf_fnlo, **kw)

    gta.extension(src_name, outfile=modelname + '_ext_gauss_ext_psfhi',
                  prefix=modelname + '_gauss_psfhi',
                  psf_scale_fn=psf_fnhi, **kw)

    gta.free_source(src_name)
    gta.fit()
    gta.update_source(src_name,reoptimize=True)
    gta.print_roi()
    gta.print_params()
    
    gta.sed(src_name,outfile=modelname + '_ext_gauss_sed',
            prefix=modelname + '_gauss',
            free_radius=1.0, make_plots=True)
    gta.sed(src_name,outfile=modelname + '_ext_gauss_sed_bin4',
            prefix=modelname + '_gauss_bin4', loge_bins=gta.log_energies[::2],
            free_radius=1.0, make_plots=True)
    gta.write_roi(modelname + '_ext_gauss_roi')

    gta.tsmap(modelname + '_ext_gauss', model=model_pl20,
              loge_bounds=loge_bounds, make_plots=True)
    gta.tsmap(modelname + '_ext_gauss', model=model_pl27,
              loge_bounds=loge_bounds, make_plots=True)

    # -----------------------------------------------------------------
    # Disk Analysis
    # -----------------------------------------------------------------
    gta.load_roi(modelname + '_roi')
    gta.reload_source(src_name)

    kw = dict(spatial_model='RadialDisk',
              free_radius=1.0, make_tsmap=False)
    
    gta.extension(src_name, outfile=modelname + '_ext_disk_ext',
                  prefix=modelname + '_disk',
                  fit_position=True, free_background=True,
                  make_plots=True, update=True, **kw)

    gta.extension(src_name, outfile=modelname + '_ext_disk_ext_psflo',
                  prefix=modelname + '_disk_psflo',
                  psf_scale_fn=psf_fnlo, **kw)

    gta.extension(src_name, outfile=modelname + '_ext_disk_ext_psfhi',
                  prefix=modelname + '_disk_psfhi',
                  psf_scale_fn=psf_fnhi, **kw)
 
    gta.free_source(src_name)
    gta.fit()
    gta.update_source(src_name,reoptimize=True)
    gta.print_roi()
    gta.print_params()
    
    gta.sed(src_name,outfile=modelname + '_ext_disk_sed',
            prefix=modelname + '_disk',
            free_radius=1.0, make_plots=True)
    gta.sed(src_name,outfile=modelname + '_ext_disk_sed_bin4',
            prefix=modelname + '_disk_bin4', loge_bins=gta.log_energies[::2],
            free_radius=1.0, make_plots=True)
    gta.write_roi(modelname + '_ext_disk_roi')
    
    gta.load_roi(modelname + '_roi')
    gta.reload_source(src_name)    
    gta.logger.info('Finished Region Fit %s'%(modelname))


def fit_halo_sed(gta,modelname,src_name,halo_width,
                 halo_index,spatial_model='RadialGaussian',
                 loge_bounds=None):

    gta.logger.info('Starting Halo SED Fit %s'%(modelname))
    
    halo_source_name = 'halo_' + spatial_model
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 1.0, 'max' : 4.5 }, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13, 'min' : 1E-5, 'max' : 1E4 },
        'SpatialModel' : spatial_model,
        'SpatialWidth' : 1.0
        }

    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

    gta.load_roi(modelname)
    if loge_bounds is not None:
        gta.set_energy_range(loge_bounds[0],loge_bounds[1])


    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
        
    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm', exclude=diff_sources)
    gta.write_xml(modelname + '_base')
    
    for i, w in enumerate(halo_width):

        halo_source_dict['SpatialWidth'] = w        
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name,halo_source_dict,free=True)

        # Do one fit with index free
        gta.set_parameter(halo_source_name,'Index',-2.0,
                          update_source=False)

        gta.fit()
        
        # SED w/ Index = 2.0
        gta.sed(halo_source_name,prefix='%s_%02i'%(modelname,i),
                fix_background=False, cov_scale=5.0)
        gta.write_roi('%s_halo_gauss_sed_%02i'%(modelname,i),
                      make_plots=False)

    gta.logger.info('Finished Halo SED Fit %s'%(modelname))

def fit_halo_scan(gta, modelname, src_name, halo_width,
                  halo_index, spatial_model='RadialGaussian',
                  loge_bounds=None, optimizer='NEWTON'):

    gta.logger.info('Starting Halo Scan %s'%(modelname))

    halo_source_name = 'halo_' + spatial_model
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 0.5, 'max' : 4.5 }, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13,
                        'min' : 1E-5, 'max' : 1E4 },
        'SpatialModel' : spatial_model,
        'SpatialWidth' : 1.0
        }

    outprefix = '%s_%s'%(modelname,halo_source_name)
    
    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

    #gta.load_roi(modelname)
    #if loge_bounds is not None:
    #    gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    skydir = gta.roi[src_name].skydir
    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(skydir=skydir,distance=1.0,pars='norm',
                     exclude=diff_sources)
    gta.write_xml(modelname + '_base')

    halo_tab = gta.roi.create_table([])
    halo_tab_idx_free = gta.roi.create_table([])
    halo_data = []
    halo_data_idx_free = []
        
    for i, w in enumerate(halo_width):

        gta.logger.info('Fitting Halo Width %.3f',w)
        
        halo_source_dict['SpatialWidth'] = w
        gta.load_xml(modelname + '_base')
        
        gta.add_source(halo_source_name, halo_source_dict, free=True)

        # Free Index
        gta.free_norm(halo_source_name)
        gta.fit(optimizer=optimizer)
        gta.sed(halo_source_name, prefix='%s_cov05_%02i'%(modelname,i),
                outfile='%s_cov05_%02i_sed'%(outprefix,i),
                free_radius=1.0, cov_scale=5.0,
                optimizer={'optimizer' : 'MINUIT'},
                make_plots=False)

        gta.sed(halo_source_name, prefix='%s_cov10_%02i'%(modelname,i),
                outfile='%s_cov10_%02i_sed'%(outprefix,i),
                free_radius=1.0, cov_scale=10.0,
                optimizer={'optimizer' : 'MINUIT'},
                make_plots=False)
        
        gta.free_parameter(halo_source_name,'Index')
        gta.fit(optimizer=optimizer)
        gta.free_parameter(halo_source_name,'Index',False)
        gta.update_source(halo_source_name,reoptimize=True,
                          optimizer={'optimizer' : optimizer})

        halo_data_idx_free += [copy.deepcopy(gta.roi[halo_source_name].data)]
        gta.roi[halo_source_name].add_to_table(halo_tab_idx_free)
        gta.write_roi('%s_%02i'%(outprefix,i),make_plots=False)

        gta.print_params(loglevel=logging.DEBUG)
        
        # Scan over fixed index
        for j, idx in enumerate(halo_index):

            gta.logger.info('Fitting Halo Index %.3f',idx)
            
            model_idx = i*len(halo_index) + j
            gta.set_norm(halo_source_name, 0.1, update_source=False)            
            gta.set_parameter(halo_source_name, 'Index', -1.0*idx,
                              update_source=False)
            
            gta.fit(update=False, optimizer=optimizer)

            gta.print_params(loglevel=logging.DEBUG)
            
            gta.update_source(halo_source_name,reoptimize=True,
                              optimizer={'optimizer' : optimizer})

            ul_flux = get_parameter_limits(gta.roi[halo_source_name]['flux_scan'],
                                           gta.roi[halo_source_name]['loglike_scan'])
            ul_eflux = get_parameter_limits(gta.roi[halo_source_name]['eflux_scan'],
                                            gta.roi[halo_source_name]['loglike_scan'])

            gta.roi[halo_source_name]['flux_err'] = ul_flux['err']
            gta.roi[halo_source_name]['eflux_err'] = ul_eflux['err']
            
            gta.logger.info('%s Halo Width: %6.3f Index: %6.2f TS: %6.2f Flux: %8.4g',
                            modelname,w,idx,
                            gta.roi[halo_source_name]['ts'],
                            gta.roi[halo_source_name]['flux'])
    
            #gta.write_roi('%s_%02i_%02i'%(outprefix,i,j),make_plots=False)
            halo_data += [copy.deepcopy(gta.roi[halo_source_name].data)]
            gta.roi[halo_source_name].add_to_table(halo_tab)
            
        gta.delete_source(halo_source_name,save_template=False) 

    np.save(os.path.join(gta.workdir,'%s_data.npy'%outprefix),halo_data)
    np.save(os.path.join(gta.workdir,'%s_data_idx_free.npy'%outprefix),
            halo_data_idx_free)

    tab_halo_width, tab_halo_index = np.meshgrid(halo_width,halo_index,indexing='ij')
    halo_tab['halo_width'] = np.ravel(tab_halo_width)
    halo_tab['halo_index'] = np.ravel(tab_halo_index)
    halo_tab_idx_free['halo_width'] = halo_width

    stack_files(sorted(glob.glob(os.path.join(gta.workdir,'%s*cov05*fits'%outprefix))),
                os.path.join(gta.workdir,'%s_cov05_sed.fits'%outprefix),
                new_cols=[Column(name='halo_width',data=halo_width, unit='deg')])

    stack_files(sorted(glob.glob(os.path.join(gta.workdir,'%s*cov10*fits'%outprefix))),
                os.path.join(gta.workdir,'%s_cov10_sed.fits'%outprefix),
                new_cols=[Column(name='halo_width',data=halo_width, unit='deg')])
    
    halo_tab.write(os.path.join(gta.workdir,'%s_data.fits'%outprefix),
                   overwrite=True)
    halo_tab_idx_free.write(os.path.join(gta.workdir,'%s_data_idx_free.fits'%outprefix),
                            overwrite=True)
    gta.logger.info('Finished Halo Scan %s'%(modelname))
    
    
def fit_halo(gta, modelname, src_name,
             spatial_model='RadialGaussian',
             loge_bounds=None, optimizer='NEWTON'):

    gta.logger.info('Starting Halo Fit %s'%(modelname))
    
    halo_source_name = 'halo_' + spatial_model
    halo_source_dict = {
        'SpectrumType' : 'PowerLaw', 
        'Index' : { 'value' : 2.0, 'scale' : -1.0, 'min' : 1.0, 'max' : 4.5 }, 
        'Scale' : 1000,
        'Prefactor' : { 'value' : 1E-5, 'scale' : 1e-13,
                        'min' : 1E-5, 'max' : 1E4 },
        'SpatialModel' : spatial_model,
        'SpatialWidth' : 1.0
        }

    outprefix = '%s_%s'%(modelname,halo_source_name)
    
    halo_source_dict['ra'] = gta.roi[src_name]['ra']
    halo_source_dict['dec'] = gta.roi[src_name]['dec']

#    gta.load_roi(modelname)
#    if loge_bounds is not None:
#        gta.set_energy_range(loge_bounds[0],loge_bounds[1])

    diff_sources = [s.name for s in gta.roi.sources if s.diffuse]
    
    gta.free_sources(False)
    gta.free_sources(distance=1.0,pars='norm',
                     exclude=diff_sources)

    # Find best-fit halo model
    halo_source_dict['SpatialWidth'] = 0.1
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_norm(halo_source_name)
    gta.extension(halo_source_name,update=True,
                  optimizer={'optimizer' : optimizer},
                  free_radius=1.0)

    # Fit spectrum
    gta.free_parameter(halo_source_name,'Index')
    gta.fit()

    # Re-fit extension
    gta.extension(halo_source_name,update=True,
                  optimizer={'optimizer' : optimizer},
                  free_radius=1.0)    

    # Re-fit Spectrum
    gta.fit()

    gta.update_source(halo_source_name,reoptimize=True,
                      optimizer={'optimizer' : optimizer})

    gta.print_params()
    
    gta.write_roi(outprefix,make_plots=False)    
    np.save(os.path.join(gta.workdir,'%s_data.npy'%outprefix),
            copy.deepcopy(gta.roi[halo_source_name].data))
    gta.delete_source(halo_source_name,save_template=False)
    
    gta.logger.info('Finished Halo Fit %s'%(modelname))
