import numpy as np
import logging
from os import path
import pyIrfLoader as irf_loader
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import fermiAnalysis
import time
try:  
    from ROOT import gROOT
    print "Loading macro"
    gROOT.LoadMacro(path.join(path.dirname(fermiAnalysis.__file__),
                                "TS_estimate_P8_C.so"))
except:
    try:
        gROOT.LoadMacro(path.join(path.dirname(fermiAnalysis.__file__),
                                "TS_estimate_P8.C") +"++")
    except:
        gROOT.LoadMacro(path.join(path.dirname(fermiAnalysis.__file__),
                                "TS_estimate_P8.C"))
try:
    from ROOT import get_evtlim
    from ROOT import initialize_TS, get_E1
    from ROOT import Glob_mult
    from ROOT import feed_sources
    from ROOT import Get_time_unc_mult, Get_flux_mult_src_minuit, Get_unc_mult
except ImportError as e:
    print "There was a problem loading the ROOT functions: {0}".format(e)
from array import array

def comp_E1(gta):
    src = gta.roi[gta.config['selection']['target']]
    emin = gta.config['selection']['emin']
    emax = gta.config['selection']['emin']
    index = src['dnde_index']
    flux = src['flux']
    irfs = gta.config['gtlike']['irfs']
    l = src['glon']
    b = src['glat']
    initialize_TS(path.dirname(fermiAnalysis.__file__) + '/')  
    logging.error("{0}".format([flux,index, l, b, emin, emax]))
    return get_E1(flux,index, l, b, emin, emax)

def comp_exposure_phi_new(gta, energy=None, nevt_max=None):
    """
    Compute the approximate exposure 
    trying to speed things up on 09/21/2020

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        The analysis object

    {options}

    energy: float or None
        The energy at which the exposure is evaluated.
        If none, central energy from selection is used
    """
    from tqdm import tqdm

    # load the effective areas
    irf_loader.Loader_go()
    irfs = gta.config['gtlike']['irfs']
    factory = irf_loader.IrfsFactory_instance()
    fronta = factory.create(irfs + "::FRONT")
    backa = factory.create(irfs + "::BACK")
    aeff_fa = fronta.aeff()
    aeff_ba = backa.aeff()

    if type(energy) == type(None):
        energy = float(np.sqrt(gta.config['selection']['emin'] * \
                            gta.config['selection']['emax']))

    # load the FT1 and FT2 files
    ft1fits = fits.open(gta.config['data']['evfile'])
    ft2fits = fits.open(gta.config['data']['scfile'])
    TSTART = float(ft1fits['EVENTS'].header["TSTART"]) - 30.
    TSTOP = float(ft1fits['EVENTS'].header["TSTOP"]) + 30.
    nevt = int(ft2fits['SC_DATA'].header['NAXIS2']) if nevt_max is None else nevt_max
    src = gta.roi[gta.config['selection']['target']]
    c0 = SkyCoord(ra = src.radec[0],dec =  src.radec[1], unit = 'deg')

    # the GTIs from the FT1 file
    ft1_gti_start = ft1fits['GTI'].data.field('START')
    ft1_gti_stop = ft1fits['GTI'].data.field('STOP')

    # central times of GTI in SC file
    tcen = 0.5 * (ft2fits['SC_DATA'].data.field("START") + \
                    ft2fits['SC_DATA'].data.field("STOP"))
    mask = (tcen > TSTART) & (tcen < TSTOP)

    # Spacecraft coordinates
    cz = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_SCZ'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_SCZ'),
                    unit = 'deg')
    cx = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_SCX'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_SCX'),
                    unit = 'deg')
    czen = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_ZENITH'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_ZENITH'),
                    unit = 'deg')

    sep_z0 = c0.separation(cz)
    sep_zenith0 = c0.separation(czen)

    id = np.cos(np.radians(sep_z0))

    # compute phi direction to src in SC coordinates
    phi = get_phi(cz = cz, cx = cx, c0 = c0)

    # livetime
    livetime = ft2fits['SC_DATA'].data.field('LIVETIME')
    
    # loop over GTIs in SC file
    af, ab = np.zeros(nevt), np.zeros(nevt)
    logging.info("Calculating exposure for E = {0:.2f} MeV".format(energy))

    mask = mask & (sep_zenith0.value <= gta.config['selection']['zmax'])
    t1 = time.time()
    for i in tqdm(range(nevt)):
        if not mask[i]: continue

        #igti = 0

        # find the right GTI for the considered SC GTI 
        #while igti < ft1_gti_start.size and ft1_gti_start[igti] < tcen[i]:
            #igti += 1
            #if igti >= ft1_gti_start.size: 
                #break
        m1 = tcen[i] >= ft1_gti_start
        m2 = tcen[i] < ft1_gti_stop
        if not (m1 & m2).sum():
            continue

#        if tcen[i] > ft1_gti_start[igti - 1] and tcen[i] < ft1_gti_stop[igti - 1]:
#            af[i] = aeff_fa.value(energy,
#                            float(sep_z0[i].value),
#                            float(phi[i])) *livetime[i]
#            ab[i] = aeff_ba.value(energy,
#                            float(sep_z0[i].value),
#                            float(phi[i])) *livetime[i]
        af[i] = aeff_fa.value(energy,
                        float(sep_z0[i].value),
                        float(phi[i])) *livetime[i]
        ab[i] = aeff_ba.value(energy,
                        float(sep_z0[i].value),
                        float(phi[i])) *livetime[i]

    print ("it took {0:.2f}s".format(time.time() - t1))
    logging.info("Done")
    return tcen,af,ab

def comp_exposure_phi(gta, energy = None, nevt_max=None):
    """
    Compute the approximate exposure 

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        The analysis object

    {options}

    energy: float or None
        The energy at which the exposure is evaluated.
        If none, central energy from selection is used
    """

    # load the effective areas
    irf_loader.Loader_go()
    irfs = gta.config['gtlike']['irfs']
    factory = irf_loader.IrfsFactory_instance()
    fronta = factory.create(irfs + "::FRONT")
    backa = factory.create(irfs + "::BACK")
    aeff_fa = fronta.aeff()
    aeff_ba = backa.aeff()

    if type(energy) == type(None):
        energy = float(np.sqrt(gta.config['selection']['emin'] * \
                            gta.config['selection']['emax']))

    # load the FT1 and FT2 files
    ft1fits = fits.open(gta.config['data']['evfile'])
    ft2fits = fits.open(gta.config['data']['scfile'])
    TSTART = float(ft1fits['EVENTS'].header["TSTART"]) - 30.
    TSTOP = float(ft1fits['EVENTS'].header["TSTOP"]) + 30.
    nevt = int(ft2fits['SC_DATA'].header['NAXIS2']) if nevt_max is None else nevt_max
    src = gta.roi[gta.config['selection']['target']]
    c0 = SkyCoord(ra = src.radec[0],dec =  src.radec[1], unit = 'deg')

    # the GTIs from the FT1 file
    ft1_gti_start = ft1fits['GTI'].data.field('START')
    ft1_gti_stop = ft1fits['GTI'].data.field('STOP')

    # central times of GTI in SC file
    tcen = 0.5 * (ft2fits['SC_DATA'].data.field("START") + \
                    ft2fits['SC_DATA'].data.field("STOP"))
    mask = (tcen > TSTART) & (tcen < TSTOP)

    # Spacecraft coordinates
    cz = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_SCZ'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_SCZ'),
                    unit = 'deg')
    cx = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_SCX'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_SCX'),
                    unit = 'deg')
    czen = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_ZENITH'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_ZENITH'),
                    unit = 'deg')

    sep_z0 = c0.separation(cz)
    sep_zenith0 = c0.separation(czen)

    id = np.cos(np.radians(sep_z0))

    # compute phi direction to src in SC coordinates
    phi = get_phi(cz = cz, cx = cx, c0 = c0)

    # livetime
    livetime = ft2fits['SC_DATA'].data.field('LIVETIME')
    
    # loop over GTIs in SC file
    af, ab = np.zeros(nevt), np.zeros(nevt)
    logging.info("Calculating exposure for E = {0:.2f} MeV".format(energy))

    mask = mask & (sep_zenith0.value <= gta.config['selection']['zmax'])
    t1 = time.time()
    for i in range(nevt):
        if not mask[i]: continue

        igti = 0

        # find the right GTI for the considered SC GTI 
        while igti < ft1_gti_start.size and ft1_gti_start[igti] < tcen[i]:
            igti += 1
            if igti >= ft1_gti_start.size: 
                break

        if tcen[i] > ft1_gti_start[igti - 1] and tcen[i] < ft1_gti_stop[igti - 1]:
            af[i] = aeff_fa.value(energy,
                            float(sep_z0[i].value),
                            float(phi[i])) *livetime[i]
            ab[i] = aeff_ba.value(energy,
                            float(sep_z0[i].value),
                            float(phi[i])) *livetime[i]

    print ("it took {0:.2f}s".format(time.time() - t1))
    logging.info("Done")
    return tcen,af,ab

def comp_exposure_phi_v2(gta, energy = None):
    """
    Compute the approximate exposure 

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        The analysis object

    {options}

    energy: float or None
        The energy at which the exposure is evaluated.
        If none, central energy from selection is used
    """

    # load the effective areas
    irf_loader.Loader_go()
    irfs = gta.config['gtlike']['irfs']
    factory = irf_loader.IrfsFactory_instance()
    fronta = factory.create(irfs + "::FRONT")
    backa = factory.create(irfs + "::BACK")
    aeff_fa = fronta.aeff()
    aeff_ba = backa.aeff()

    if type(energy) == type(None):
        energy = float(np.sqrt(gta.config['selection']['emin'] * \
                            gta.config['selection']['emax']))

    # load the FT1 and FT2 files
    ft1fits = fits.open(gta.config['data']['evfile'])
    ft2fits = fits.open(gta.config['data']['scfile'])
    TSTART = float(ft1fits['EVENTS'].header["TSTART"]) - 30.
    TSTOP = float(ft1fits['EVENTS'].header["TSTOP"]) + 30.
    nevt = int(ft2fits['SC_DATA'].header['NAXIS2'])
    src = gta.roi[gta.config['selection']['target']]
    c0 = SkyCoord(ra = src.radec[0],dec =  src.radec[1], unit = 'deg')

    # the GTIs from the FT1 file
    ft1_gti_start = ft1fits['GTI'].data.field('START')
    ft1_gti_stop = ft1fits['GTI'].data.field('STOP')

    # central times of GTI in SC file
    tcen = 0.5 * (ft2fits['SC_DATA'].data.field("START") + \
                    ft2fits['SC_DATA'].data.field("STOP"))
    mask = (tcen > TSTART) & (tcen < TSTOP)

    # Spacecraft coordinates
    cz = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_SCZ'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_SCZ'),
                    unit = 'deg')
    cx = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_SCX'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_SCX'),
                    unit = 'deg')
    czen = SkyCoord(ra = ft2fits['SC_DATA'].data.field('RA_ZENITH'),
                    dec = ft2fits['SC_DATA'].data.field('DEC_ZENITH'),
                    unit = 'deg')

    sep_z0 = c0.separation(cz)
    sep_zenith0 = c0.separation(czen)

    id = np.cos(np.radians(sep_z0))

    # compute phi direction to src in SC coordinates
    phi = get_phi(cz = cz, cx = cx, c0 = c0)

    # livetime
    livetime = ft2fits['SC_DATA'].data.field('LIVETIME')
    
    # loop over GTIs in SC file
    af, ab = np.zeros(nevt), np.zeros(nevt)
    logging.info("Calculating exposure for E = {0:.2f} MeV".format(energy))

    mask = mask & (sep_zenith0.value <= gta.config['selection']['zmax'])

    ft1_gti_bins = np.array([ft1_gti_start, ft1_gti_stop])

    igti = np.digitize(tcen, ft1_gti_start)

    # make 1d array with start0,stop0,start1,stop0,...
    ft1_gti_bins = ft1_gti_bins.T.reshape(-1)
    # check for all times of GTI file in which bin they are. 
    # if even, then inside a GTI if odd then outside 

    #igti = np.digitize(tcen, ft1_gti_bins)
    #mask = (igti > 0) & (igti < len(ft1_gti_bins)) & \
            #mask & np.invert((igti - 1 % 2).astype(np.bool))

    t1 = time.time()
    for i in range(nevt):
        if not mask[i]: continue

        if tcen[i] > ft1_gti_start[igti[i] - 1] and tcen[i] < ft1_gti_stop[igti[i] - 1]:
            af[i] = aeff_fa.value(energy,
                            float(sep_z0[i].value),
                            float(phi[i])) *livetime[i]
            ab[i] = aeff_ba.value(energy,
                            float(sep_z0[i].value),
                            float(phi[i])) *livetime[i]

    print ("it took {0:.2f}s".format(time.time() - t1))
    logging.info("Done")
    return tcen,af,ab

def get_phi(cz, cx, c0):
    """Calculate Phi direction to source given space craft direction

    Parameters
    ----------
    cz: `~astropy.coordinates.SkyCoord`
        coordinates of space craft z direction
    cx: `~astropy.coordinates.SkyCoord`
        coordinates of space craft x direction
    c0: `~astropy.coordinates.SkyCoord`
        coordinates of source

    Returns
    -------
    """
    vx_z = np.cos(np.radians(cz.dec.value)) * \
            np.cos(np.radians(cz.ra.value))
    vy_z = np.cos(np.radians(cz.dec.value)) * \
            np.sin(np.radians(cz.ra.value))
    vz_z = np.sin(np.radians(cz.dec.value))

    vx_0 = np.cos(np.radians(c0.dec.value)) * \
            np.cos(np.radians(c0.ra.value))
    vy_0 = np.cos(np.radians(c0.dec.value)) * \
            np.sin(np.radians(c0.ra.value))
    vz_0 = np.sin(np.radians(c0.dec.value))

    vx_x = np.cos(np.radians(cx.dec.value)) * \
            np.cos(np.radians(cx.ra.value))
    vy_x = np.cos(np.radians(cx.dec.value)) * \
            np.sin(np.radians(cx.ra.value))
    vz_x = np.sin(np.radians(cx.dec.value))

    z0_vect_x = vz_0 * vy_z - vz_z * vy_0
    z0_vect_y = vz_z * vx_0 - vz_0 * vx_z
    z0_vect_z = vy_0 * vx_z - vy_z * vx_0

    sin_phi = z0_vect_x * vx_x + \
                z0_vect_y* vy_x + \
                z0_vect_z * vz_x

    zx_vect_x = vz_x * vy_z - vz_z * vy_x
    zx_vect_y = vz_z * vx_x - vz_x * vx_z
    zx_vect_z = vy_x * vx_z - vy_z * vx_x

    cos_phi = zx_vect_x * z0_vect_x + \
                zx_vect_y * z0_vect_y + \
                z0_vect_z * zx_vect_z

    phi = np.degrees(np.arccos(cos_phi))

    m = sin_phi > 0.
    phi[m] = 360. - phi[m]
    return phi

def Get_flux_unc(index,glon,glat,tstart,tstop,
                    tcen, exp, 
                    El,normeg,normgal,indexgal, nevt):
    """
    Get uncertainty on flux
    """
    # get the times of exposure within tstart and tstop
    m = (tcen >= tstart) & ( tcen < tstop)
    # add up the total exposure in this time range
    exp_tot = float(exp[m].sum())
    # output array
    out=array('d',[0])                 

    ff=Get_flux_mult_src_minuit(float(index),float(glon),float(glat),
                    float(normeg),float(normgal),float(indexgal),
                    float(tstart),float(tstop),float(exp_tot),out)
    if (ff<1.1e-9):
        return -1,-1
    unc=Get_unc_mult(index,ff,glon,glat,
                tstart,tstop,nevt,
                normeg,normgal,indexgal,El)

    logging.debug("{0} {1}".format(ff, unc))
    return ff,unc


#def time_bins(gta, macroC, tstart, tcen, front, back, rev = 0, sources = None):
def time_bins(gta, tcen, exp,
    critval = 20., 
    tstart = 0, rev = 0, sources = None, crit = 1,
    tstop = 1e40, 
    forcecatalog = False,
    normeg=1, normgal=1, indexgal=0, Epivot = None):
    """
    Run the adaptive time binning

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis`
        analysis object
    tcen: `~numpy.ndarray`
        Times for exposure calculation
    exp: `~numpy.ndarray`
        exposure (either front or back)

    {options}

    critval: float
        if crit == 0: value of TS, elif crit == 1, value of 
        desired flux uncertainty in percent. Default: 20.
    crit: int
        if 0, based on constant TS, if 1, based on constant flux unc.
        Default: 1.
    tstart: float
        Alternative starting value in MET,
        maximum between TSTART in ft1 in this value is used
    tstop: float
        Alternative stop value in MET,
        minimum between TSTOP in ft1 and this value is used
    Epivot: float
        Pivot energy. if not given, it will be calculated. 
        gta.optimize has to be run prior to this.
    sources: None or list
        If none, calculate light curve for central source, 
        otherwise for sources whose names match the gta source names.
        First source in list will be used.
    forcecatalog: bool
        if true, use catalog values for flux and index
    """

    irfs = gta.config['gtlike']['irfs']

    # get the source coordinates and the flux 
    # and index(es)
    if type(Epivot) == type(None):
        Epivot = comp_E1(gta)
        logging.info("Epivot = {0} MeV".format(Epivot))
        if np.isnan(Epivot):
            Epivot = np.max([300.,gta.config['selection']['emin']])
            logging.error("Epivot is nan, " \
                    "using {0:.2f} MeV instead!".format(Epivot))

    if type(sources) == type(None):
        src = gta.roi[gta.config['selection']['target']]
        c = [SkyCoord(ra = src.radec[0],dec =  src.radec[1], unit = 'deg')]
        if np.isnan(src['flux']) or forcecatalog:
            flux = array('d',[src['catalog']['Flux1000']])
            index = array('d',[src['catalog']['Spectral_Index']])
        else:
            flux = array('d',[src['flux']])
            index = array('d',[src['dnde_index']])
        glon_all = array('d', [c[-1].galactic.l.value])
        glat_all = array('d', [c[-1].galactic.b.value])
    elif type(sources) == list:
        if len(sources) > 2:
            raise Exception("too many sources given, max is 2!")
        c = []
        index = []
        flux = []
        glon_all = []
        glat_all = []
        for s in sources:
            src = gta.roi[s]
            c.append(SkyCoord(ra = src.radec[0],dec =  src.radec[1], unit = 'deg'))
            if np.isnan(src['flux']):
                flux.append(src['catalog']['Flux1000'])
                index.append(src['catalog']['Spectral_Index'])
            else:
                flux.append(src['flux'])
                index.append(src['dnde_index'])
            irfs.append(gta.config['gtlike']['irfs'])

            glon_all.append(c[-1].galactic.l.value)
            glat_all.append(c[-1].galactic.b.value)

        flux = array('d',flux)
        index = array('d',index)
    else:
        raise Exception("sources keyword not understood")


    logging.info("Source properties: glon = {0}, glat = {1}, flux = {2}, index = {3}".format(
                    glon_all, glat_all, flux, index))
    ft1fits = fits.open(gta.config['data']['evfile'])
    ft1 = Table.read(gta.config['data']['evfile'], hdu = 'EVENTS')

    roicut = ft1fits['EVENTS'].header['DSVAL2']
    nevt = ft1fits['EVENTS'].header['NAXIS2']

    irfs = gta.config['gtlike']['irfs']
    if (irfs=="P8R2_SOURCE_V6"):
        initialize_TS(path.dirname(fermiAnalysis.__file__) + '/')  
    else:
        raise Exception("Only IRF P8R2_SOURCE_V6 implemented!")

    if nevt > get_evtlim():
        raise Exception("number of events exceeds limits in TS_estimate (evtlim = {0:n})".format(get_evtlim()))

    tstart = np.max([tstart, ft1fits['EVENTS'].header['TSTART']])
    Tstop = np.min([tstop, ft1fits['EVENTS'].header['TSTOP']])
    logging.info("Using start time {0:.2f}".format(tstart))

    # reversed time arrow?
    if rev:
        ft1 = ft1[::-1]

    if len(c) > 1:
        sep_all = np.empty(ft1['RA'].data.size * len(c))
        sep = []

        for i, ci in enumerate(c):
            sep.append(ci.separation(SkyCoord(ra = ft1['RA'].data, dec = ft1['DEC'].data, unit = 'deg')).value)

        sep_all[::2] = sep[0]
        sep_all[1::2] = sep[1]
    else:
        sep_all = c[0].separation(SkyCoord(ra = ft1['RA'].data, dec = ft1['DEC'].data, unit = 'deg')).value
    sep_all = np.radians(sep_all)

    vglob = Glob_mult(array('d',ft1['ENERGY'].data),
                array('d', ft1['TIME'].data),
                array('d', ft1['CONVERSION_TYPE'].data),
                array('d', sep_all),
                ft1['TIME'].data.size, len(c))

    feed_sources(glon_all,glat_all,index,flux,len(c))    

    tstop=tstart
    stop=0
    logging.info("tstop={0}".format(tstop))
    logging.info("Tstop={0}".format(Tstop))
    end=0
    #m=0
    ff=flux[0]

    result = {}
    result['tstart'] = []
    result['tstop'] = []
    result['flux'] = []
    result['dflux'] = []
    result['index'] = []
    result['error'] = []

    # bins with constant TS
    if (crit==0):
        TS0=critval
        scrit="TS"
        vcrit=TS0
        raise Exception("Constant TS not implemented")
    # bins with constant uncertainty
    else:
        unc0=critval
        scrit="unc"
        vcrit=unc0
        while (tstop<Tstop and end==0):
            stop=0
            tmax=Tstop
            tmin=tstart
            ff0=ff
            kj=0
            unc1=unc0
            iav=0
            uncmin=100
            uncmax=0
            while (stop==0 and kj<20):
                if (iav==0):
                    tstop=Get_time_unc_mult(
                            unc0,index[0],ff,
                            glon_all[0],glat_all[0],
                            tstart,normeg,normgal,indexgal,Epivot)
                else:    
                    iav=0
                    tstop=(tmin+tmax)/2.
                    logging.debug("{0}".format([tmin,uncmin,tmax,uncmax]))
                    logging.info("in patch iav")
                logging.info("tstop {0} {1}".format(tstop, ff))
                if (tstop<0):
                    ff,unc=Get_flux_unc(index[0],glon_all[0],glat_all[0],
                                        tstart,Tstop,tcen,exp,
                                        Epivot,normeg,normgal,indexgal, nevt)
                    logging.info("ff = {0}, unc = {1}".format(ff,unc))
                    if (ff<0 or unc>unc0):
                        end=1
                        tstop=Tstop
                        break
                    tstop=(tmin+tmax)/2.
                if (tstop>0 and tstop<=tmin):
                    logging.info("taking mid time 1")
                    stop=(tmin+tmax)/2.
                if (tstop>0 and tstop>=tmax):
                    logging.info("taking mid time 2")
                    tstop=(tmin+tmax)/2.     
                ff,unc=Get_flux_unc(index[0],glon_all[0],glat_all[0],
                            tstart,tstop,tcen, exp,
                            Epivot,normeg,normgal,indexgal, nevt)
                ki=0
                while (ff<0 and ki<10):
                    tstop=(tmax+tmin)/2.
                    ff,unc=Get_flux_unc(index[0],glon_all[0],glat_all[0],
                            tstart,tstop,tcen, exp,
                            Epivot,normeg,normgal,indexgal, nevt)
                    if (ff<0):
                        tmin=tstop
                    ki += 1
                if (ki==10):
                    logging.error("ki reached 10")
                    break
                if (unc>unc0):
                    if (tstop>tmin):
                        tmin=tstop
                        uncmin=unc
                        logging.info("new tmin: {0} {1}".format(tmin,unc))
                else:
                    if (tstop<tmax):   
                        tmax=tstop
                        uncmax=unc
                        logging.info("new tmax: {0} {1}".format(tmax,unc))
                df= 2*np.abs(ff-ff0)/(ff+ff0)
                logging.info("unc, df: {0} {1}".format(unc,df))
                kj=kj+1
                if (np.abs(unc-unc1)<0.02):
                    iav=1
                    logging.info("iav set to 1")
                if (np.abs(unc-unc0)<1.0 or (tmax-tmin)<10):
                    #if (np.abs(unc-unc0)<0.5 or np.abs(unc-unc1)<0.05 ):
                    stop=1
                    break
                unc1=unc
                ff0=ff
            if (end==1):
                break
            logging.info("unc={0}".format(unc))
            logging.info("tstart={0}".format(tstart))
            logging.info("tstop={0}".format(tstop))
            logging.info("ff={0}".format(ff))
            error=0
            result['tstart'].append(tstart)
            result['tstop'].append(tstop)
            result['flux'].append(ff)
            result['dflux'].append(ff * unc / 100.)
            result['index'].append(index[0])
            result['error'].append(error)

            # continue with next step
            tstart=tstop

    result['tstart'] = np.array(result['tstart'])
    result['tstop']  = np.array(result['tstop'])
    result['flux']   = np.array(result['flux'])
    result['dflux']  = np.array(result['dflux'])
    result['index']  = np.array(result['index'])
    result['error']  = np.array(result['error'])
    return result
