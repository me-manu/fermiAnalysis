import numpy as np
import logging
import matplotlib.pyplot as plt


def hop(x):
    """
    Run HOP algorithm
    
    Parameters
    ----------
    x: `~numpy.ndarray`
        array with density (fluxes, ...)
    
    Returns
    -------
    Three `~numpy.ndarrays` with ids of peaks,
    the working id_hop vector, and a list with the 
    ids of the sheds

    Notes
    -----
    Adopted from Jeff Scargle's implementation 
    in MatLab
    """
    id_hop_vec = range(len(x))
    
    # Exceptions for end points:

    if x[1] > x[0]:
        id_hop_vec[0] = 1

    if x[-2] > x[-1]:
        id_hop_vec[-1] = len(x) - 2
    
    for ii in range(2, len(x) - 1):

        x_left = x[ii-1]
        x_this = x[ii]
        x_rite = x[ii+1]
        
        id_next = np.argmax([x_left, x_rite])
        x_next = [x_left, x_rite][id_next]
    
        if x_next > x_this:
            if not id_next:
                id_hop_vec[ii] = ii - 1
            else:
                id_hop_vec[ii] = ii + 1
    
    #===================================
    #         now do HOP
    #===================================
    id_hop_vec = np.array(id_hop_vec)
    id_hop_work = np.array(id_hop_vec)
    id_hop_old = np.array(id_hop_work)


    while True:
    
    # update pointers via local hill climbing
        
        id_hop_work = id_hop_vec[ id_hop_work ]
        delt = np.sum( id_hop_work - id_hop_old )
        if delt == 0:
            break
        id_hop_old = id_hop_work
    
    id_peaks_vec = np.unique(id_hop_work)
    id_sheds = []
    for ii,idp in enumerate(id_peaks_vec):
        id_sheds.append(np.where(id_hop_work == idp)[0])
    return id_peaks_vec,id_hop_work, id_sheds

def data_in_sheds(bbf,xmin, xmax, xref, y, id_sheds, yerr = None):
    """
    Get raw data points in bayesian block belonging to one shed 
    
    Parameters
    ----------
    bbf: `~numpy.ndarray`
        Bayesian block (BB) ids
    
    xmin: `~numpy.ndarray`
        left bin bounds of x values 
    
    xmax: `~numpy.ndarray`
        right bin bounds of x values 
    
    xref: `~numpy.ndarray`
        central x value of data bins 
    
    y: `~numpy.ndarray`
        y value of data points 
    
    id_sheds: list
        list containing the BB ids of the sheds
    
    yerr: `~numpy.ndarray` or None (optional)
        y errors value of data points 
        
    Returns
    -------
    list with x and y values for each shed
    """
    data_x = []
    data_y = []
    if not type(yerr) == type(None):
        data_yerr = []
    
    for ip, idp in enumerate(id_sheds):
        if len(idp) == 1:
            continue
        if ip < len(id_sheds) - 1:
            x = np.concatenate([xmin[bbf[idp]],[xmin[bbf[idp[-1]+1]]]] )
        else:
            x = np.concatenate([xmin[bbf[idp]],[xmax[-1]]] )
        m = (xref >= x.min()) & (xref < x.max())
        data_x.append(xref[m])
        data_y.append(y[m])
        if not type(yerr) == type(None):
            data_yerr.append(yerr[m])
    if not type(yerr) == type(None):
        return data_x,data_y,data_yerr
    else:
        return data_x,data_y

def plot_sheds(xmin, xmax, y, id_sheds, bbf = None, ax = None, 
               col = [plt.cm.tab10(0.8),plt.cm.tab10(1.0)], y2 = 0.,
               hatch = ['////','\\\\\\\\'], exp = 0., **kwargs):
    """
    Plot the sheds of flares
    
    Parameters
    ----------
    xmin: `~numpy.ndarray`
        left bin bounds of x values 
    
    xmax: `~numpy.ndarray`
        right bin bounds of x values 
    
    id_sheds: `~numpy.ndarray`
        id vector from hop algorithm

    y: `~numpy.ndarray`
        y value of data points / of BB
       
    bbf: `~numpy.ndarray` or None (optional)
        Bayesian block (BB) ids
    """
    if type(ax) == type(None):
        ax = plt.gca()
    for ip, idp in enumerate(id_sheds):
        if len(idp) == 1:
            continue
    
        if ip < len(id_sheds) - 1:
            if not type(bbf) == type(None):
                x = np.concatenate([xmin[bbf[idp]],[xmin[bbf[idp[-1]+1]]]] )
            else:
                x = np.concatenate([xmin[idp],[xmin[idp[-1]+1]]] )
        else:
            if not type(bbf) == type(None):
                x = np.concatenate([xmin[bbf[idp]],[xmax[-1]]] )
            else:
                x = np.concatenate([xmin[idp],[xmax[-1]]] )
        f = np.concatenate([y[idp],[y[idp[-1]]]]) / 10.**exp    
        if hatch is None:
            ax.fill_between(x,f,
                     y2 = y2,
                     color = col[ip % len(col)],
                     step = 'post', **kwargs)
        else:
            ax.fill_between(x,f,
                     y2 = y2, 
                     edgecolor = col[ip % len(col)],
                     hatch = hatch[ip % len(hatch)],
                     facecolor = 'none',
                     step = 'post', **kwargs)
    return

def hop_flares(xmin, xmax, y, id_sheds, threshold_flux, bbf = None, min_flux = None, extend = 0., ts = None):
    """
    Determine flare start and end time
    from HOP sheds above some threshold flux
    
    Parameters
    ----------
    xmin: `~numpy.ndarray`
        left bin bounds of x values 
    
    xmax: `~numpy.ndarray`
        right bin bounds of x values 
    
    id_sheds: `~numpy.ndarray`
        id vector from hop algorithm

    y: `~numpy.ndarray`
        y value of data points / of BB

    threshold_flux: float
        threshold flux above which HOPs are selected
        (Or TS value if ts array is provided)
       
    bbf: `~numpy.ndarray` or None (optional)
        Bayesian block (BB) ids

    min_flux: float or None (optional)
        if given, include only data points above this flux in 
        the flare time range

    extend : float, 
        extend each flare range by +/- this value
    """
    if not type(min_flux) == type(None):
        if min_flux > threshold_flux:
            raise ValueError("Min flux > threshold flux!")

    xhop, fhop_start, fhop_stop = [],[], []
    for ip, idp in enumerate(id_sheds):
        if len(idp) == 1:
            continue
    
        if ip < len(id_sheds) - 1:
            if not type(bbf) == type(None):
                #x = np.concatenate([xmin[bbf[idp]],[xmin[bbf[idp[-1]+1]]]] )
                x = np.concatenate([xmin[bbf[idp]],[xmin[bbf[idp[-1]+1]]]] )
            else:
                x = np.concatenate([xmin[idp],[xmin[idp[-1]+1]]] )
        else:
            if not type(bbf) == type(None):
                x = np.concatenate([xmin[bbf[idp]],[xmax[-1]]] )
            else:
                x = np.concatenate([xmin[idp],[xmax[-1]]] )
        f = np.concatenate([y[idp],[y[idp[-1]]]])
        fhop_start.append(np.concatenate([y[idp],[y[idp[-1]]]]))
        fhop_stop.append(np.concatenate([[y[idp[0]]],y[idp]]))
        xhop.append(x)
    # select hops with f > threshold
    mask = np.array([(f >= threshold_flux).sum() for f in fhop_start]) > 0
    tflare = []
    for i,m in enumerate(mask):
        if m:
            if type(min_flux) == type(None):
                xhopi = xhop[i]
            else:
                #xhopi = xhop[i][fhop[i] >= min_flux]
                xhopi = xhop[i][(fhop_start[i] >= min_flux) | (fhop_stop[i] >= min_flux)] 

            if not xhopi.size:
                #raise ValueError("HOP interval has no entries")
                logging.warning("No flux values above min flux! Trying next HOP.")
                continue

            if len(tflare):
                if xhopi.min() - extend <= tflare[-1][-1]:
                    tflare[-1][-1] = xhopi.max() + extend # extend time range
                else:
                    tflare.append([xhopi.min() - extend , xhopi.max() + extend])
            else:
                tflare.append([xhopi.min() - extend, xhopi.max() + extend])
    return np.array(tflare)

def estimate_time(f,dt,fpeak):
    """Estimate rise / decay time"""
    return (f * dt).sum() / fpeak

def calc_flare_duration_rise_decay(xmin, xmax, y, id_sheds, bbf = None):
    """
    Calculate flare duration, rise, and decay times 
    using the flares defined with HOP algorithm
    and simple formulas (no fitting of flare profile)
    
    Parameters
    ----------
    xmin: `~numpy.ndarray`
        left bin bounds of x values 
    
    xmax: `~numpy.ndarray`
        right bin bounds of x values 
    
    xref: `~numpy.ndarray`
        central x value of data bins 
    
    y: `~numpy.ndarray`
        y value of data points / of BB
       
    bbf: `~numpy.ndarray` or None (optional)
        Bayesian block (BB) ids
    """
    tr, td, dtflare, fpeak,fintegral = [],[],[],[],[]

    # loop through water sheds
    for ip, idp in enumerate(id_sheds):
        if len(idp) == 1:
            continue
    
        if ip < len(id_sheds) - 1:
            if not type(bbf) == type(None):
                x = np.concatenate([xmin[bbf[idp]],[xmin[bbf[idp[-1]+1]]]] )
            else:
                x = np.concatenate([xmin[idp],[xmin[idp[-1]+1]]] )
        else:
            if not type(bbf) == type(None):
                x = np.concatenate([xmin[bbf[idp]],[xmax[-1]]] )
            else:
                x = np.concatenate([xmin[idp],[xmax[-1]]] )

        dt = np.diff(x)
        f = y[idp]
        tcen = 0.5 * (x[:-1] + x[1:])
        idmax = np.argmax(f)

        if idmax == f.size - 1:
        #    trise = estimate_time(f, dt, f[idmax])
            idmin = np.argmin(f)
            trise = tcen[idmax] - tcen[idmin]
            tdecay = 0. # no decay

        elif idmax == 0:
            idmin = np.argmin(f)
            trise = 0.
            #tdecay = estimate_time(f, dt, f[idmax])
            tdecay = tcen[idmin] - tcen[idmax]

        else:
            #trise = estimate_time(f[:idmax + 1], dt[:idmax + 1], f[idmax])
            #tdecay = estimate_time(f[idmax:], dt[idmax:], f[idmax])
            idmin_pre = np.argmin(f[:idmax])
            idmin_post = np.argmin(f[idmax:])
            trise = tcen[idmax] - tcen[:idmax][idmin_pre]
            tdecay = tcen[idmax:][idmin_post] - tcen[idmax] 


        tr.append(trise)
        td.append(tdecay)
        dtflare.append(dt.sum())
        fpeak.append(f[idmax])
        fintegral.append((f * dt).sum())

    return np.array(tr), np.array(td), np.array(dtflare), np.array(fpeak), np.array(fintegral)
