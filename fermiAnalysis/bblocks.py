import numpy as np
import logging
from scipy.interpolate import interp1d
from astropy.table import Table

mjd_to_met = lambda mjd: (mjd - 54682.65) * 86400. +  239557414.0
met2mjd = lambda met: 54682.65 + (met - 239557414.0) / (86400.)
mjd2met = lambda mjd: (mjd - 54682.65) * 86400. + 239557414.0


def nlogn (n, dt):
    # I really feel like there must be a cleverer way to do this
    # scalar-or-vector possible-bad-value masking.

    if np.isscalar (n):
        if n == 0:
            return 0.
        return n * (np.log (n) - np.log (dt))

    n = np.asarray (n)
    mask = (n == 0)
    r = n * (np.log (np.where (mask, 1, n)) - np.log (dt))
    
    return np.where (mask, 0, r)

def binblock_unbinned(times, exp = np.array([]), p0 = 0.05, itr = 1):
    """
    Wrapper to run the binblock algorith in unbinned mode
    """
    widths = np.diff(times)
    counts = np.ones_like(widths)

    if exp.size: 
        expw = exp[:-1]
    else:
        expw = exp

    bb, avg = binbblock (widths, counts, p0=p0, itr = itr, datatype = 'unbinned',
                            exp =expw)

    return np.concatenate([times[:-1][bb], [times[-1]]])

def binbblock (widths, counts, p0=0.05, itr = 1, datatype = 'unbinned',
                errs = np.array([]), exp = np.array([]), return_davg = False):
    """
    Calculate Bayesian blocks 
    
    Parameters
    ----------
    widths: `~numpy.ndarray`
        n-dim array with bin widths
    counts: `~numpy.ndarray`
        n-dim array with bin contents
    
    kwargs
    ------
    datatype: string
        either 'unibinned', 'binned', or 'point', determines 
        type of histogram
    errs: `~numpy.ndarray`
        if datatype = point, gives measurement error in each bin
    exp: `~numpy.ndarray`
        Gives exposure in each bin / for each event
    p0: float
        false alarm probability for new change point, 
        relevant for prior
    itr: int
        number of iterations, each iteration changes 
        p0 to p0^1/(number of change points)
    return_davg: bool
        if true and datatype == point,
        also return the uncertainty in the average flux
        (default: false)
    
    Returns
    -------
    tuple with `~numpy.ndarrays`
    containing the block starts and the average block content
    """
    widths = np.asarray (widths)
    counts = np.asarray (counts)
    ncells = widths.size
    origp0 = p0

    if np.any (widths <= 0):
        raise ValueError ('bin widths must be positive')
    if widths.size != counts.size:
        raise ValueError ('widths and counts must have same size')
    if p0 < 0 or p0 >= 1.:
        raise ValueError ('p0 must lie within [0, 1)')
        
    if datatype == 'unbinned':
        prior = lambda p0: 4. - np.log (p0 / (0.0136 * ncells**0.478))
        fitness = lambda n,t: nlogn(n,t)
    elif datatype == 'point':
        if errs.size != counts.size:
            raise ValueError ('errors and counts must have same size')
        prior = lambda p0: 1.32 + 0.577 * np.log10(ncells) # should change with p0 not implemented 
        fitness = lambda a,b: b**2. / 4. / a

    elif datatype == 'binned':
        prior = lambda p0: 4. - np.log (p0 / (0.0136 * ncells**0.478))
        fitness = lambda n,t: nlogn(n,t)
    else:
        raise ValueError ("datatype must be either 'binned', 'unbinned', or 'point'.")

    if exp.size:
        counts /= exp

    if not datatype == 'point':
        vedges = np.cumsum (np.concatenate (([0], widths))) # size: ncells + 1
        block_remainders = vedges[-1] - vedges # size: nedges = ncells + 1
        ccounts = np.cumsum (np.concatenate (([0], counts)))
        count_remainders = ccounts[-1] - ccounts

    prev_blockstarts = None
    best = np.zeros (ncells, dtype=np.float)
    last = np.zeros (ncells, dtype=np.int)
    
    # iterate over prior
    for _ in xrange (itr):
        ncp_prior = prior(p0)
        #print "Prior on blocks is ", ncp_prior

        for r in xrange (ncells):
            if datatype == 'point':
                sum_x_1 = np.cumsum( (counts / errs / errs)[:r+1][::-1] )[::-1] # sum( x / sig^2 ) = b
                sum_x_0 = np.cumsum( (1. / errs / errs)[:r+1][::-1] )[::-1] # sum( 1 / sig^2 ) = a
                fit_vec = fitness (sum_x_0, sum_x_1) 
            else:
                tk = block_remainders[:r+1] - block_remainders[r+1]
                nk = count_remainders[:r+1] - count_remainders[r+1]

                try:
                    fit_vec = fitness (nk, tk)
                except:
                    raise

            # This incrementally penalizes partitions with more blocks:
            tmp = fit_vec - ncp_prior
            tmp[1:] += best[:r]

            imax = np.argmax (tmp)
            last[r] = imax
            best[r] = tmp[imax]


        # different semantics than Scargle impl: our blockstarts is similar to
        # their changepoints, but we always finish with blockstarts[0] = 0.

        work = np.zeros (ncells, dtype=int)
        workidx = 0
        ind = last[-1]

        while True:
            work[workidx] = ind
            workidx += 1
            if ind == 0:
                break
            ind = last[ind - 1]

        blockstarts = work[:workidx][::-1]

        if prev_blockstarts is not None:
            if (blockstarts.size == prev_blockstarts.size and
                (blockstarts == prev_blockstarts).all ()):
                break # converged

        if blockstarts.size == 1:
            break # can't shrink any farther

        # Recommended ad-hoc iteration to favor fewer blocks above and beyond
        # the value of p0:
        p0 = 1. - (1. - p0)**(1. / (blockstarts.size - 1))
        prev_blockstarts = blockstarts

    # make sure that 0-th block starts
    # with first time interval
    assert blockstarts[0] == 0
    
    # number of new blocks
    nblocks = blockstarts.size
    
    # calculate new counts and widths 
    # of rebinned histogram
    rcounts = np.empty (nblocks, dtype=np.int)
    rwidths = np.empty (nblocks)
    avg = np.empty (nblocks)
    if datatype == 'point':
        davg = np.empty (nblocks)
    
    for iblk in xrange (nblocks):
        cellstart = blockstarts[iblk]
        if iblk == nblocks - 1:
            cellend = ncells - 1
        else:
            cellend = blockstarts[iblk+1] - 1

        rwidths[iblk] = widths[cellstart:cellend+1].sum ()
        rcounts[iblk] = counts[cellstart:cellend+1].sum ()
        avg[iblk] = counts[cellstart:cellend+1].sum () / (cellend + 1 - cellstart)
        if datatype == 'point':
            davg[iblk] = np.sqrt((errs[cellstart:cellend+1] ** 2.).sum ()) / (cellend + 1 - cellstart)

    rates = rcounts / rwidths
    
    if return_davg and datatype == 'point':
        return blockstarts, avg, davg
    else:
        return blockstarts, avg

# ----- an alternative implementation ------------------------ #
def bayesian_blocks(t, p0 = 0.05, exp = np.array([])):
    """Bayesian Blocks Implementation

    By Jake Vanderplas.  License: BSD
    Based on algorithm outlined in http://adsabs.harvard.edu/abs/2012arXiv1207.5578S

    Parameters
    ----------
    t : ndarray, length N
        data to be histogrammed

    Returns
    -------
    bins : ndarray
        array containing the (N+1) bin edges

    {options}
    exp: scalar or `~numpy.ndarray`
        exposure, corrects the counts by 1 / exp, 
        if array, needs to be of length N

    Notes
    -----
    This is an incomplete implementation: it may fail for some
    datasets.  Alternate fitness functions and prior forms can
    be found in the paper listed above.
    """
    # copy and sort the array
    t = np.sort(t)
    N = t.size

    # create length-(N + 1) array of cell edges
    edges = np.concatenate([t[:1],
                            0.5 * (t[1:] + t[:-1]),
                            t[-1:]])
    block_length = t[-1] - edges
    logging.debug(block_length)
    logging.debug(block_length.shape)

    # arrays needed for the iteration
    nn_vec = np.ones(N)
    if exp.size:
         nn_vec /= exp
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)
    
    ncp_prior = 4 - np.log (p0 / (0.0136 * N**0.478))
    #-----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    #-----------------------------------------------------------------
    for K in range(N):
        # Compute the width and count of the final bin for all possible
        # locations of the K^th changepoint
        width = block_length[:K + 1] - block_length[K + 1]
        count_vec = np.cumsum(nn_vec[:K + 1][::-1])[::-1]
        


        # evaluate fitness function for these possibilities
        fit_vec = count_vec * (np.log(count_vec) - np.log(width))
        fit_vec -= ncp_prior  # 4 comes from the prior on the number of changepoints
        fit_vec[1:] += best[:K]
        

        # find the max of the fitness: this is the K^th changepoint
        i_max = np.argmax(fit_vec)
        last[K] = i_max
        best[K] = fit_vec[i_max]

    #-----------------------------------------------------------------
    # Recover changepoints by iteratively peeling off the last block
    #-----------------------------------------------------------------
    change_points =  np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while True:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    change_points = change_points[i_cp:]
    logging.debug(change_points)
    return edges[change_points]

def remove_gti_blocks(block_bins, times, gtistart, gtistop):
    """
    Remove Bayesian block starting times that co-incide 
    with the start or end of a GTI
    
    Parameters
    ----------
    block_bins: `~numpy.ndarray`
    n+1 dim array with the bayesian block bin (the ones 
    you would provide for plt.hist)
    
    times: `~numpy.ndarray`
    event times
    
    gtistart: `~numpy.ndarray`
    array with gtistart times except the first start time
    
    gtistop: `~numpy.ndarray`
    array with gtistop times except the last stop time
    """
    delete = []
    for tmax in gtistop:
        # index of last data point where exposure > 0
        idx = np.where(times == np.sort(times[times < tmax])[-1])[0]
        # check if this point is in bin definitions
        # and if so, remove it
        if times[idx] in block_bins:
            delete.append(np.where(block_bins == times[idx])[0])
    for tmin in gtistart:
        # index of first data point where exposure > 0
        idx = np.where(times == np.sort(times[times > tmin])[0])[0]
        # check if this point is in bin definitions
        # and if so, remove it
        if times[idx] in block_bins:    
            delete.append(np.where(block_bins == times[idx])[0])
    # delete the blocks
    bbins_new = np.delete(block_bins, delete)
    # calculate the likelihood difference
    h = np.histogram(times, bins=block_bins)
    hnew = np.histogram(times, bins=bbins_new)
    dlog = bblocks.nlogn(hnew[0],hnew[1][:1]).sum() - bblocks.nlogn(h[0],h[1][:1]).sum()
    if dlog > 0: 
        return bbins_new
    else:
        return block_bins



def flare_search(flux, tmin, thr_flux, avg_flux,
    min_fl_flux, decr_thr = False, append = 0.):
    """
    Search for flares

    Parameters
    ----------
    flux: `~numpy.ndarray`
        n-dim array with time series of flux in bayesian block, with last flux appended 
    tmin: `~numpy.ndarray`
        n-dim array with start time of bayesian block, with last end time appended
    thr_flux: float
        number by which average flux has to increase so that 
        flux is considered a flare
    avg_flux: float
        Average flux
    min_fl_flux: float
        minimum flux so that block will still be considered a flare
    decr_thr: bool
        if True, decrease threshold if no flare is detected
        default: False
    append: float
        for each flare, extend it by this amount in days. Combine flares that overlap.
        default: 0.
    """
    idx = np.where(flux >= thr_flux * avg_flux)
    logging.debug(idx)
    logging.debug("Extending flares by {0} days".format(append))

    # decrease threshold flux if there's no flare detected
    if decr_thr:
        while (not len(idx[0])) and thr_flux >= 1.5 and thr_flux > min_fl_flux:
            idx = np.where(flux > avg_flux)
            thr_flux -= 0.1
    #print "idx above thr_flux: ", idx, thr_flux # these are the indices of the blocks with flux > thr_flux
    
    # create an interval around each flare: include everything up to first block where flux is 
    # down to some level (e.g. 2 * average)
    # only works if there are more than one bayesian block
    if flux.size == 1:
        logging.warning("Only one bayesian block in time series!")
        return np.array([])

    flare_ids = []
    
    for ij,ii in enumerate(idx[0]):
        flare_blocks = [ii]
        if ij > 0: # check if flare is already contained in last interval
            if ii in flare_ids[-1]:
                continue
        db = 1
        # go back in time from flare
        while flux[ii - db] > avg_flux * min_fl_flux:
                flare_blocks.insert(0,ii - db)
                db += 1
                if ii - db <= 0: 
                    db -= 1
                    break
        
        if ii == flux.size - 1:
            #print "flare blocks: ", ij, flare_blocks
            continue
        # go forward in time from flare
        db = 1
        while flux[ii + db] > avg_flux * min_fl_flux: 
            flare_blocks.append(ii + db)
            db += 1
            if ii + db >= flux.size: 
                db -= 1
                break
            
        flare_blocks.append(ii + db)    
        #print "flare blocks: ", ij, flare_blocks

        flare_ids.append(np.array(flare_blocks))
    
    t_flare = []
    #print 'flare_ids'
    #print flare_ids
    #for ij, fi in enumerate(flare_ids):
        #t_flare.append([tmin[fi].min(), tmin[np.array(fi).max()]])

    for ij, fi in enumerate(flare_ids):

        if len(t_flare) and t_flare[-1][-1] >= tmin[np.array(fi).max()]:
            continue

        tf_min = tmin[fi].min() - append
        tf_max = tmin[np.array(fi).max()] + append

        add = 1
        while ij < len(flare_ids) - add and tf_max >= tmin[flare_ids[ij+add]].min() - append:
            tf_max = tmin[np.array(flare_ids[ij + add]).max()] + append
            add += 1
        t_flare.append([tf_min, tf_max])

    return np.array(t_flare), flare_ids

def calc_bb_unbinned(gta, energies, times, conv, prob,
                        tmin = 0., tmax = 1e20):
    """
    Calculate unbinned bayesian blocks from a 
    FT1 file with src probabilities

    Parameters
    ----------
    gta: `~fermipy.GTAnalysis` object
        The fermipy analysis object
    energies: `~numpy.ndarray`
        array with photon energies
    times: `~numpy.ndarray`
        array with photon arrival times
    conv: `~numpy.ndarray`
        array with photon conversion types
    prob: `~numpy.ndarray`
        array with photon src probabilities

    {options}

    tmin: float
        mininum time for events in MET
    tmax: float
        maxinum time for events in MET
        
    Returns
    -------
    tuple with bins, times, exposure , src probability, and energies
    """
    from fermiAnalysis import adaptivebinning as ab

    # cut on times
    m = (times >= tmin ) & ( times < tmax )
    energies = energies[m]
    times = times[m]
    conv = conv[m]
    prob = prob[m]

    # calculate exposure in bins of energy
    EMeVbins = gta.energies
    EMeVbins = np.array([1e2, 1e3])
    EMeVbins = np.logspace(
                np.log10(energies.min()),
                np.log10(energies.max()), 4)
    EMeV = np.sqrt(EMeVbins[1:] * EMeVbins[:-1])
    EMeV = [500.]

    exp = np.zeros_like(prob)
    for ie, e in enumerate(EMeV):
        me = (energies >= EMeVbins[ie]) & (energies < EMeVbins[ie+1])

        texp, f, b = ab.comp_exposure_phi(gta, energy = e)
        mt = (texp >= times.min() - 60.) & (texp < times.max() + 60.)
        splinef = interp1d(texp[mt], f[mt], kind = 'nearest')
        splineb = interp1d(texp[mt], b[mt], kind = 'nearest')

        mef = me & conv.astype(np.bool) & \
                (times >= texp[mt].min()) & \
                (times <= texp[mt].max())
        meb = me & (~conv.astype(np.bool)) & \
                (times >= texp[mt].min()) & \
                (times <= texp[mt].max())

        exp[meb] = splineb(times[meb]) / b.max()
        exp[mef] = splinef(times[mef]) / f.max()

    if np.sum(exp > 0.):
        logging.info("Calculating unbinned BBs")
        bins = bayesian_blocks(
            times[exp > 0.], exp = (exp / prob)[exp > 0.]
            )
    else:
        logging.error("Exp > 0. everywhere, selecting one bin over entire time range")
        bins = np.array([times.min(), times.max()])
    return bins, times, exp, prob, energies
