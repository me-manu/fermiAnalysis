"""
utility functions 
"""

# --- Imports -------------- #
import sys
from glob import glob
from subprocess import call,check_call
from os.path import *
from os import access,R_OK
from numpy import isscalar,array,where,sqrt,where,log10,linspace,cumsum,argmin,abs,empty,pi
import numpy as np
from time import sleep
from numpy import string_
from subprocess import check_call
import logging
import pickle 
import gzip
# -------------------------- #

def save(object, filename, protocol = -1):
    """
    Save an object to a compressed disk file.
    Works well with huge objects.

    Parameters
    ----------
    object:        some python object
    filename:        str, full path to file

    kwargs
    ------
    protocol:        int, see Pickler documentation
    """
    file = gzip.GzipFile(filename, 'wb')
    pickle.dump(object, file, protocol)
    file.close()

def load(filename):
    """
    Loads a compressed object from disk

    Parameters
    ----------
    filename:        str, full path to file

    Returns
    -------
    The loaded object.
    """
    file = gzip.GzipFile(filename, 'rb')
    object = pickle.load(file)
    file.close()

    return object

def set_njobs(njobs, missing, irf = "None", clobber = False):
    """
    determine number of jobs

    Parameters
    ----------
    njobs:        int, total number of jobs
    missing:        list, list with job numbers

    kwargs
    ------
    irf:        str, irf name
    clobber: bool, if true, overwrite existing files

    Returns
    -------
    list or int with job ids to be submitted to lsf cluster
    """
    if len(missing) < njobs:
        if not len(missing):
            logging.debug('all files present for irf {0:s}.'.format(irf))
            if clobber:
                logging.warning('Clobber is set to yes, will overwrite files')
            else:
                njobs = 0
        else:
            njobs = missing
            logging.info('there are {0:n} files missing'.format(len(missing)))

    logging.debug('missing: {0}'.format(missing))
    return njobs

def mkdir(directory, log = True, dry = False):
    """
    Create a new directory

    Paramters
    ---------
    directory: string of new directory

    kwargs
    ------
    log:        bool, if True, report new directory 
    dry:        bool, if True, don't create dir

    Returns
    -------
    directory: string of new directory
    """
    if not exists(directory):
        if not dry:
            call(['mkdir','-p',directory])
            sleep(1.)
        if log:
            logging.info('Set up directory {0:s}'.format(directory))

    return directory

def rm(path, dry = False):
    """
    Remove a file or a directory

    Parameters
    ----------
    path:        str, full path to file or directory that should be deleted

    kwargs
    ------
    dry:        bool, if true, do not remove but pring messages

    Returns
    -------
    int, if zero, deletion successfull.
    """

    if type(path) == str or type(path) == string_:
        if not exists(path):
            logging.error("Path {0:s} does not exist".format(path))
            return -1
        if isdir(path):
            logging.warning("*** Deleting {0:s}".format(path))
            sleep(0.5)
            if not dry:
                call(['rm','-r',path])
        elif isfile(path):
            logging.warning("*** Deleting {0:s}".format(path))
            sleep(0.5)
            if not dry:
                call(['rm',path])
        else:
            logging.error("Path {0:s} does not point to a file or directory")
            return -2

    elif type(path) == list:
        logging.warning("*** Deleting {0:s}".format(path))
        sleep(0.5)
        if isdir(path[0]):
            if not dry:
                call(['rm','-r'] + path)
        elif isfile(path[0]):
            if not dry:
                call(['rm'] + path)
    return 0

def update_diff(dict1,dict2):
    """
    Update dict1 with values of dict2 without adding additional keys from dict2 to dict1

    Parameters
    ----------
    dict1:        old dictionary
    dict2:        dictionary with updated values and possibly additional keys

    Returns
    -------
    updated dict1 and dict2 with updated items removed
    """
    diff_keys1        = set(dict2) - set(dict1)
    diff_keys2        = set(dict1) - set(dict2)
    dict1.update(dict2)

    for k in diff_key1s: del dict1[k]        # remove unwanted / additional keywords from dict1
    for k in diff_key2s: del dict2[k]        # remove unwanted / additional keywords from dict1

    return dict1,dict2


def missing_files(fname,fn, missing = True, minimum = 0, num = 5,
                    startnum = 1, split = '.dat.gz', folder = False):
    """
    Check for missing files.

    Parameters
    ----------
    fname:        string, file names (with wildcards), last num spaces need to be file numbers, e.g., if num = 5, file00001.file
    fn:                int, number of files that should be present

    kwargs
    ------
    minimum:        int, minimum job id that will be resubmitted (default: 0)
    missing:        bool, if True, look for missing files, if false, look for present files.
    num:        int, number of digits of file identifier (default: 5)
    split:        str, str used to split file names, default: '.dat.gz'
    folder:         bool, if true, numerals with job numbers are given in directory as /some/path/00001/output.file
    startnum:   int, starting number for file ids (default: 1)

    Returns
    -------
    list with missing file numbers
    """
    logging.debug(fname)
    files = glob(fname)
    
    if folder:
        if split == '':
            idxs  = array(map( lambda f: int(dirname(f).split('/')[-1][-num:]), files))
        else:
            idxs  = array(map( lambda f: int(dirname(f).split('/')[-1].split(split)[0][-num:]), files))
    else:
        idxs  = array(map( lambda f: int(basename(f).split(split)[0][-num:]), files))

    miss = []
    logging.debug('number of files that should be there: {0:n}'.format(fn))
    logging.debug('file ids: {0}'.format(idxs))
    if not fn == len(files):
        for idxi in range(startnum,fn + startnum):
            if missing:
                if not len(where(idxi == idxs)[0]):
                    miss.append(idxi)
            else:
                if len(where(idxi == idxs)[0]):
                    miss.append(idxi)
    if fn == len(files) and not missing:
        miss = range(startnum,fn + startnum)

    logging.debug('missing: {0}{1}, min: {2}'.format(missing,miss,minimum))

    if minimum:
        miss = np.array(miss)
        miss = miss[miss > minimum]
        logging.debug('2nd missing: {0}{1}, min: {2}'.format(missing,miss,minimum))
        return list(miss)

    else:
        return miss
def parse_logfiles(path2log, string, nostring = None, filenames = 'err.*'):
    """
    Parse log files of job array and return list with job numbers
    where string was found in log file. 
    It's assumed that the job ids are giving by file endings
    of log files

    Parameters
    ----------
    path2log: str
        full path to logfiles
    string: str
        string to be searched for in log file

    kwargs
    ------
    filenames: str
        wild card for log files, default is "err.*"
    nostring: str or none
    if given, do not exlude logfile if line contains this string
    """
    logfiles = glob(join(path2log,filenames))
    logfiles = sorted(logfiles, key = lambda f: int(basename(f).split('.')[-1]))
    ids = array([int(basename(f).split('.')[-1]) for f in logfiles])
    logging.info("Found {0:n} log files".format(len(ids)))
    # list of indices where string was found
    idx = []
    # loop through files
    for i,lf in enumerate(logfiles):
        for line in open(lf):
            if string in line:
                if not type(nostring) == type(None):
                    if not nostring in line:
                        idx.append(i)
                else:
                    idx.append(i)
    logging.info("Found {0:n} log files that contain {1:s} which will be removed from missing files".format(len(idx),string))
    return ids[idx]




def filesExists(*files):
    """
    test if file(s) exists.

    Parameters
    ----------
    *files: list of file names

    Returns
    -------
    bool, 0 if one file is missing, 1 if all files are present.
    """
    res = 0
    for f in files:
        if exists(f):
            logging.info('Found file {0}'.format(f))
            res += 1
    if res == len(files): 
        return True
    else: 
        return False

def unzipfiles(filename, outdir = None):
    """
    unzip a file

    Paramters
    ---------
    filename:        string, the file to unzip

    kwargs
    ------
    outdir:        str or None (default), path where files are unzipped in. If None, unzip in current directory
    """
    if outdir == None:
        call(['tar']+['-xzvf']+[filename])
    else:
        call('tar -xzvf {0:s} -C {1:s}'.format(filename,outdir),shell = True)
    return

def zipfiles(filename,output,nodir = False):
    """
    zip a file or a list of files

    Paramters
    ---------
    filename:        string, list, or dict of the file(s) to zip
    output:        string, output file name,R_OK

    kwargs
    ------
    nodir:        bool, if True, do not retain directory structure

    Returns
    -------
    output
    """
    if type(filename) == list:
        if nodir:
            for i,f in enumerate(filename):
                if not i:
                    call(['tar']+['-cf']+[output]+['--directory={0:s}'.format(dirname(f))]+[basename(f)])
                else:
                    call(['tar']+['--append']+['--file={0:s}'.format(output)]+['--directory={0:s}'.format(dirname(f))]+[basename(f)])
            call(['gzip']+['-f']+[output])
        else:
            call(['tar']+['-czf']+[output]+filename)
    elif type(filename) == dict:
        if nodir:
            for i,f in enumerate(filename.values()):
                if not i:
                    call(['tar']+['-cf']+[output]+['--directory={0:s}'.format(dirname(f))]+[basename(f)])
                else:
                    call(['tar']+['--append']+['--file={0:s}'.format(output)]+['--directory={0:s}'.format(dirname(f))]+[basename(f)])
            call(['gzip']+['-f']+[output])
        else:
            call(['tar']+['-czf']+[output]+filename.values())
    elif type(filename) == str or type(filename) == string_:
        if nodir:
            call(['tar']+['-czf']+[output]+['--directory={0:s}'.format(dirname(filename))] + [basename(filename)])
        else:
            call(['tar']+['-czf']+[output]+[filename])
    return output

def copy2scratch(filename,dir,zipping = False, dry = False, stime = 10.):
    """
    Copy a file or a list of files 

    Parameters
    ----------
    filename:        string, list, or dict of the file(s) to copy
    tmpdir:        string, target to directory for file

    kwargs
    ------
    zipping:        bool, if True, zip files before copying and unzip at target location
    dry:        bool, if True, only set path name, do not copy (default = False)

    Return
    ------
    if filename is string: returns new path to filename, else returns just dir
    """
    if not check_read_access(filename):
        #raise ValueError('*** no read access to all files requested for copying: {0}'.format(filename))
        logging.error('*** no read access to all files requested for copying: {0}'.format(filename))

    if type(filename) == list:
        if not dry:
            if zipping:
                zipfiles(filename,'tmp.tar.gz')
                check_call(['cp'] + ['tmp.tar.gz'] + [dir])
                check_call(['tar']+['-xzf']+[join(dir,'tmp.tar.gz')])
                check_call(['rm']+['tmp.tar.gz'])
            else:
                check_call(['cp'] + filename + [dir])
            logging.info('copied {0:s} to {1:s}'.format(filename,dir))
            sleep(stime)
        return dir
    elif type(filename) == dict:
        if not dry:
            check_call(['cp'] + filename.values() + [dir])
            logging.info('copied {0:s} to {1:s}'.format(filename,dir))
            sleep(stime)
        return dir
    elif type(filename) == str or type(filename) == string_:
        if not dry:
            check_call(['cp',filename,dir])
            logging.info('copied {0:s} to {1:s}'.format(filename,dir))
            sleep(stime)
        return join(dir,basename(filename))
    else:
        raise TypeError('*** filename must be of type str or list, is {0}'.format(type(filename)))
        return

def check_read_access(path):
    """
    Check if user has read access to one or multiple files

    Parameters
    ----------
    path:        str, list, or dict: path to file(s)

    Returns
    -------
    bool:        True if read access (to all files if path is list or dict), else False
    """
    if type(path) == list:
        for p in path:
            if not access(p, R_OK): 
                logging.warning('*** no read access to {0:s}'.format(p))
                return False
    elif type(path) == dict:
        for p in path.values():
            if not access(p, R_OK): 
                logging.warning('*** no read access to {0:s}'.format(p))
                return False
    elif type(path) == str or type(path) == string_:
        if not access(path, R_OK): 
                logging.warning('*** no read access to {0:s}'.format(path))
                return False
    return True

def init_logging(level, color = False):
    """
    Setup logger.

    Parameters 
    ----------
    level:        string, level of logging: DEBUG,INFO,WARNING,ERROR. (default: INFO).

    kwargs
    ------
    color:        bool, if true, enable colored output for bash output

    Notes
    -----
    for color see
        stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
        https://wiki.archlinux.org/index.php/Color_Bash_Prompt
    """
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if level.upper() == 'INFO':
        level = logging.INFO
    elif level.upper() == 'DEBUG':
        level = logging.DEBUG
    elif level.upper() == 'WARNING':
        level = logging.WARNING
    elif level.upper() == 'ERROR':
        level = logging.ERROR


    if color:
        logging.basicConfig(level=level,stream = sys.stderr, format='\033[0;36m%(filename)10s:\033[0;35m%(lineno)4s\033[0;0m --- %(levelname)7s: %(message)s')
        logging.addLevelName( logging.DEBUG, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
        logging.addLevelName( logging.INFO, "\033[1;36m%s\033[1;0m" % logging.getLevelName(logging.INFO))
        logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
        logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))
    else:
        logging.basicConfig(level=level,stream = sys.stderr, format='%(filename)10s:%(lineno)4s --- %(levelname)7s: %(message)s')

    return

def printOneLine(string):
    """
    Print an output on the same line

    Parameter
    ---------
    string:        str,the string to be printed
    """
    sys.stdout.write(string)
    sys.stdout.flush()
    lenStr = len(string)
    sys.stdout.write("\b"*lenStr)
    return
