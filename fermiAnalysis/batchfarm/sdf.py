"""
Utility functions for sumbmission to an lsf cluster
"""

# --- imports ------------------ #
#import yaml
import shlex
import logging
import fermiAnalysis.batchfarm as bf
import re
import sys
from fermiAnalysis.batchfarm.utils import mkdir
from subprocess import call,check_call,Popen,PIPE,check_output,CalledProcessError
from time import sleep
from os import environ,system,access,R_OK,listdir
from os.path import *
from numpy import string_
from numpy.random import randint
from astropy.io.misc import yaml
from datetime import timedelta
# ------------------------------ #
sdf_defaults = {
    'queue': 'time',
    'time': 1439,
    'jname': 'sdf',
    'sleep': 10.,
    'lsb_steps': 1,
    'concurrent': 0,
    'dependency': None,
    'minimumJID': 1,
    'forceJob': '0',
    'nolog': True,
    'extraDelay': False,
    'lsb_steps': 1,
    'logdir': './log',
    'tmpdir': './tmp',
    'max_rjobs': None,
    'nminjob': None, 
    'nmaxjob': None,
    'ntasks_per_node': 1,
    'nodes': 1,
    'partition': 'milano',
    'no_resubmit_running_jobs': True
}

def set_sdf(func):
    """
    Read in default lsf keywords and pass to function
    """
    def init(*args, **kwargs):
        for k in sdf_defaults.keys():
            kwargs.setdefault(k, sdf_defaults[k])
        if kwargs['nolog']:
            kwargs.setdefault('log','/dev/null')
            kwargs.setdefault('err','/dev/null')
        else:
            kwargs.setdefault('log',join(kwargs['logdir'],'{0:s}.out'.format(kwargs['jname'])))
            kwargs.setdefault('err',join(kwargs['logdir'],'{0:s}.err'.format(kwargs['jname'])))
        kwargs.setdefault('dry',False)
        if kwargs['extraDelay']:
            kwargs['extraDelay'] = randint(10,3*60)
        else:
            kwargs['extraDelay'] = 0.
        
        return func(*args, **kwargs)
    return init

def get_jobs(user="mmeyer"):
    p = Popen(['sacct'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    try:
        out = check_output(("grep", user), stdin=p.stdout)
        return len(out.split('\n')) - 1
    except CalledProcessError:
        return 0

def get_running_jobs_list(user="mmeyer"):
    """retrun output list of running jobs with full name"""
    command = "squeue -u {0:s} -r --Format JobArrayID,NAME:100".format(user)
    p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    lines = p.stdout.readlines()
    if sys.version_info[0] >= 3:
        lines = [line.decode('utf-8') for line in lines]
    return lines

def remove_running_job_from_list(job_name, job_list, user="mmeyer"):
    """Remove a job from a job if it is currently running"""
    out = get_running_jobs_list(user=user) 

    # get the job id and array id
    x = [job.split()[0] for job in out if job_name in job] 

    # from this array id, which is separated from job id by underscore
    array_ids = [int(j.split("_")[1]) for j in x]

    # remove the array_ids of the running jobs from the job list
    # and return the new job list
    for i in array_ids:
        try:
            job_list.remove(i)
        except ValueError:  # job is not in list 
            pass

    return job_list


def init_sdf(local_id=0):
    """
    Init sdf cluster jobs: set up tmpdir on cluster scratch, determine job_id and set pfiles to tmpdir on scratch

    kwargs
    ------
    local_id: int
        if not on lsf, return this value as job_id

    Returns
    -------
    tuple with tmpdir on cluster scratch and lsf job id
    """

    try:
        list(environ.keys()).index("LSB_JOBNAME")
        job_id = int(environ["SLURM_ARRAY_TASK_ID"])
        tmpdir = environ["LSCRATCH"]
        logging.info('os.listdir: {0}'.format(listdir(tmpdir)))

        sleep(20.)
        tmpdir = join(tmpdir,'{0:s}.XXXXXX'.format(environ["SLURM_JOB_ID"]))
        p = Popen(["mktemp", "-d", tmpdir], stdout=PIPE)  # make a temporary directory and get its name
        sleep(20.)
        out, err = p.communicate()

        logging.info('out: {0}'.format(out))
        logging.info('err: {0}'.format(err))

        tmpdir = join(environ["LSCRATCH"], out.decode('ascii').split()[0])
        sleep(20.)
    except ValueError:
        job_id = local_id
        tmpdir = mkdir(join(environ["PWD"],'tmp/'))

    logging.info('tmpdir is {0:s}.'.format(tmpdir))

    if not exists(tmpdir):
        logging.error('Tmpdir does not exist: {0}. Exit 14'.format(tmpdir))
        sys.exit(14)

    return tmpdir,job_id


def make_bash(script,config,logdir,**kwargs):
    """
    return string with bash script

    Parameters
    ----------
    strlist: list with strings, containing 0: python script, 1: yaml file 2: dir for logging 

    kwargs
    ------
    bash_template: shell script template. 
    add_opt: string with additional options passed to python script in strlist[0]
    sleep: seconds to sleep between commands (default: 10)

    Returns
    -------
    bash script string
    """
    kwargs.setdefault('add_opt', '')
    kwargs.setdefault('sleep', 10.)
    kwargs.setdefault('extrasleep', kwargs['sleep'] + 0.)
    kwargs.setdefault('bash_template', join(dirname(bf.__file__),'std_sdf_bash.sh'))

    kwargs.update({'script' : script, 
            'config' : config, 
            'logdir' : logdir,
            'conda_env' : environ['CONDA_DEFAULT_ENV'],
            'scriptfile' : basename(script), 
            'configfile' : basename(config)})

    f = open(kwargs['bash_template'], 'r')

    bash = f.read()
    f.close()
    bash = bash.format(kwargs)

    return bash

@set_sdf
def submit_sdf(script, config, option, njobs, **kwargs):
    """
    Create a bash script and submit it to the sdf cluster

    Parameters
    ----------
    script: string, full path to python script that is executed on lsf cluster
    config: dict with options to be stored in a yaml file and parsed by python script
    option: some additional option to be parsed by the python script
    njobs:  either int: number of cluster array jobs; or list with job numbers 

    kwargs
    ------
    queue:  string, queue of lsf cluster (default: long)
    n: int, number of processor requested for mpi jobs
    ptile: int, number of processors per hosts for mpi jobs 
    lsb_steps: int, step width for job indeces, only applies if njobs is int (default: 1)
    time: string, if queue == time then this determines the cpu running time asked for. 
          Format has to be 'hh:mm'
    jname: string, name of lsf cluster job (default: lsf)
    sleep: float, seconds to sleep (default: 10s)
    nolog: bool, if True (default), send standard output and stderr of cluster job to /dev/null
           (note that the output from the python script will still be saved)
    log: string, name of log file (default: tmpdir/jname.out)
    err: string, name of err file (default: tmpdir/jname.err)
    dry: bool, if false submit job to lsf cluster (default: False)
    concurrent: int, limit the number of simultaneously running jobs. If zero (default): no limit.
    dependency: string, if given, job id of job that needs to have ended before current job is started. 
                Should be of the form "myjob1", "myjob[1-10]", or "1234", "1234[1-10]", where 1234 is the job id.
    forceJob: str, if not '0', njobs is over-written with the value of this keyword. 
    extraDelay: bool, if true, add an extra random delay (between 10 seconds and 3 minutes) 
                to sleep in the bash script before python is called.
    max_rjobs: int or None
        maximum jobs that are allowed to run. If exceeded, wait 20s and try again.
    tmpdir: string, local directory for temporary storage of bash script
    logdir: string, local directory for log files
    no_resubmit_running_jobs: bool
        if job array job is running, don't resubmit
    """
    mkdir(kwargs['logdir'])
    mkdir(kwargs['tmpdir'])
    yamlfile = join(kwargs['tmpdir'],'{0[jname]:s}_{1[configname]:s}.yaml'.format(kwargs,config))
    yaml.dump(config, stream = open(yamlfile, 'w'), default_flow_style=False)
    # test yaml file
    par = yaml.load(open(yamlfile))

    bash = make_bash(script, yamlfile, kwargs['logdir'], 
                     add_opt = option, 
                     sleep = kwargs['sleep'], 
                     extrasleep = kwargs['extraDelay']
                     )
    bashScript = join(kwargs['tmpdir'],'{0[jname]:s}_{1[configname]:s}.sh'.format(kwargs,config))

    f = open(bashScript,'w')
    f.write(bash)
    f.close()

    call(['chmod','u+x',bashScript])

    if kwargs['forceJob'] == '0':
        if kwargs['no_resubmit_running_jobs'] and (isinstance(njobs, list) or isinstance(njobs, int)):
            # check if jobs with same name are already running and remove them
            njobs = remove_running_job_from_list(kwargs['jname'], njobs if isinstance(njobs, list) else list(range(1,njobs + 1,1)))
            if not len(njobs):
                logging.warning("all jobs requested for submission are currently running / pending! Returning without submission")
                return

        if isinstance(njobs, int):
            if not kwargs['minimumJID']: kwargs['minimumJID'] = 1
            njobs = '{0[minimumJID]:d}-{1:d}:{0[lsb_steps]:d}'.format(kwargs,njobs)
            nsubmit = 1
        elif isinstance(njobs, list):
            if len(njobs) > 100:
                nsubmit = len(njobs) / 100 if len(njobs) % 100 == 0 else len(njobs) / 100 + 1
            else:
                nsubmit = 1
    else:
        njobs = kwargs['forceJob']
        nsubmit = 1

    nsubmit = int(nsubmit)

    for i in range(nsubmit):
        command = """sbatch -o {0[log]:s} -e {0[err]:s} -J "{0[jname]:s}" """.format(kwargs)

        if nsubmit == 1:
            if not njobs == "1-1":
                if isinstance(njobs, list):
                    njobs_str = ','.join(map(str,njobs))
                else:
                    njobs_str = njobs
                if not kwargs['concurrent']:
                    command += """ -a {0} """.format(njobs_str)
                else:
                    command += """ -a {1}%{0[concurrent]:d} """.format(kwargs, njobs_str)
        else:
            if not kwargs['concurrent']:
                command += """-a {0:s} """.format(','.join(map(str,njobs[i * 100:(i + 1) * 100
                                                   if (i + 1) * 100 < len(njobs) else -1]))
                                                   )
            else:
                command += """-a {0:s}%{1[concurrent]:d} """.format(','.join(map(str,njobs[i * 100:(i + 1) * 100
                                                                    if (i + 1) * 100 < len(njobs) else -1])),
                                                                    kwargs
                                                                    )

        if kwargs['queue'] == 'time':
            try:
                time = int(kwargs['time'])
            except ValueError:
                delta = timedelta(hours=int(kwargs['time'].split(":")[0]),
                                  minutes=int(kwargs['time'].split(":")[1]))
                time = int(delta.total_seconds() / 60)

            command += """--time {0:d} """.format(time)
        else:
            command += """-q {0[queue]:s} """.format(kwargs)

        if kwargs['ntasks_per_node'] > 1 or kwargs['nodes'] > 1:
            command += """--ntasks-per-node {0[ntasks_per_node]:d} --nodes {0[nodes]:d} """.format(kwargs)

        if not kwargs['dependency'] == None:
            command += """--depend={0[dependency]:s} """.format(kwargs)

        if "partition" in kwargs.keys():
            command += """-p {0[partition]:s} """.format(kwargs)

        if "mem" in kwargs.keys():
            command += """--mem {0[mem]:d} """.format(kwargs)

        if not 'login' in environ['HOSTNAME']:
            command += """ --account fermi:default"""

        # check if we're submitting 
        # from S3DF, if so, add the account option 
        # see https://confluence.slac.stanford.edu/display/SAS/Get+Started+on+S3DF+-+Cheat+Sheet
        command += """ {0:s} """.format(bashScript)

        # get the current number of running jobs
        # and wait 

        if type(kwargs['max_rjobs']) == int:
            rjobs = get_jobs()
            while rjobs >= kwargs['max_rjobs']:
                logging.info('{0:d} jobs running, max number of running jobs allowed: {1:d}'.format(rjobs, kwargs['max_rjobs']))
                logging.info('Sleep for {0:.2f} s ...'.format(kwargs['sleep'] * 3.))
                sleep(kwargs['sleep'] * 3.)
                rjobs = get_jobs()


        if not kwargs['dry']:
            logging.info('Sending command\n\t{0:s}\nto sdf cluster'.format(command))
            call(shlex.split(command))
        else:
            logging.info('Dry run for command\n\t{0:s}\nto sdf cluster'.format(command))

        logging.info('Going to sleep for {0[sleep]:.2f} s ...'.format(kwargs))
        sleep(kwargs['sleep'])

    return 
