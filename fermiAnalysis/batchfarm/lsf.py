"""
Utility functions for sumbmission to an lsf cluster
"""

# --- imports ------------------ #
import yaml
import shlex
import logging
import haloanalysis.batchfarm as bf
from haloanalysis.batchfarm.utils import mkdir
from subprocess import call,check_call,Popen,PIPE
from time import sleep
from os import environ,system,access,R_OK,listdir
from os.path import *
from numpy import string_
from numpy.random import randint
# ------------------------------ #
lsfDefaults = {
    'queue': 'time',
    'time': '05:59',
    'n': 0,
    'span': 'span[ptile=1]', # could also be span[hosts=1]
    'jname': 'lsf',
    'sleep': 10.,
    'lsb_steps': 1,
    'concurrent': 0,
    'dependency': None,
    'minimumJID':1,
    'forceJob':'0',
    'nolog':True,
    'extraDelay':False,
    'lsb_steps':1,
    'logdir':'./log',
    'tmpdir':'./tmp'
}

def setLsf(func):
    """
    Read in default lsf keywords and pass to function
    """
    def init(*args, **kwargs):
	for k in lsfDefaults.keys():
	    kwargs.setdefault(k,lsfDefaults[k])
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


def init_lsf(local_id = 0):
    """
    Init lsf cluster jobs: set up tmpdir on cluster scratch, determine job_id and set pfiles to tmpdir on scratch

    kwargs
    ------
    local_id:	int, if not on lsf, return this value as job_id

    Returns
    -------
    tuple with tmpdir on cluster scratch and lsf job id
    """

    try:
	environ.keys().index("LSB_JOBNAME")
	job_id		= int(environ["LSB_JOBINDEX"])
	tmpdir		= mkdir(join('/scratch/{0:s}.{1:s}/'.format(environ['USER'],environ["LSB_JOBID"])))
	logging.info('os.listdir: {0}'.format(listdir(tmpdir)))

	sleep(20.)
	tmpdir		= join(tmpdir,'{0:s}.XXXXXX'.format(environ["LSB_JOBID"]))
	p		= Popen(["mktemp", "-d", tmpdir], stdout=PIPE)  # make a temporary directory and get its name
	sleep(20.)
	out, err	= p.communicate()

	logging.info('out: {0}'.format(out))
	logging.info('err: {0}'.format(err))

	tmpdir		= join('/scratch',out.split()[0])
	sleep(20.)
    except ValueError:
	job_id 	    	= local_id
	#tmpdir		= environ["PWD"]
	tmpdir		= mkdir(join(environ["PWD"],'tmp/'))

    #system('export PFILES={0:s}:$PFILES'.format(tmpdir))
    #sleep(10.)

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
    kwargs.setdefault('bash_template', join(dirname(bf.__file__),'std_lsf_bash.sh'))

    kwargs.update({'script' : script, 
	    'config' : config , 
	    'logdir' : logdir ,
	    'scriptfile' : basename(script), 
	    'configfile' : basename(config)})

    f = open(kwargs['bash_template'], 'r')

    bash = f.read()
    f.close()
    bash = bash.format(kwargs)

    return bash

@setLsf
def submit_lsf(script,config,option,njobs,**kwargs):
    """
    Create a bash script and submit it to the lsf cluster

    Parameters
    ----------
    script:		string, full path to python script that is executed on lsf cluster
    config:		dict with options to be stored in a yaml file and parsed by python script
    option:		some additional option to be parsed by the python script
    njobs:		either int: number of cluster array jobs; or list with job numbers 

    kwargs
    ------
    queue:		string, queue of lsf cluster (default: long)
    n:			int, number of processor requested for mpi jobs
    ptile:		int, number of processors per hosts for mpi jobs 
    lsb_steps:		int, step width for job indeces, only applies if njobs is int (default: 1)
    time:		string, if queue == time then this determines the cpu running time asked for. 
			Format has to be 'hh:mm'
    jname:		string, name of lsf cluster job (default: lsf)
    sleep:		float, seconds to sleep (default: 10s)
    nolog:		bool, if True (default), send standard output and stderr of cluster job to /dev/null
    			(note that the output from the python script will still be saved)
    log:		string, name of log file (default: tmpdir/jname.out)
    err:		string, name of err file (default: tmpdir/jname.err)
    dry:		bool, if false submit job to lsf cluster (default: False)
    concurrent:		int, limit the number of simultaneously running jobs. If zero (default): no limit.
    dependency:		string, if given, job id of job that needs to have ended before current job is started. 
    			Should be of the form "myjob1", "myjob[1-10]", or "1234", "1234[1-10]", where 1234 is the job id.
    forceJob:		str, if not '0', njobs is over-written with the value of this keyword. 
    extraDelay:		bool, if true, add an extra random delay (between 10 seconds and 3 minutes) 
			to sleep in the bash script before python is called.
    tmpdir:		string, local directory for temporary storage of bash script
    logdir:		string, local directory for log files
    """
    mkdir(kwargs['logdir'])
    mkdir(kwargs['tmpdir'])
    yamlfile = join(kwargs['tmpdir'],'{0[jname]:s}_{1[configname]:s}.yaml'.format(kwargs,config))
    yaml.dump(config, stream = open(yamlfile, 'w'), default_flow_style=False)
    # test yaml file
    par	= yaml.load(open(yamlfile))

    bash	= make_bash(script,yamlfile,kwargs['logdir'], 
				add_opt = option, 
				sleep = kwargs['sleep'], 
				extrasleep = kwargs['extraDelay']
				)
    bashScript	= join(kwargs['tmpdir'],'{0[jname]:s}_{1[configname]:s}.sh'.format(kwargs,config))

    f = open(bashScript,'w')
    f.write(bash)
    f.close()

    call(['chmod','u+x',bashScript])

    if kwargs['forceJob'] == '0':
	if type(njobs) == int:
	    if not kwargs['minimumJID']: kwargs['minimumJID'] = 1
	    njobs = '[{0[minimumJID]:n}-{1:n}:{0[lsb_steps]:n}]'.format(kwargs,njobs)
	    nsubmit = 1
	elif type(njobs) == list:
	    if len(njobs) > 100:
		nsubmit = len(njobs) / 100 if len(njobs) % 100 == 0 else len(njobs) / 100 + 1
	    else:
		nsubmit = 1
    else:
	njobs	= kwargs['forceJob']
	nsubmit	= 1

    for i in range(nsubmit):
	if nsubmit == 1:
	    if not kwargs['concurrent']:
		command	= """bsub -oo {0[log]:s} -eo {0[err]:s} -J "{0[jname]:s}{1}" """.format(kwargs,njobs)
	    else:
		command	= """bsub -oo {0[log]:s} -eo {0[err]:s} -J "{0[jname]:s}{1}%{0[concurrent]:n}" """.format(kwargs,njobs)
	else:
	    if not kwargs['concurrent']:
		command	= """bsub -oo {0[log]:s} -eo {0[err]:s} -J "{0[jname]:s}{1}" """.format(kwargs,
		    njobs[i * 100:(i + 1) * 100 if (i + 1) * 100 < len(njobs) else -1]
		    )
	    else:
		command	= """bsub -oo {0[log]:s} -eo {0[err]:s} -J "{0[jname]:s}{1}%{0[concurrent]:n}" """.format(kwargs,
		    njobs[i * 100:(i + 1) * 100 if (i + 1) * 100 < len(njobs) else -1]
		    )
	if kwargs['queue'] == 'time':
	    command += """-W {0[time]:s} """.format(kwargs)
	else:
	    command += """-q {0[queue]:s} """.format(kwargs)

	if kwargs['n'] > 0:
	    command += """-n {0[n]:n} -R "{0[span]:s}" """.format(kwargs)

	if not kwargs['dependency'] == None:
	    command += """-w "ended({0[dependency]:s})" """.format(kwargs)

	command += """ {0:s} """.format(bashScript)

	if not kwargs['dry']:
	    logging.info('Sending command\n\t{0:s}\nto lsf cluster'.format(command))
	    call(shlex.split(command))
	else:
	    logging.info('Dry run for command\n\t{0:s}\nto lsf cluster'.format(command))

	logging.info('Going to sleep for {0[sleep]:.2f} s ...'.format(kwargs))
	sleep(kwargs['sleep'])

    return 
