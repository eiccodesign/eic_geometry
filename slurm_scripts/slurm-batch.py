#! /usr/bin/env python2.7
# NOTE: THE COMMANDS NEED TO BE MODIFIED FOR YOUR USE CASE
from __future__ import print_function
import os
import shutil
import argparse
import subprocess

# -------------------------- ARGUMENT PARSING ------------------------ #

prog_description = 'Creates files needed for slum jobs and optionally submits jobs'
prog_epilogue    = 'And that\'s all, folks!'

# Initialize argument parser
parser = argparse.ArgumentParser(
    description=prog_description,
    usage='%(prog)s [options]',
    epilog=prog_epilogue,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# Arguments 
parser.add_argument('-d', '--dir',type=str,
                    default='/p/lustre2/milton3/epic_zdc_lambda/',
                    help='top level scratch directory')
parser.add_argument('-j', '--jobset',type=str,default='zdc_neutron_log10continuous_01_25_2024',
                    help='an identifier for the jobset, jobs submitted from a subdir of this name')
parser.add_argument('-nj', '--njob', type=int, default='10',
                    help='number of jobs to submit')
parser.add_argument('-s', '--submit',action='store_true',
                    help='submit jobs with sbatch after creating')
parser.add_argument('-m', '--merge',action='store_true',
                    help='submit separate merge job with dependencies on this one')
parser.add_argument('-v', '--verbose',action='store_true',
                    help='turn on diagnostic print statements')
parser.add_argument('-n', '--nevents',type=str,default='1000',
                    help='number of events to run per subjob')
parser.add_argument('-part', '--particle',type=str,default='pi+',
                    help='particle species: pi+, pi-, pi0, mu+, mu-, etc')
parser.add_argument('-pmin', '--energy_min',type=str,default='10.0',
                    help='min generated momentum')
parser.add_argument('-pmax', '--energy_max',type=str,default='300.0',
                    help='max generated momentum')

# Parse argument
args = parser.parse_args()

if not os.path.isdir(args.dir):
    raise ValueError('directory {} not found.'.format(args.dir))

#parameters for this example
bank='qcdtq'
subjob_script='run_gun.sh'
merge_script='run-merge.sh'
time_est_subjob='10:00:00' # Ten hour
time_est_merge='00:01:00' #two hours

#create the submission directory args.dir/args.jobset
submit_dir = os.path.join(args.dir,args.jobset)
if (args.verbose):
    print('Submission directory: ', submit_dir)

if not os.path.isdir(submit_dir):
    os.mkdir(submit_dir) 
    os.chmod(submit_dir,0o2770)

#also create subdirectory for subjob log files
logdir = os.path.join(submit_dir,'log')
if not os.path.isdir(logdir):
    os.mkdir(logdir) 
    os.chmod(logdir,0o2770) 

#...and output
outputdir = os.path.join(submit_dir,'output')
if not os.path.isdir(outputdir):
    os.mkdir(outputdir) 
    os.chmod(outputdir,0o2770) 

#copy scripts to submission directory
#could alternatively put them in your path
extra_files=[subjob_script,merge_script]
for script in extra_files:
    script_src = os.path.join("./",script)
    script_dest = os.path.join(submit_dir,script)
    shutil.copyfile(script_src,script_dest)
    shutil.copystat(script_src,script_dest)

#Create a job script for slurm
basefile = os.path.join(submit_dir,'job')
batchfile = basefile + '.sh'
f = open(batchfile,'w')
f.write('#!/bin/bash\n')
f.write('#SBATCH -n 1\n')
f.write('#SBATCH -t %s\n' % time_est_subjob)
f.write('#SBATCH --job-name=%s\n' % args.jobset)
f.write('#SBATCH -p pbatch\n')
# f.write('#SBATCH -p pdebug\n')
f.write('#SBATCH -A %s\n' % bank)
f.write('#SBATCH --array=0-%d\n' % (args.njob-1))
#the %A's are slurm syntax to insert the numeric job ID and subjob index
f.write('#SBATCH -o %s%s.%s.out\n' % (os.path.join(submit_dir,'log/'),'%A','%a')) 
f.write('#SBATCH -e %s%s.%s.err\n' % (os.path.join(submit_dir,'log/'),'%A','%a'))
f.write('%s/%s -n %s -d %s -part %s -j %s -pmin %s -pmax %s' %(submit_dir,subjob_script,
                                                     args.nevents, submit_dir, args.particle, args.jobset, args.energy_min, args.energy_max) )
# IMPORTANT FOR EDITORS: ^This line parses all the arguments, and interfaces to the run_gun.sh script
# The flags here are similar between slurm-batch.py and run_gun.sh, but there are important differences.

f.close()

#submit job
if (args.submit):
    submit_command = 'sbatch ' + batchfile
    print (submit_command)
    os.system(submit_command)
else:
    print('%s/%s -n %s -d %s -part %s -j %s -pmin  %s -pmax %s' %(submit_dir,subjob_script,
                                                                 args.nevents, submit_dir, args.particle, args.jobset,args.energy_min, args.energy_max) )

    os.system('%s/%s -n %s -d %s -part %s -j %s -pmin  %s -pmax %s' %(submit_dir,subjob_script,
                                                                     args.nevents, submit_dir, args.particle, args.jobset, args.energy_min, args.energy_max) )
#Now optionally create job script for merge job and possibly submit
if ( args.merge ):
    batchfile = basefile + '.merge.sh'
    f = open(batchfile,'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -n 1\n')
    f.write('#SBATCH -t %s\n' % time_est_merge)
    f.write('#SBATCH --job-name=%s\n' % args.jobset)
    f.write('#SBATCH -p pbatch\n')
    f.write('#SBATCH -A %s\n' % bank)
    f.write('#SBATCH -o %smerge.out\n' % os.path.join(submit_dir,'log/'))
    f.write('#SBATCH -e %smerge.err\n' % os.path.join(submit_dir,'log/'))
    f.write("#SBATCH --dependency=singleton\n")
    f.write('%s/%s %s' % (submit_dir,merge_script,submit_dir))
    f.close()
    if (args.submit):
        submit_command = 'sbatch ' + batchfile
        print (submit_command)
        os.system(submit_command)
