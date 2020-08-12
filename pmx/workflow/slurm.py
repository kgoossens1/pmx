#!/usr/bin/env python
# coding: utf-8

"""slurm.py: Handles SLURM queue system commands."""

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def scriptHeader( fname, jobname, simtime=1, simcpu=4):
    fp = open(fname,'w')

    jobline = f'#! /usr/bin/bash\n\
#SBATCH --job-name={jobname}\n\
#SBATCH --error={jobname}.e\n\
#SBATCH --partition=nes2.8\n\n\
#SBATCH --nodes=1\n\
#SBATCH --constraint=wes2.8\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --cpus-per-task={simcpu}\n\
#SBATCH --time={simtime}-00:00:00\n\
\n\
#Function to call to run the actual code\n\
slurm_startjob(){{\n\
#----------------- Actual calculation command goes here: ---------------------------\n\
\n\
module purge\n\
module load gnu/8.2.0  \n\
module load openmpi/3.1.2\n\
\n\
source /export/home/dfhahn/bin/gromacs_2019.4/bin/GMXRC\n\
PATH=/export/home/dfhahn/bin/gromacs_2019.4/bin/:$PATH\n\
\n\
WORKDIR=$(pwd)\n\
mkdir -p $SLURM_TMPDIR/$SLURM_JOB_NAME\n\
cd $SLURM_TMPDIR/$SLURM_JOB_NAME\n\n'
    fp.write(jobline)
    fp.close()


def scriptAppend( fname, commands ):
    fp = open(fname, 'a')
    fp.write(commands)
    fp.close()

def scriptMain( fname, gromppline, simpath, runtype, run, simtime=4, simcpu=1 ):
    fp = open(fname, 'a')

    jobline = f'{gromppline} \n\n\
gmx mdrun -ntomp {simcpu} -ntmpi 1 -s {runtype}{run}.tpr -deffnm {runtype}{run}\n\
rsync -avzIi ./ {simpath}/\n\n'

    fp.write(jobline)
    fp.close()


def scriptFooter( fname ):
    fp = open(fname, 'a')

    jobline = f'\n\
}}\n\
\n\
slurm_info_out(){{\n\
\n\
echo "=================================== SLURM JOB ==================================="\n\
date\n\
echo\n\
echo "The job will be started on the following node(s):"\n\
echo $SLURM_JOB_NODELIST\n\
echo\n\
echo "Slurm User:         $SLURM_JOB_USER"\n\
echo "Run Directory:      $(pwd)"\n\
echo "Job ID:             ${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"\n\
echo "Job Name:           $SLURM_JOB_NAME"\n\
echo "Partition:          $SLURM_JOB_PARTITION"\n\
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"\n\
echo "Number of tasks:    $SLURM_NTASKS"\n\
echo "Submitted From:     $SLURM_SUBMIT_HOST:$SLURM_SUBMIT_DIR"\n\
echo "=================================== SLURM JOB ==================================="\n\
echo\n\
echo "--- SLURM job-script output ---"\n\
}}\n\
\n\
slurm_info_out\n\
\n\
slurm_startjob\n'

    fp.write(jobline)
    fp.close()


def scriptArrayjob(fname, gromppline, simpath, jobname, runtype, run, simtime=4, simcpu=1):
    fp = open(fname,'w')
    
    runargs = f'-ntmpi 1 -ntomp {simcpu} -pin on'

    jobline = f'#! /usr/bin/bash\n\
#SBATCH --job-name={jobname}\n\
#SBATCH --error={jobname}.e\n\
#SBATCH --array=1-80\n\
#SBATCH --partition=nes2.8\n\n\
#SBATCH --nodes=1\n\
#SBATCH --constraint=wes2.8\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --cpus-per-task={simcpu}\n\
#SBATCH --time={simtime}-00:00:00\n\
\n\
#Function to call to run the actual code\n\
slurm_startjob(){{\n\
#----------------- Actual calculation command goes here: ---------------------------\n\
\n\
module purge\n\
module load gnu/8.2.0  \n\
module load openmpi/3.1.2\n\
\n\
source /export/home/dfhahn/bin/gromacs_2019.4/bin/GMXRC\n\
PATH=/export/home/dfhahn/bin/gromacs_2019.4/bin/:$PATH\n\
\n\
WORKDIR=$(pwd)\n\
mkdir -p $SLURM_TMPDIR/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID\n\
cd $SLURM_TMPDIR/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID\n\
{gromppline}\n\n\
gmx mdrun {runargs} -s tpr.tpr \n\n\
rsync -avzIi dhdl.xvg {simpath}/dgdl$SLURM_ARRAY_TASK_ID.xvg\n\n\
rm -v $SLURM_TMPDIR/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID/*.*\n\n\
rmdir -v $SLURM_TMPDIR/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID\n\
}}\n\
\n\
slurm_info_out(){{\n\
\n\
echo "=================================== SLURM JOB ==================================="\n\
date\n\
echo\n\
echo "The job will be started on the following node(s):"\n\
echo $SLURM_JOB_NODELIST\n\
echo\n\
echo "Slurm User:         $SLURM_JOB_USER"\n\
echo "Run Directory:      $(pwd)"\n\
echo "Job ID:             ${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"\n\
echo "Job Name:           $SLURM_JOB_NAME"\n\
echo "Partition:          $SLURM_JOB_PARTITION"\n\
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"\n\
echo "Number of tasks:    $SLURM_NTASKS"\n\
echo "Submitted From:     $SLURM_SUBMIT_HOST:$SLURM_SUBMIT_DIR"\n\
echo "=================================== SLURM JOB ==================================="\n\
echo\n\
echo "--- SLURM job-script output ---"\n\
}}\n\
\n\
slurm_info_out\n\
\n\
slurm_startjob\n'

    fp.write(jobline)    
    fp.close()
