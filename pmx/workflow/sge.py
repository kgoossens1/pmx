#!/usr/bin/env python
# coding: utf-8

"""sge.py: Handles SGE queue system commands."""


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
#$ -N {jobname}\n\
#$ -cwd\n\
#$ -q all.q\n\
#$ -j yes\n\
#$ -l slot_type=gromacsp3,affinity_group=default\n\n\
source /home/dhahn3/bin/gromacs-patched/bin/GMXRC\n\n\
mkdir -p $TMPDIR/$JOB_ID/\n\
cd $TMPDIR/$JOB_ID/\n\n'
    fp.write(jobline)
    fp.close()

def scriptMain( fname, gromppline, simpath, runtype, run, simtime=4, simcpu=1 ):
    fp = open(fname, 'a')

    jobline = f'{gromppline}\n\
gmx mdrun -ntmpi 1 -ntomp {simcpu} -pme gpu -pin on -s {runtype}{run}.tpr -deffnm {runtype}{run}\n\
rsync -avzIi ./ {simpath}/\n\
rm -v $TMPDIR/$JOB_ID/*.*\n\n'
    fp.write(jobline)
    fp.close()


def scriptFooter( fname ):
    fp = open(fname, 'a')

    jobline = f'\n\
rmdir -v $TMPDIR/$JOB_ID\n\
            '
    fp.write(jobline)
    fp.close()
