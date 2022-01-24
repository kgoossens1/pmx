#!/usr/bin/env python
# coding: utf-8

"""workflow2_make_perturbations.py: Creates hybrid structures and topologies."""

import os
from pmx.workflow import pmxworkflow
import argparse
import subprocess

from openforcefield.utils import toolkits

### OpenEye version: uncomment the following if you have and if you want to use the OpenEye toolkit, then RDKit and Ambertools toolkits
### Attention: there could be problems in the workflow due to atom order changes, this is related to this issue: https://github.com/openforcefield/openforcefield/issues/475
#toolkit_precedence = [toolkits.OpenEyeToolkitWrapper, toolkits.RDKitToolkitWrapper, toolkits.AmberToolsToolkitWrapper]

### Non-OpenEye version: uncomment the following if you want to use the rdkit and ambertools
toolkit_precedence = [toolkits.RDKitToolkitWrapper, toolkits.AmberToolsToolkitWrapper]

toolkits.GLOBAL_TOOLKIT_REGISTRY = toolkits.ToolkitRegistry(toolkit_precedence=toolkit_precedence)

from openforcefield.topology import Molecule

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def atomsToMorph(pwf):
    pmxworkflow.printInfo(runtype='atomsToMorph',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        lig1 = item[0]
        lig2 = item[1]

        # make temporary directory tmp
        os.makedirs(f'{pwf.hybPath}/{edge}/water/tmp/', exist_ok=True)

        # input
        pdb1 = f'{pwf.hybPath}/{edge}/water/tmp/{lig1}.pdb'  # pdb structure lig1
        pdb2 = f'{pwf.hybPath}/{edge}/water/tmp/{lig2}.pdb'  # pdb structure lig2

        # conversion from sdf to pdb
        ligPath1 = f'{pwf.ligPath}/{lig1}/crd/{lig1}.sdf'
        ligand = Molecule.from_file(f'{ligPath1}', allow_undefined_stereo=True)
        ligand.to_file(pdb1, 'pdb')

        ligPath2 = f'{pwf.ligPath}/{lig2}/crd/{lig2}.sdf'
        ligand = Molecule.from_file(f'{ligPath2}', allow_undefined_stereo=True)
        ligand.to_file(pdb2, 'pdb')

        # output
        o1 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb1.pdb'  # processed structure lig1; temporary file
        o2 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb2.pdb'  # processed structure lig2; temporary file
        pairs = f'{pwf.hybPath}/{edge}/water/crd/pairs.dat'  # atom mapping
        score = f'{pwf.hybPath}/{edge}/water/crd/score.dat'  # score of the mapping

        # make coordinate directory crd
        os.makedirs(f'{pwf.hybPath}/{edge}/water/crd/', exist_ok=True)

        command = ['python3', pwf.scriptpath + '/atoms_to_morph.py',
                                    '-i1', pdb1,
                                    '-i2', pdb2,
                                    '-opdb1', o1,
                                    '-opdb2', o2,
                                    '-score', score,
                                    '-o', pairs,
                                    '-H2H',
                                    '-timeout', str(30)]
        print(f'COMMAND: {" ".join(command)}')
        # run atoms_to_morph.py script
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


# create hybrid structures/topologies
def makeHybrid(pwf):
    pmxworkflow.printInfo(runtype='makeHybrid',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        lig1 = item[0]
        lig2 = item[1]
        # make temporary directory tmp
        os.makedirs(f'{pwf.hybPath}/{edge}/water/tmp/', exist_ok=True)
        # input
        l1 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb1.pdb'  # structure lig1; temporary file
        l2 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb2.pdb'  # structure lig2; temporary file
        pairs = f'{pwf.hybPath}/{edge}/water/crd/pairs.dat'  # atom mapping
        itp1 = f'{pwf.ligPath}/{lig1}/top/{pwf.forcefield}/{lig1}.itp'  # topology lig1
        itp2 = f'{pwf.ligPath}/{lig2}/top/{pwf.forcefield}/{lig2}.itp'  # topology lig2
        # output
        opdbA = f'{pwf.hybPath}/{edge}/water/crd/mergedA.pdb'  # hybrid structure based on lig1
        opdbB = f'{pwf.hybPath}/{edge}/water/crd/mergedB.pdb'  # hybrid structure based on lig2
        oitp = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/merged.itp' # hybrid topology
        offitp = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffmerged.itp'  # force field parameters for dummies
        olog = f'{pwf.hybPath}/{edge}/water/tmp/hybrid.log'  # output log


        if not os.path.exists(l1):
            ligPath1 = f'{pwf.ligPath}/{lig1}/crd/{lig1}.sdf'
            ligand = Molecule.from_file(f'{ligPath1}', allow_undefined_stereo=True)
            ligand.to_file(l1, 'pdb')
        if not os.path.exists(l2):
            ligPath2 = f'{pwf.ligPath}/{lig2}/crd/{lig2}.sdf'
            ligand = Molecule.from_file(f'{ligPath2}', allow_undefined_stereo=True)
            ligand.to_file(l2, 'pdb')

        # make top directory
        os.makedirs(f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/', exist_ok=True)

        command = ['python', pwf.scriptpath + '/make_hybrid.py',
                                    '-pairs', pairs,
                                    '-l1', l1,
                                    '-l2', l2,
                                    '-itp1', itp1,
                                    '-itp2', itp2,
                                    '-oitp', oitp,
                                    '-oa', opdbA,
                                    '-ob', opdbB,
                                    '-ffitp', offitp,
                                    '-log', olog,
                                    '-scDUMm', str(0.001)]

        print(f'COMMAND: {" ".join(command)}')
        # creates hybrid structure/topology
        process = subprocess.Popen(command, 
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


# deal with the atomtype topology files
def oneffFile(pwf):
    pmxworkflow.printInfo(runtype='oneffFile',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        lig1 = item[0]
        lig2 = item[1]

        # input
        ff1 = f'{pwf.ligPath}/{lig1}/top/{pwf.forcefield}/ff{lig1}.itp'  # atomtypes lig1
        ff2 = f'{pwf.ligPath}/{lig2}/top/{pwf.forcefield}/ff{lig2}.itp'  # atomtypes lig2
        ffmerged = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffmerged.itp'  # force field parameters for dummies
        # output
        ffout = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffMOL.itp'  # atomtypes out

        command = ['python', pwf.scriptpath + '/one_ff_file.py',
                                    '-ffitp', ff1, ff2, ffmerged,
                                    '-ffitp_out', ffout]

        print(f'COMMAND: {" ".join(command)}')
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


# clean unnecessary files
def cleanUp(pwf):
    pmxworkflow.printInfo(runtype='Cleaning up',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        # clean up
        toclean = os.listdir(f'{pwf.hybPath}/{edge}/water/tmp/')
        for file in toclean:
            os.remove(f'{pwf.hybPath}/{edge}/water/tmp/{file}')
        os.rmdir(f'{pwf.hybPath}/{edge}/water/tmp/')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',
                        '--target',
                        metavar='TARGET',
                        type=str,
                        default='01_jnk1',
                        help='The target protein.')
    parser.add_argument('-f',
                        '--forcefield',
                        metavar='FORCEFIELD',
                        type=str,
                        default='smirnoff99Frosst-1.1.0.offxml',
#                        choices=['smirnoff99Frosst-1.1.0.offxml', 'openff-1.0.0.offxml', 'openff-1.2.0.offxml', 'gaff2'],
                        help='The force field used.')
    parser.add_argument('-p',
                        '--path',
                        metavar='PATH',
                        type=str,
                        default='./',
                        help='The path to the data.')
    parser.add_argument('-e',
                        '--edges',
                        metavar='EDGES',
                        nargs='+',
                        type=str,
                        default=['all'],
                        help='Either "all" or a list of edges to be calculated.')
    parser.add_argument('-v',
                        '--verbose',
                        type=bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Turn on verbose output.')
    args = parser.parse_args()

    pwf = pmxworkflow.pmxvariables(target=args.target,
                                   forcefield=args.forcefield,
                                   path=args.path,
                                   replicates=[1],
                                   verbose=args.verbose)

    if args.edges[0] != 'all':
        edges = dict()
        for edge in args.edges:
            edges[edge] = pwf.edges[edge]
        pwf.edges = edges

    atomsToMorph(pwf)
    makeHybrid(pwf)
    oneffFile(pwf)
    cleanUp(pwf)



