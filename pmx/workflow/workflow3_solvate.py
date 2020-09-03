#!/usr/bin/env python
# coding: utf-8

"""workflow3_solvate.py: Solvates simulation scripts and add ions."""

import os
import glob
import shutil
from pmx.workflow import pmxworkflow
import argparse
import subprocess

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

################################################
# prepare water boxes for simulations in water #
################################################
def solvateLigand(pwf):
    pmxworkflow.printInfo(runtype='solvateLigand',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        lig1 = item[0]
        lig2 = item[1]
        # make temporary directory tmp
        os.makedirs(f'{pwf.hybPath}/{edge}/water/tmp/', exist_ok=True)

        # input
        topology = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/topol.top'
        pmxworkflow.create_top(fname=topology, itp=['ffMOL.itp', 'merged.itp'])

        # create box
        # input
        inp = f'{pwf.hybPath}/{edge}/water/crd/mergedA.pdb'  # input ligand structure
        # output
        box = f'{pwf.hybPath}/{edge}/water/tmp/box.pdb'  # ligand placed in box; temporary file
        process = subprocess.Popen(['gmx', 'editconf',
                                    '-f', inp,
                                    '-o', box,
                                    '-bt', 'dodecahedron',
                                    '-d', str(1.5)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if pwf.verbose:
            out = process.communicate()
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))

        # solvate
        # output
        water = f'{pwf.hybPath}/{edge}/water/tmp/water.pdb'  # solvated ligand; temporary file
        process = subprocess.Popen(['gmx', 'solvate',
                                    '-cp', box,
                                    '-o', water,
                                    '-cs', 'spc216.gro',
                                    '-scale', str(1.0),
                                    '-p', topology],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if pwf.verbose:
            out = process.communicate()
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))

        # add ions, 2 sub-steps: a) grompp, b) add ions for 3 replicas independently
        # grompp
        mdpemA = f'{pwf.mdpPath}/em_stateA.mdp'
        mdout = f'{pwf.hybPath}/{edge}/water/tmp/mdout.mdp'  # temporary mdout file (not important)
        tprfile = f'{pwf.hybPath}/{edge}/water/tmp/tpr.tpr'  # temporary tpr file
        process = subprocess.Popen(['gmx', 'grompp',
                                    '-p', topology,
                                    '-c', water,
                                    '-o', tprfile,
                                    '-f', mdpemA,
                                    '-po', mdout,
                                    '-maxwarn', str(2)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if pwf.verbose:
            out = process.communicate()
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


################################################
# add ions for water boxes #
################################################
def addIonsLigand(pwf):
    pmxworkflow.printInfo(runtype='addIonsLigand',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge in pwf.edges.keys():
        topology = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/topol.top'
        tprfile = f'{pwf.hybPath}/{edge}/water/tmp/tpr.tpr'  # temporary tpr file

        # add ions
        for i in pwf.replicates:
            ions = f'{pwf.hybPath}/{edge}/water/crd/ions{i}.pdb'  # ion file
            topfile = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/topol{i}.top'  # top file
            shutil.copy(topology, topfile)
            process = subprocess.Popen(['gmx', 'genion',
                                        '-s', tprfile,
                                        '-p', topfile,
                                        '-conc', str(0.15),
                                        '-neutral',
                                        '-nname', 'ClJ',
                                        '-pname', 'NaJ',
                                        '-o', ions],
                                        stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            process.stdin.write('SOL'.encode())

            if pwf.verbose:
                out = process.communicate()
                print('STDERR{} '.format(out[1].decode("utf-8")))
                print('STDOUT{} '.format(out[0].decode("utf-8")))

def assembleComplex(file, molfiles=[]):
    fpout = open(file, 'w')
    for molfile in molfiles:
        fp = open(molfile, 'r')
        lines = fp.readlines()
        fp.close()

        for l in lines:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                fpout.write(l)
    fpout.close()


################################################
# assemble protein and ligand                  #
################################################
def combineProteinLigand(pwf):
    pmxworkflow.printInfo(runtype='Combine Protein and Ligand',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        # make directory for complex
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/crd/', exist_ok=True)
        # make temporary directory
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/tmp/', exist_ok=True)
        # make directory for complex topology
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}', exist_ok=True)

        # assemble coordinates
        # input
        proteinFile = f'{pwf.protPath}/crd/protein.pdb'  # input protein structure with ions, cofactors and crystal waters
        ligandFile = f'{pwf.hybPath}/{edge}/water/crd/mergedA.pdb'  # input ligand structure
        waterFile = f'{pwf.protPath}/crd/cofactors_crystalwater.pdb'  # input protein structure with ions, cofactors and crystal waters
        #output
        complexFile = f'{pwf.hybPath}/{edge}/complex/crd/complex.pdb' # complex of the former two structures

        if os.path.isfile(ligandFile):
            assembleComplex(complexFile, molfiles=[proteinFile, ligandFile, waterFile])
        else:
            print(f"No hybrid structure available. File {ligandFile} missing.")
            continue
 
        # assemble topologies
        proteinTop = f'{pwf.protPath}/top/amber99sb-star-ildn-mut.ff/protein.top' # protein topol        
        cofactors_crystalwaterTop = f'{pwf.protPath}/top/amber99sb-star-ildn-mut.ff/cofactors_crystalwater.top' # protein topology
        # output
        complexTop = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top'  # complex topology


        # included itp files
        # merged molecules itp
        includes = ['ffMOL.itp', 'merged.itp']
        # protein
        includes += list(filter(lambda x: x.endswith('.itp'), os.listdir(f'{pwf.protPath}/top/amber99sb-star-ildn-mut.ff/')))


        # molecules to be included
        molecules = []
        # protein system parts
        with open(proteinTop, 'r') as fp:
            lines = fp.readlines()
            for i, l in enumerate(lines):
                if l.startswith('[ molecules ]'):
                    break
            for l in lines[i+1:]:
                if l.startswith(';'):
                    continue
                molecules.append(l.split())
        # add ligand molecule
        molecules.append(['MOL', 1])
        # add cofactor and water system parts
        with open(cofactors_crystalwaterTop, 'r') as fp:
            lines = fp.readlines()
            for i, l in enumerate(lines):
                if l.startswith('[ molecules ]'):
                    break
            for l in lines[i+1:]:
                if l.startswith(';'):
                    continue
                molecules.append(l.split())

        pmxworkflow.create_top(fname=f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top', ff='amber99sb-star-ildn-mut.ff', water='tip3p',
                   itp=includes, mols=molecules,
                   destination=f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/', toppaths=[f'{pwf.protPath}/top/amber99sb-star-ildn-mut.ff/', f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/'])

################################################
# prepare water boxes for complex simulations  #
################################################
def solvateComplex(pwf):
    pmxworkflow.printInfo(runtype='solvateComplex',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        # make directory for complex
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/crd/', exist_ok=True)
        # make temporary directory
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/tmp/', exist_ok=True)

        # create box
        # input
        inp = f'{pwf.hybPath}/{edge}/complex/crd/complex.pdb'  # input complex structure
        topology = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top'
        # output
        box = f'{pwf.hybPath}/{edge}/complex/tmp/box.pdb'  # complex placed in box, temporary file
        process = subprocess.Popen(['gmx', 'editconf',
                                    '-f', inp,
                                    '-o', box,
                                    '-bt', 'dodecahedron',
                                    '-d', str(1.5)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if pwf.verbose:
            out = process.communicate()
            print('STDOUT{} '.format(out[0].decode("utf-8")))
            print('STDERR{} '.format(out[1].decode("utf-8")))

        # solvate
        water = f'{pwf.hybPath}/{edge}/complex/tmp/water.pdb'  # solvated complex, temporary file
        process = subprocess.Popen(['gmx', 'solvate',
                                    '-cp', box,
                                    '-o', water,
                                    '-cs', 'spc216.gro',
                                    '-scale', str(1.0),
                                    '-p', topology],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if pwf.verbose:
            out = process.communicate()
            print('STDOUT{} '.format(out[0].decode("utf-8")))
            print('STDERR{} '.format(out[1].decode("utf-8")))

        # grompp
        mdpemA = f'{pwf.mdpPath}/em_stateA.mdp'
        mdout = f'{pwf.hybPath}/{edge}/complex/tmp/mdout.mdp'  # temporary mdout file (not important)
        tprfile = f'{pwf.hybPath}/{edge}/complex/tmp/tpr.tpr'  # temporary tpr file
        process = subprocess.Popen(['gmx', 'grompp',
                                    '-p', topology,
                                    '-c', water,
                                    '-o', tprfile,
                                    '-f', mdpemA,
                                    '-po', mdout,
                                    '-maxwarn', str(3)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if pwf.verbose:
            out = process.communicate()
            print('STDOUT{} '.format(out[0].decode("utf-8")))
            print('STDERR{} '.format(out[1].decode("utf-8")))


################################################
# add ions for protein boxes #
################################################
def addIonsComplex(pwf):
    pmxworkflow.printInfo(runtype='addIonsComplex',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge in pwf.edges.keys():
        # make directory for complex
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/', exist_ok=True)
        # make temporary directory
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/tmp/', exist_ok=True)

        topology = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top'
        tprfile = f'{pwf.hybPath}/{edge}/complex/tmp/tpr.tpr'  # temporary tpr file
        if not os.path.isfile(topology):
            print(f"Topology file {topfile} does not exist.")
            continue

        # add ions
        for i in pwf.replicates:
            ions = f'{pwf.hybPath}/{edge}/complex/crd/ions{i}.pdb'  # ion file
            topfile = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol{i}.top'  # top file
            shutil.copy(topology, topfile)

            process = subprocess.Popen(['gmx', 'genion',
                                        '-s', tprfile,
                                        '-p', topfile,
                                        '-conc', str(0.15),
                                        '-neutral',
                                        '-nname', 'ClJ',
                                        '-pname', 'NaJ',
                                        '-o', ions],
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.stdin.write('SOL'.encode())

            if pwf.verbose:
                out = process.communicate()
                print('STDOUT{} '.format(out[0].decode("utf-8")))
                print('STDERR{} '.format(out[1].decode("utf-8")))

# clean unnecessary files
def cleanUp(pwf):
    pmxworkflow.printInfo(runtype='Cleaning up',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')

        # clean up
        toclean = glob.glob(f'{pwf.hybPath}/{edge}/water/tmp/*')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/complex/tmp/*')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/water/crd/*#')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/*#')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/complex/crd/*#')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/*#')
        for clean in toclean:
            os.remove(clean)
        os.rmdir(f'{pwf.hybPath}/{edge}/water/tmp/')
        os.rmdir(f'{pwf.hybPath}/{edge}/complex/tmp/')



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
                        choices=['smirnoff99Frosst-1.1.0.offxml', 'openff-1.0.0.offxml', 'gaff2'],
                        help='The force field used.')
    parser.add_argument('-p',
                        '--path',
                        metavar='PATH',
                        type=str,
                        default='./',
                        help='The path to the data.')
    parser.add_argument('-r',
                        '--replicates',
                        metavar='REPLICATES',
                        nargs='+',
                        type=str,
                        default='1',
                        help='Comma or space separated list of integers or range of integers denoted by hyphen')
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

    replicates = args.replicates
    replicates = [n for r in replicates for n in r.split(',') if n is not '']
    replicates = [n for r in replicates for n in ' - '.join(r.split('-')).split()]
    replicatesToUse = []
    for i, r in enumerate(replicates):
        if r == '-':
            replicatesToUse += range(int(replicates[i - 1]) + 1, int(replicates[i + 1]), 1)
        else:
            replicatesToUse.append(int(r))

    pwf = pmxworkflow.pmxvariables(target=args.target,
                                   forcefield=args.forcefield,
                                   path=args.path,
                                   replicates=replicates,
                                   verbose=args.verbose)

    if args.edges[0] != 'all':
        edges = dict()
        for edge in args.edges:
            edges[edge] = pwf.edges[edge]
        pwf.edges = edges

    solvateLigand(pwf)
    addIonsLigand(pwf)

    combineProteinLigand(pwf)

    solvateComplex(pwf)
    addIonsComplex(pwf)

    cleanUp(pwf)


