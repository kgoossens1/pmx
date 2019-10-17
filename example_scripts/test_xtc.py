from __future__ import print_function
from pmx.parser import *
from pylab import *
from pmx import *
from pmx.xtc import *
import sys
import os

m = Model("out.pdb")


atom1 = m.residues[0]["CA"]  # CA-atom of first residues
atom2 = m.residues[-1]["CA"]  # CA-atom of last residue
# function that returns the distance between atom1 and atom2


def distance_12(a, b):
    return a - b


fp = open("analysis.dat", "w")  # open output file


trj = Trajectory("concoord.xtc")  # open xtc file and read first frame
for frame in trj:  # go over each frame
    print(frame)  # print some info
    trj.update(m)  # update coords in model
    d = distance_12(atom1, atom2)  # calculate observable
    print("%8.3f %8.3f" % (frame.time.value, d), file=fp)  # store data
fp.close()  # close output file
trj.close_xtc()  # close xtc
data = read_and_format("analysis.dat", "ff")
time = [a[0] for a in data]
dist = [a[1] for a in data]
plot(time, dist, "r-", lw=2)
xlabel("time [ps]")
ylabel(r"end to end distance [$\AA$]")
savefig("plot.png")
