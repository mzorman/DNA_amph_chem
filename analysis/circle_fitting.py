from schrodinger.application.desmond.packages import traj, topo
from schrodinger.application.desmond.packages import msys
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from analysis_utils import *
from os.path import exists
from numpy import *
import pickle
import sys

#-------------------------------------------------------------------------------
## Script used to generate monomeric units from DNA amphiphile conformations ##
## User must define filenames/paths, script will output time-plot of scores ##
## and csv file of top five scores and respective time (in ps) ##
#-------------------------------------------------------------------------------

# define job names and csv file
to_save = open("data/monomeric_unit_scores.csv", "w")
files = ["yourjobname"]
path = "/home/you/wherever/you/run/jobs/"
job_len = 100


#-------------------------------------------------------------------------------
# iterate through jobs and get top candidates
for file in files:
    # check for pickle file
    if exists("data/" + file + "_scores.p"):
        rmsds, times = pickle.load(open("data/" + file + "_scores.p", "rb"))
    # else generate data
    else:
        # location_ins
        location_in = path + file + "/" + file + "-in.cms"
        location_traj = path + file + "/" + file + "_trj/"

        # load cms and traj
        _, cms_model = topo.read_cms(location_in)
        tr = traj.read_traj(location_traj)

        # extract alkyl-portion (chain A) phosphorous backbone atoms
        aids = cms_model.select_atom('atom.ele P and chain.name A')
        st = cms_model.extract(aids)
        gids = topo.aids2gids(cms_model, aids, include_pseudoatoms=False)
        
        # iterate through trajectory
        rmsds = []
        times = []
        for fr in tr:
            # update and get positions
            st.setXYZ(fr.pos(gids))
            posns = []
            for atom in st.atom:
                posns.append(atom.xyz)

            # append score to list
            rmsds.append(get_rmsd(posns))
            times.append(fr.time)

        #---------------------------------------------------------------------------
        # calculate minimums
        rmsds_copy = [i for i in rmsds]
        to_save.write(file+"\n")
        n = 5

        aids = cms_model.select_atom('chain.name A or chain.name B')
        st = cms_model.extract(aids)
        gids = topo.aids2gids(cms_model, aids, include_pseudoatoms=False)

        # top n minimums
        for i in range(n):
            # get min and index
            rmsd_min = min(rmsds_copy)[0]
            rmsds_copy.pop(rmsds_copy.index(rmsd_min))
            index = rmsds.index(rmsd_min)

            # save structure and write time:score to csv
            to_save.write(str(index) + ": " + str(round(rmsd_min, 2)) + "\n")
            st.setXYZ(tr[index].pos(gids))
            st.write("../structs/monomeric_units/"+file+"_frame_"+str(index)+".pdb")
        to_save.write("\n")

        #---------------------------------------------------------------------------
        # pickle data
        pickle.dump([rmsds, times], open("data/" + file + "_scores.p", "wb"))

    #---------------------------------------------------------------------------
    # plot relative scores
    rel_rmsds = [i/rmsds[0] for i in rmsds]
    plt.plot([i/1000 for i in times], rel_rmsds)
    plt.xlabel("Time (ns)")
    plt.ylabel("Relative Score")
    plt.title(file + " circle scores")
    plt.savefig("../figs/" + file + "_circle_score.png")
    plt.close()

to_save.close()
