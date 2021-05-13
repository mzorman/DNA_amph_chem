from schrodinger.application.desmond.packages import analysis, traj, topo, traj_util, msys
from schrodinger.structutils.analyze import calculate_sasa
from schrodinger import structure
import matplotlib.pyplot as plt
from statistics import mean
from os.path import exists
import pickle
import sys

#-------------------------------------------------------------------------------
## Script to generate SASA plots for DNA amph chem paper ##
## user must suply job filenames and paths ##
#-------------------------------------------------------------------------------

# define jobnames and path
files = ["RS_jobname", "R_jobname"]
path = "/home/you/wherever/you/run/jobs/"
fig_name = "test"
job_len = 100 #ns

#-------------------------------------------------------------------------------
# iterate through files and get SASA
all_results = []
for file in files:
    print(file)
    # check for pickle file
    if exists("data/" + file + "_SASA.p"):
        results = pickle.load(open("data/" + file + "_SASA.p", "rb"))
    else:
        location_in = path + file + "/" + file + "-in.cms"
        location_traj = path + file + "/" + file + "_trj/"

        msys_model, cms_model = topo.read_cms(location_in)
        trj = traj.read_traj(location_traj)[::10] # full trajectories take a long time, for testing try skipping frames
        # trj = trj[::int(len(trj)/300)]

        sasa = analysis.SolventAccessibleSurfaceArea(msys_model, cms_model, asl="chain A", resolution=0.2)
        trj_results = analysis.analyze(trj, sasa)

        #-----------------------------------------------------------------------
        # generate and add init sasa
        # we use the unrelaxes SASA here because it is a better representation of the max sasa
        location = path + file + "/" + file + ".cms"
        st = next(structure.StructureReader(location))
        indeces = []
        for chain in st.chain:
            if chain.name == "A":
                for atom in chain.atom:
                    indeces.append(atom.index)
        init_sasa = calculate_sasa(st, atoms=indeces, resolution=0.2, exclude_water=True)
        results = [init_sasa] + trj_results

        #-----------------------------------------------------------------------
        # pickle results
        pickle.dump(results,  open("data/" + file + "_SASA.p", "wb"))
    all_results.append(results)

#-------------------------------------------------------------------------------
# plot relative results
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
fig = plt.figure(dpi=300, figsize=(15,11))
ax1 = fig.add_subplot(111)
time = [i*(job_len/len(all_results[0])) for i in range(len(all_results[0]))]

for results in all_results:
    results = [i/results[0] for i in results]
    plt.plot(time, results, linewidth=3.0)

plt.xlabel("Time (ns)", weight='bold', fontsize=20)
plt.ylabel("Relative SASA", WEIGHT='bold', fontsize=20)
plt.legend(files, loc='best', prop={'size': 10})
plt.title(fig_name)
plt.savefig("../figs/"+fig_name+".png", bbox_inches='tight')
