# from schrodinger.application.desmond.packages import analysis, traj, topo, traj_util, msys
from statistics import mean, stdev
import matplotlib.pyplot as plt
from os.path import exists
import pickle
import sys

#-----------------------------------------------------------------------
## Script to generate SASA replica std-dev rabge plots for DNA amph chem paper ##
## user must suply job filenames and paths ##
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# function to generate moving average of list
def moving_avg(list):
    N = 3 # moving average of N
    cumsum, moving_aves = [0], []
    for val in list[:6]:
        moving_aves.append(val)

    for i, x in enumerate(list[N:], 1):
        cumsum.append(cumsum[i-1] + x)
        if i>=N:
            moving_ave = (cumsum[i] - cumsum[i-N])/N
            moving_aves.append(moving_ave)
    return moving_aves

#-----------------------------------------------------------------------
# function to get relative results
def relative_sasa(list):
    to_return = []
    for item in list:
        to_return.append(item/list[0])
    return(to_return)

#-----------------------------------------------------------------------
# define jobnames and path
files = ["R", "RS"]
reps = ["", "_rep", "_rep_2"]
path = "/home/marlo/HDD/desmond_jobs/jobs/hanadi/C8/"
fig_name = "test"
job_len = 100 #ns

# initalize lists
avgd_results = []
std_devs = []

#-----------------------------------------------------------------------
# iterate thru stereos
for i, file in enumerate(files):
        replica_results = []
        tmp = []
        dev_tmp = []

        #-----------------------------------------------------------------------
        # iterate over replicas
        for rep in reps:
            filename = file + rep
            # check for pickle file
            #-----------------------------------------------------------------------
            if exists("data/" + filename + "_SASA.p"):
                results = pickle.load(open("data/" + filename + "_SASA.p", "rb"))
                replica_results.append([i/results[0] for i in results])
            else:
                location_in = path + filename + "/" + filename + "-in.cms"
                location_traj = path + filename + "/" + filename + "_trj/"

                msys_model, cms_model = topo.read_cms(location_in)
                trj = traj.read_traj(location_traj)[::10] # full trajectories take a long time, for testing try skipping frames
                # trj = trj[::int(len(trj)/300)]

                sasa = analysis.SolventAccessibleSurfaceArea(msys_model, cms_model, asl="chain A", resolution=0.2)
                trj_results = analysis.analyze(trj, sasa)

                #-----------------------------------------------------------------------
                # generate and add init sasa
                # we use the unrelaxes SASA here because it is a better representation of the max sasa
                location = path + filename + "/" + filename + ".cms"
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
                pickle.dump(results,  open("data/" + filename + "_SASA.p", "wb"))
                replica_results.append([i/results[0] for i in results])

        # # iterate thru replicates
        # for rep in reps:
        #     # get data
        #     label = file + rep + "_SASA"
        #     results = pickle.load(open("data/" + label + ".p", "rb"))
        #     replica_results.append([i/results[0] for i in results])

        # get standard deviations and averages
        for i, val in enumerate(replica_results[0]):
            tmp.append((val + replica_results[1][i])/2)
            dev_tmp.append(stdev([val, replica_results[1][i]]))
        avgd_results.append(tmp)
        std_devs.append(dev_tmp)

#-----------------------------------------------------------------------
# get high and lows from std devs
high_results = []
low_results = []
for i, set in enumerate(avgd_results):
    tmp_high = []
    tmp_low = []
    for j, val in enumerate(set):
        tmp_high.append(val + std_devs[i][j]/2)
        tmp_low.append(val - std_devs[i][j]/2)
    high_results.append(tmp_high)
    low_results.append(tmp_low)

#-----------------------------------------------------------------------
# smooth data with a moving avg
RS_high = moving_avg(high_results[0])
RS_low = moving_avg(low_results[0])
R_high = moving_avg(high_results[1])
R_low = moving_avg(low_results[1])

# construct x axis
time_per_frame = job_len/len(RS_high)
time = []
for i in range(len(RS_high)):
  time.append(i*time_per_frame)

#-----------------------------------------------------------------------
# plot settings
plt.rc('xtick', labelsize=30)
plt.rc('ytick', labelsize=30)
fig = plt.figure(dpi=300, figsize=(15,11))
ax1 = fig.add_subplot(111)

# plot lines
ax1.plot(time, RS_high,
        time, RS_low, color = "dimgrey")

ax1.plot(time, R_high,
         time, R_low, color = "coral")

# fill between lines
RS_label = plt.fill_between(time, RS_high, RS_low, color="dimgrey", label="RS$_{6}$ SASA")
R_label = plt.fill_between(time, R_high, R_low, color="coral", label="R$_{12}$ SASA")

#-----------------------------------------------------------------------
# finalize graphs
# plt.xlabel("Time (ns)", weight='bold', size="xx-large")
# plt.ylabel("Surface Area (A$^2$)", WEIGHT='bold', size='xx-large')
lgnd = plt.legend(handles=[RS_label, R_label], labels=["RS$_{6}$ SASA", "R$_{12}$ SASA"], loc='upper right', prop={'size': 35})
lgnd.get_frame().set_edgecolor('black')
plt.savefig("../figs/"+fig_name+".png", transparent=True, bbox_inches='tight')
# plt.show()
