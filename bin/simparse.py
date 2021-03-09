# -*- coding: utf-8 -*-
"""
#parsing treesequence results from SliM spatial simulation
#Tristan Dennis 05.02.21
#tristanpwdennis@gmail.com
Liberally adapted from the tskit spatial vignette at https://pyslim.readthedocs.io/en/stable/vignette_space.html#
"""
import pyslim 
import numpy as np
import msprime
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import plotnine as p9
#load spatial trees
slim_ts = pyslim.load("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/sims/spatial_sim.trees")
print(f"The tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
      f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
      f" and {slim_ts.num_mutations} mutations.")

#now we use msprime to recapitate, mutate, and save updated treesequence that contains mutation data
recap_ts = slim_ts.recapitate(recombination_rate=1e-8, Ne=1000)
ts = pyslim.SlimTreeSequence(
      msprime.mutate(recap_ts, rate=1e-8, keep=True))
ts.dump("spatial_sim.recap.trees")
print(f"The tree sequence now has {ts.num_trees} trees,"
      f" and {ts.num_mutations} mutations.")

#now let's plot our individuals over their region
#set random seed
np.random.seed(23)

#get individuals alive 30 gens ago (before the disaster)
old_ones = ts.individuals_alive_at(30)
#get individuals alive today
new_ones = ts.individuals_alive_at(0)

#make a dictionary of the groups we are sampling
#can extend to multiple groups of any kind 
#here, random sample across all our individuals

groups = {
    'ancient' : np.random.choice(old_ones, size=70),
    'today' : np.random.choice(new_ones, size=70)
    }

#print some group info
for g in groups:
   print(f"We have {len(groups[g])} individuals in the {g} group.")

#order groups   
group_order = ['ancient', 'today']

#assign colours to groups
ind_colors = np.repeat(0, ts.num_individuals)
for j, k in enumerate(group_order):
   ind_colors[groups[k]] = 1 + j
   
   
   
   
   

#get locations of old and new ones
old_locs = ts.individual_locations[old_ones, :]
today_locs = ts.individual_locations[new_ones, :]



fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("today")
ax.scatter(today_locs[:,0], today_locs[:,1], s=20, c=ind_colors[new_ones])
ax = fig.add_subplot(122)
ax.set_title("long ago")
ax.scatter(old_locs[:, 0], old_locs[:, 1], s=20, c=ind_colors[old_ones])
fig.savefig("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/plots/spatial_sim_locations.png")


#now to calculate relatedness and geographic distance between each individual
#create list of nodes where each entry is the nodes belonging to that individual, being mindful of which group
#belongs to which
ind_nodes = []
ind_group = []
ind_ids = []
for j, group in enumerate(group_order):
   for ind in groups[group]:
      ind_ids.append(ind)
      ind_nodes.append(ts.individual(ind).nodes)
      ind_group.append(group_order[j])

nind = len(ind_ids)
pairs = [(i, j) for i in range(nind) for j in range(i, nind)]
ind_div = ts.divergence(ind_nodes, indexes=pairs)
#then create a list of paired individuals
#and use this as input to query the nodes with ts.divergence
nind = len(ind_ids)
pairs = [(i, j) for i in range(nind) for j in range(i, nind)]
ind_div = ts.divergence(ind_nodes, indexes=pairs)

#likewise with geographic distances
geog_dist = np.repeat(0.0, len(pairs))
locs = ts.individual_locations
for k, (i, j) in enumerate(pairs):
   geog_dist[k] = np.sqrt(np.sum((locs[ind_ids[i], :] - locs[ind_ids[j], :])**2))
   
#check to make sure the distance of individuals from themselves is zero
for (i, j), x in zip(pairs, geog_dist):
  if i == j:
    assert(x == 0)
    
#create big df of all our data
geodf = []
geodf = pd.DataFrame(data=[pairs, ind_div, geog_dist]).T
geodf.columns = ['pair_id', 'ind_div', 'geog_dist']
geodf[['ind1', 'ind2']] = pd.DataFrame(geodf['pair_id'].tolist(), index=geodf.index)
geodf['ind1'] = 'tsk_' + geodf['ind1'].astype(str)
geodf['ind2'] = 'tsk_' + geodf['ind2'].astype(str)
geodf.to_csv("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/metadata/isolation_by_distance.txt")  

#export ind metadata to txt file
indivlist = []
indivnames = []
with open("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/metadata/spatial_sim_individuals.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label", "tskit_id", "slim_id"]
                               + ["birth_time_ago", "age", "x", "y"]) + "\n")
  for group in group_order:
     for i in groups[group]:
        indivlist.append(i)
        ind = ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(ind.id), str(ind.metadata["pedigree_id"]), str(ind.time),
                str(ind.metadata["age"]), str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")


pair_colors = np.repeat(0, len(pairs))
for k, (i, j) in enumerate(pairs):
   if ind_group[i] == "ancient" or ind_group[j] == "ancient":
      pair_colors[k] = 1

fig = plt.figure(figsize=(6, 6), dpi=300)
ax = fig.add_subplot(111)
ax.scatter(geog_dist, 1e3 * ind_div, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("genetic distance (diffs/Kb)")
fig.savefig("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/plots/spatial_sim_ibd.png")


