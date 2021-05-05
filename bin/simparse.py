# -*- coding: utf-8 -*-
"""
#parsing treesequence results from SliM spatial simulation
#Tristan Dennis 05.02.21
#tristanpwdennis@gmail.com
"""
import pyslim 
import pandas as pd
import tskit 
import numpy as np
import msprime
import os
import plotnine as p9
import itertools as it
import seaborn as sns

os.chdir("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/")

#load spatial trees
slim_ts = pyslim.load("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/everythingisterribleafter30.trees")
print(f"The tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
      f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
      f" and {slim_ts.num_mutations} mutations.")

for t in np.unique(slim_ts.individual_times):
  print(f"There are {np.sum(slim_ts.individual_times == t)} individuals from time {t}.")


for t in [0, 1, 2, 3, 30]:
   alive = slim_ts.individuals_alive_at(t)
   print(f"There were {len(alive)} individuals alive {t} time steps in the past.")
  
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
#get all living inds
alive = ts.individuals_alive_at(0)
old = ts.individuals_alive_at(30)

locs = ts.individual_locations[alive, :]
old_locs = ts.individual_locations[old, :]

W = 35
w = 15
groups = {
   'topleft' : alive[np.logical_and(locs[:, 0] < w, locs[:, 1] < w)],
   'topright' : alive[np.logical_and(locs[:, 0] < w, locs[:, 1] > W - w)],
   'bottomleft' : alive[np.logical_and(locs[:, 0] > W - w, locs[:, 1] < w)],
   'bottomright' : alive[np.logical_and(locs[:, 0] > W - w, locs[:, 1] > W - w)],
   'center' : alive[np.logical_and(np.abs(locs[:, 0] - W/2) < w/2, np.abs(locs[:, 1] - W/2) < w/2)],
   'oldtopleft' : old[np.logical_and(old_locs[:, 0] < w, old_locs[:, 1] < w)],
   'oldtopright' : old[np.logical_and(old_locs[:, 0] < w, old_locs[:, 1] > W - w)],
   'oldbottomleft' : old[np.logical_and(old_locs[:, 0] > W - w, old_locs[:, 1] < w)],
   'oldbottomright' : old[np.logical_and(old_locs[:, 0] > W - w, old_locs[:, 1] > W - w)],
   'oldcenter' : old[np.logical_and(np.abs(old_locs[:, 0] - W/2) < w/2, np.abs(old_locs[:, 1] - W/2) < w/2)]
   }


groups = {
    'pretx' : np.random.choice(alive, size=100),
    'posttx' : np.random.choice(old, size=100)
    }

for k in groups:
   print(f"We have {len(groups[k])} individuals in the {k} group.")
   
#get all our sampled inds and their group names
group_memb = []
for a, group in enumerate(groups):
    for ind in groups[group]:
        row = [ind, group]
        group_memb.append(row)
        
group_memb = pd.DataFrame(group_memb, columns = ['ind', 'group'])

#ok now let's concatenate the np arrays and coerce to dataframe for old and new
now = pd.DataFrame(np.column_stack((locs, alive)))
then = pd.DataFrame(np.column_stack((old_locs, old)))

#add the times in a col
now['time'] = "now"
then['time'] = 'then'

metadata = now.append(then) #merge dfs
metadata.columns = ['x', 'y', 'z', 'ind', 'time'] #add colnames
metadata['ind'] = pd.to_numeric(metadata['ind'], downcast='integer') #coerce ind id to integer from float
metadata = metadata.merge(group_memb, how = 'left') #jopin with group data
metadata = metadata.fillna('unsampled') #turn nans into 'unsampled'

plot = (p9.ggplot(metadata, p9.aes('x', 'y', color='group'))
 + p9.geom_point()
  + p9.theme_bw()
  + p9.stat_density_2d()
 + p9.facet_wrap('~time'))
plot


metadata.to_csv("everythingisterribleafter30.csv")

group_order = ['pretx', 'posttx']
sampled_nodes = [[] for _ in groups]
for j, k in enumerate(group_order):
   for ind in groups[k]:
      sampled_nodes[j].extend(ts.individual(ind).nodes)
      
#create lists of ids, nodes and groups for our sampled inds
ind_nodes = []
ind_group = []
ind_ids = []
for j, group in enumerate(group_order):
   for ind in groups[group]:
      ind_ids.append(ind)
      ind_nodes.append(ts.individual(ind).nodes)
      ind_group.append(group_order[j])

#for each pair calculate genetic distance
nind = len(ind_ids)
pairs = [(i, j) for i in range(nind) for j in range(i, nind)] #create a list of indexes to query the list of nodes with
ind_div = ts.divergence(ind_nodes, indexes=pairs)

#for each pair calculate the geographic distance
geog_dist = np.repeat(0.0, len(pairs))
locs = ts.individual_locations
for k, (i, j) in enumerate(pairs):
   geog_dist[k] = np.sqrt(np.sum((locs[ind_ids[i], :] - locs[ind_ids[j], :])**2))

#sanity check to see if any ind is not in the same place as itself (distance from n to n should be zero)
for (i, j), x in zip(pairs, geog_dist):
    if i == j:
        assert(x == 0)


geo= [] #empty list
for (a, b), c, d in zip(pairs, geog_dist, ind_div):   #iterate over zipped pairids, geo and ind dis, make row and append to geo list
    row = [a,b,c,d]
    geo.append(row)    
geo = pd.DataFrame(geo, columns = ['ind1', 'ind2', 'dist', 'div']) #list of lists to df

lookup = dict(zip(np.arange(0, len(ind_ids)), ind_ids)) #create lookup with indexes and corresponding ind_ids from in_ids above

geo['ind1'] = geo['ind1'].replace(lookup)
geo['ind2'] = geo['ind2'].replace(lookup)

finalgeodf = geo.merge(metadata, left_on = 'ind1', right_on = 'ind', suffixes = '_ind1').merge(metadata, left_on = 'ind2', right_on = 'ind')

finalgeodf = finalgeodf.replace(np.inf, np.nan).dropna() #REMOVE FUCKING ANNOYING INF
finalgeodf['temp'] = 0
finalgeodf.loc[(finalgeodf['time_x'] == 'now') & (finalgeodf['time_y'] == 'now'), 'temp'] = 'now'
finalgeodf.loc[(finalgeodf['time_x'] == 'then') & (finalgeodf['time_y'] == 'then'), 'temp'] = 'then'
finalgeodf = finalgeodf.loc[(finalgeodf['temp'] != 0)]
finalgeodf.to_csv('everythingisterribleafter30_geodf.csv')

(p9.ggplot(finalgeodf, p9.aes('dist', 'div', colour='temp'))
 + p9.geom_point()
 + p9.theme_bw()
 + p9.facet_wrap('~temp'))


