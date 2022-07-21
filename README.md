# move_review
code supporting MOVE Review

Simulation implementing (crudely):

* Spatial map - two nice habitats connected by a corridor of kind-of-nice habitats

* Temporal landscape change - after 2/3 generations have elapsed, the corridor will disappear

* Interventions - either at a specified generation (e.g. after landscape change) or every nth generation, individual fitness can be changed. Say, implement a mapwide reduction of fitness to represent global spraying, or targeted reduction in an area to represent localised population suppression

* Changes in dispersal (σ)


- Simulations (in SLiM) are contained in *sims*
- R and Python Scripts for parsing the output and performing various analyses are contained in *bin*
- Treesequences (cool!) and other sim outputs are contained in sim_outputs
