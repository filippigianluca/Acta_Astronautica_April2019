# LambTAN v1.0

`LambTAN` \(Lambert to Target Asteroids at Nodal points\) is a deterministic optimal sequence finder algorithm. It has been implemented based on lambert and to target asteroids at its nodal points. This algorithm was inspired by the combination of the Branch-And-Prune and the incremental pruning. Incremental pruning is based on the idea of constructing sets of trajectories, one arc at a time; the arcs are therefore independent, allowing to prune out full transfer subsets that do not satisfy a given criterion e.g., the ∆V of the arc is greater than a given maximum value. By constructing and assessing the arcs one after the other, the space of acceptable transfers is pruned out incrementally in such a way that the computational complexity is now polynomial with respect to the number of arcs. The solutions are composed adding branches and nodes incrementally with the iteration of the algorithm.

## Example

In oreder to see how to set and run LambTAN, check the `mainLambTAN.m` script. This is script is the one that has been used for the ATIRA publication.

## Reference

* Di Carlo, M, Ortiz, N., Romero M., J. M., Tardioli, C., Gachet, F., Kumar, K., Vasile, M., "[Optimized low-thrust mission to the atira asteroids](http://strathprints.strath.ac.uk/51194/)"

  25th AAS/AIAA Space Flight Mechanics Meeting, AAS15-299, 2015/1/11

