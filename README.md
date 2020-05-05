# rhoReactingCentralFoam
An OpenFOAM solver that combines rhoCentralFoam and rhoReactingFoam for high-speed reactive flows. Tested extensively for ability to model detonations. 

This solver was originally created by Caelan Lapointe at the University of Colorado Boulder. 

---
**Detonation with adaptive meshing:**

![detonation with AMR](https://github.com/duncanam/thesis/blob/master/doc/figs/amr_cells.png)



## Example Usage
Look in the example directory for an example of a detonation being modeled, with adaptive meshing (AMR), adaptive time stepping tracking both central and acoustic time stepping (see `controlDict`), and parallelized with MPI. 

Additionally, an [MS Aerospace engineering thesis](https://github.com/duncanam/thesis) with example detonation runs for static and adaptive meshing can be found [here.](https://github.com/duncanam/thesis/tree/master/sim/analysis)

---
## Notes

Feel free to create pull requests with questions or concerns. While I am not the primary author, I have rigorously tested this solver for detonation modeling for my master's degree work. I may be able to assist with usage or point you towards someone who can help. 


