# Simple Detonation Example

1. copy the `0.orig` directory into a new `0` directory
2. run `blockMesh`
3. run either:
	* `funkySetFields -time 0` for gradient ignition region (preferred method)
	* `setFields` for block detonation ignition region
4. run `decomposePar`, and make sure to edit the file to reflect the number of processors you are capable of using (default 20)
5. run `rhoReactingCentralFoam`

To post-process, you'll need to:
1. `reconstructParMesh`, add `-latestTime` flag for only last result
2. `reconstructPar`, add `-latestTime` flag for only last result
3. `postProcess -func sample` to use the `sample` line sampling utility
	
