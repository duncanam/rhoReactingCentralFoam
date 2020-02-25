bool useAcousticCo = runTime.controlDict().lookupOrDefault<Switch>("useAcousticCo",false);

scalar maxAcousticCo
(
    readScalar(runTime.controlDict().lookup("maxAcousticCo"))
);

#include "createTimeControls.H"
