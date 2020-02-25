{
    volScalarField& rDeltaT = trDeltaT.ref();

    // Cache old reciprocal time scale field
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Maximum flow Courant number
    scalar maxCo(runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.5));

    // Maximum time scale
    scalar maxDeltaT(runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff(runTime.controlDict().lookupOrDefault<scalar>("rDeltaTSmoothingCoeff",0.02));

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff(runTime.controlDict().lookupOrDefault<scalar>("rDeltaTDampingCoeff",0.9));

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(runTime.controlDict().lookupOrDefault("alphaTemp",1));

    {
        // Set the reciprocal time-step from the local Courant number
        rDeltaT.ref() =
        (
            fvc::surfaceSum(amaxSf)()
           /((2*maxCo)*mesh.V())
        );
        
        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);
    }

    // Reaction source time scale
    if (alphaTemp < 1.0)
    {
        volScalarField::Internal rDeltaTT
        (
            mag(Qdot)/(alphaTemp*rho*thermo.Cp()*T)
        );

        rDeltaT.ref() = max
        (
            rDeltaT(),
            rDeltaTT
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/(rDeltaT.primitiveField() + 1.0e-12))
        << ", " << gMax(1/(rDeltaT.primitiveField() + 1.0e-12)) << endl;
}
