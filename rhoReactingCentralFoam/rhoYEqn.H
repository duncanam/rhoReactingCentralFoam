{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& rhoYi = rhoY[i];
	        volScalarField& Yi = Y[i];

	        surfaceScalarField Yi_pos = fvc::interpolate(Yi, pos, "reconstruct(Yi)");
	        surfaceScalarField Yi_neg = fvc::interpolate(Yi, neg, "reconstruct(Yi)");
	        surfaceScalarField phiYi("phiYi", aphiv_pos*rho_pos*Yi_pos + aphiv_neg*rho_neg*Yi_neg);

	        solve(fvm::ddt(rhoYi) + fvc::div(phiYi));

	        Yi.ref() = rhoYi()/rho();
	        Yi.correctBoundaryConditions();
	        Yi.max(0.0);
	        rhoYi.boundaryFieldRef() == rho.boundaryField()*Yi.boundaryField();

            if (!inviscid)
	    {
                solve
                (
                     fvm::ddt(rho,Yi) - fvc::ddt(rho,Yi)
                   + fvm::div(phi,Yi) - fvc::div(phi,Yi)
                   - fvm::laplacian(turbulence->muEff(), Yi)
                    ==
                      reaction->R(Yi)
                );
                rhoYi = rho*Yi;
            }

            Yi.correctBoundaryConditions();
            Yi.max(0.0);
            Yt += Yi;
        }
    }
    
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
