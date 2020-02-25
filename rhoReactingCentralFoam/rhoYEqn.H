tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            //volScalarField& rhoYi = rhoY[i];
	    volScalarField& Yi = Y[i];

	    //surfaceScalarField Yi_pos = interpolate(Yi, pos, "reconstruct(Yi)");
	    //surfaceScalarField Yi_neg = interpolate(Yi, neg, "reconstruct(Yi)");
            //surfaceScalarField phiYi("phiYi", aphiv_pos*rho_pos*Yi_pos + aphiv_neg*rho_neg*Yi_neg);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
	      - reaction->R(Yi)
            );
	    
            if (!inviscid)
	    {
                YiEqn -= fvm::laplacian(turbulence->muEff(), Yi);
            }

	    YiEqn.solve("Yi");
            Yi.max(0.0);
	    Yi.correctBoundaryConditions();
            Yt += Yi;
        }
    }
    
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
