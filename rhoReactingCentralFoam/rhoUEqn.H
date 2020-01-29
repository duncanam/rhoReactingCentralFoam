solve
(
    fvm::ddt(rhoU)
  + fvc::div(phiUp)
);

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

if (!inviscid)
{
    solve
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
      - fvm::laplacian(turbulence->muEff(), U)
      - fvc::div(tauMC)
    );
    rhoU = rho*U;
}

Info << "\nmin / max mag(U) ; " << gMin(mag(U)()()) << " / " << gMax(mag(U)()()) << nl << endl;
