surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(turbulence->muEff())*mesh.magSf()*fvc::snGrad(U)
      + (mesh.Sf() & fvc::interpolate(tauMC))
    )
    & (a_pos*U_pos + a_neg*U_neg)
);

solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  - fvc::div(sigmaDotU)
  ==
    reaction->Qdot()
);

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();
thermo.correct();
rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

if (!inviscid)
{
    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
    );
    eEqn.solve("e");

    thermo.correct();
    rhoE = rho*(e + 0.5*magSqr(U));
}

Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

p.ref() = rho()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
