volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

//momentum solution without pressure gradient
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

// HbyA flux, corrected to be globally conservative and 
// ensure a solution for peqn
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

// Update the fixedFluxPressure BCs to ensure flux consistency
setSnGrad<fixedFluxPressureFvPatchScalarField>
(
    p.boundaryField(),
    (
        phiHbyA.boundaryField()
        - (mesh.Sf().boundaryField() & U.boundaryField())
    )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
);
// Non-orthogonal pressure corrector loop  
while (piso.correctNonOrthogonal())
{
   // Pressure corrector 
   fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}
#include "continuityErrs.H"

