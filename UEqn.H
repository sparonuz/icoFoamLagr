fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U) 
  == 
    (1.0/rhoc)*cloudSU
);

UEqn.relax();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

surfaceScalarField phicForces
(
   (fvc::interpolate(rAU*cloudVolSUSu/rhoc) & mesh.Sf())
   + rAUf*(g & mesh.Sf())
);

if (piso.momentumPredictor())
{
    solve
    (
       UEqn
    ==
       fvc::reconstruct
       (
            phicForces/rAUf - fvc::snGrad(p)*mesh.magSf()
       )
   );
}

