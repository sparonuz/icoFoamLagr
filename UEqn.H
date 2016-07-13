fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U) == (1.0/rhoc)*cloudSU
);

UEqn.relax();

volScalarField rAUc(1.0/UEqn.A());
surfaceScalarField rAUcf("Dp", fvc::interpolate(rAUc));

surfaceScalarField phicForces
(
   (fvc::interpolate(rAUc*cloudVolSUSu/rhoc) & mesh.Sf())
   + rAUcf*(g & mesh.Sf())
);

if (piso.momentumPredictor())
{
    solve
    (
       UEqn
    ==
       fvc::reconstruct
       (
            phicForces/rAUcf - fvc::snGrad(p)*mesh.magSf()
       )
   );
}

