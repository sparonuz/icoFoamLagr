/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "pisoControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#ifdef MPPIC
    #include "basicKinematicMPPICCloud.H"
    #define basicKinematicTypeCloud basicKinematicMPPICCloud
#else
    #include "basicKinematicCollidingCloud.H"
    #define basicKinematicTypeCloud basicKinematicCollidingCloud
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
       "cloudName",
       "name",
       "specify alternative cloud name. default is 'kinematicCloud'"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //include presence of cloud of particles, and evolve
        continuousPhaseTransport.correct();
        muc = rhoc*continuousPhaseTransport.nu();

        Info<< "Evolving " << kinematicCloud.name() << endl;
        kinematicCloud.evolve();

        fvVectorMatrix cloudSU(kinematicCloud.SU(U));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector
            (
                "0",
                cloudSU.dimensions()/dimVolume,
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.internalField() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = vector::zero;

        // correct UEqn taking into account particles
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

        Info<< "Writing U matrix " << endl;
        if (piso.momentumPredictor())
        {
            Info<< "Solving U Eqn " << endl;
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

         Info<< "Entering PISO loop " << endl;
        // --- PISO loop
        while (piso.correct())
        {
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
//            adjustPhi(phiHbyA, U, p);
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
                Info<< "Writing P matrix " << endl;
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                Info<< "SOlving P Equation " << endl;
                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"
            //Correct the approximate velocity field 
            //using the corrected pressure gradient
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
