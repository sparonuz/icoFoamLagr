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

    #include "createControl.H"

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
                Zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = Zero;

        #include "UEqn.H"

        // --- PISO loop
        while (piso.correct())
        {
            #include "pEqn.H"
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
