#include "readGravitationalAcceleration.H"

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar rhocValue
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);


volScalarField rhoc
(
    IOobject
    (
        rhocValue.name(),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhocValue
);

volScalarField muc
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoc*nu

);

Info << "Creating field alphac\n" << endl;

// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName", kinematicCloudName);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
basicKinematicTypeCloud kinematicCloud
(
    kinematicCloudName,
    rhoc,
    U,
    muc,
    g
);

// Particle fraction upper limit
scalar alphacMin
(
    1.0
  - readScalar
    (
        kinematicCloud.particleProperties().subDict("constantProperties")
       .lookup("alphaMax")
    )
);

// Update alphac from the particle locations
// alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
alphac = 1.0;
alphac.correctBoundaryConditions();

