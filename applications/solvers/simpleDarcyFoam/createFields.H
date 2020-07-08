//- Read pressure and velocity
Info << nl << "Reading field p\n" << endl;
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
mesh.setFluxRequired(p.name());

//- Read fluid properties
Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


const dimensionedScalar rho(transportProperties.lookup("rho"));
const dimensionedScalar  mu(transportProperties.lookup("mu"));;

//#include "createPhi.H"

#include "createFvOptions.H"

//- Read permeability
Info << "Reading field K\n" << endl;
volTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
K.write();


//- Create mobility surface tensor field
surfaceTensorField Mf("Mf",fvc::interpolate(K,"K")/mu);

//- Create gravitational flux
surfaceScalarField phiG("phiG",(rho * (Mf & g)) & mesh.Sf());

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(phiG)
);
