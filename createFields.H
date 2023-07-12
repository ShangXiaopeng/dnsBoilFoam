#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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

#include "createPhi.H"

#include "J.H"

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

dimensionedScalar D23("D23", dimViscosity, mixture);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);

surfaceScalarField alpha1Phi
(
    IOobject
    (
        "alpha1Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1)*phi
);

surfaceScalarField alpha2Phi
(
    IOobject
    (
        "alpha2Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha2)*phi
);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha2.name());

// #include "createMRF.H"
#include "createFvOptions.H"

/////////////////////////////////////////////////////////////////////////////

volScalarField Xv
(
    IOobject
    (
        "Xv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField mDotEvap
(
    IOobject
    (
        "mDotEvap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mDotEvap", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
);

volScalarField mDotEvapFlux
(
    IOobject
    (
        "mDotEvapFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mDotEvapFlux", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);


volScalarField gradAlpha1
(
    IOobject
    (
        "gradAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mag(fvc::grad(alpha1))
);

surfaceScalarField gradAlpha1f
(
    IOobject
    (
        "gradAlpha1f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(gradAlpha1)
);


volScalarField hDotEvap
(
    IOobject
    (
        "hDotEvap",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("hDotEvap", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
);

dimensionedScalar SMALL("SMALL", dimensionSet(0,0,0,0,0,0,0), scalar(1e-12));


dimensionedScalar mDotZero("mDotZero", dimensionSet(1,-3,-1,0,0,0,0), 0.0);

volScalarField gradXv
(
    IOobject
    (
        "gradXv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gradXv", dimensionSet(0,-1,0,0,0,0,0), scalar(0.0))
);

volScalarField rhogp
(
    IOobject
    (
        "rhogp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhogp", dimensionSet(1,-3,0,0,0,0,0), scalar(0.0))
);

volScalarField flagField
(
    IOobject
    (
        "flagField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("flagField", dimensionSet(0,0,0,0,0,0,0), scalar(0.0))
);

volScalarField liquidField
(
    IOobject
    (
        "liquidField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("liquidField", dimensionSet(0,0,0,0,0,0,0), scalar(0.0))
);


volVectorField nHat
(
    IOobject
    (
        "nHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("nHat", dimensionSet(0,0,0,0,0,0,0), vector(0.0, 0.0, 0.0))
);


volScalarField countField
(
    IOobject
    (
        "countField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("countField", dimensionSet(0,0,0,0,0,0,0), scalar(0.0))
);

Info<< "Reading transportProperties\n" << endl;
 
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
dimensionedScalar EvapThresh = transportProperties.lookup("EvapThresh");



volVectorField fSigma
(
    IOobject
    (
        "fSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("fSigma", dimensionSet(1,-2,-2,0,0,0,0), vector(0.0, 0.0, 0.0))
);

/*
volScalarField Kappa
(
    IOobject
    (
        "Kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Kappa", dimensionSet(0,-1,0,0,0,0,0), scalar(0.0))
);


volScalarField KappaS
(
    IOobject
    (
        "KappaS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("KappaS", dimensionSet(0,-1,0,0,0,0,0), scalar(0.0))
);
*/













