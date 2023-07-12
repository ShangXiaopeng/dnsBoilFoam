/*---------------------------------------------------------------------------*\

Description

    Creates and initialises the relative face-flux field Ji.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading/calculating face flux field phi\n" << endl;


Info<< "Reading Flux\n" << endl;

	surfaceScalarField J
	    (
		IOobject
		(
		    "J",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar ("nameJ", dimensionSet(0,3,-1,0,0,0,0), scalar(0))
	    );
	

// ************************************************************************* //
// ************************************************************************* //
