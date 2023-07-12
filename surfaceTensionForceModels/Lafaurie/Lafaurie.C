/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
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

\*---------------------------------------------------------------------------*/

#include "Lafaurie.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(Lafaurie, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel, Lafaurie, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::Lafaurie::Lafaurie
(
    const word& name,
    const dictionary& surfaceTensionForceProperties,
    const incompressibleThreePhaseMixture& mixture,
    const threePhaseInterfaceProperties& interface
)
:
    surfaceTensionForceModel
    (
        name,
        surfaceTensionForceProperties,
	mixture,
	interface
    ),
    mesh_(mixture.alpha1().mesh()),
    Fstffv
    (
        IOobject
        (
            "SurfaceTensionForce",
            mixture.alpha1().time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(1,-2,-2,0,0,0,0), 0 )
    ),
    nSmoothingAlpha_
    (
        surfaceTensionForceProperties.lookupOrDefault<label>
        ("nSmoothingAlpha", 5)
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceTensionForceModels::Lafaurie::correct()
{
    //Step 1: smoothing the phase fraction field (2 passes)
    scalar CSK = 0.5;
    volScalarField alpha1s(this->mixture_.alpha1());

    for (label jj=1; jj<=nSmoothingAlpha_; jj++)
    {
	    alpha1s = 
        CSK*(fvc::average(fvc::interpolate(alpha1s))) 
	+ (1.0 - CSK)*alpha1s;
    }

    const volVectorField gradAlpha = fvc::grad(alpha1s);
    
//    const dimensionedScalar deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);

    const volVectorField ns
    (

	gradAlpha/
	(
	    mag(gradAlpha) + this->interface_.deltaN()
	)
    );

    tensor I
    (
	1, 0, 0,
	0, 1, 0,
	0, 0, 1
    );

    volVectorField Fst = 
    fvc::div
    (
	this->interface_.sigma()*mag(gradAlpha)*(I - ns*ns)
    );

    Fstffv = fvc::interpolate(Fst) & (mesh_.Sf()/mesh_.magSf());

/*
    volTensorField T = 
	this->interface_.sigma()*mag(gradAlpha)*(I - ns*ns);

    Fstffv = nHatfv & fvc::surfaceIntegrate(fvc::interpolate(T));
*/

}



bool Foam::surfaceTensionForceModels::Lafaurie::
read(const dictionary& surfaceTensionForceProperties)
{
    surfaceTensionForceModel::read(surfaceTensionForceProperties);

    return true;
}


// ************************************************************************* //
