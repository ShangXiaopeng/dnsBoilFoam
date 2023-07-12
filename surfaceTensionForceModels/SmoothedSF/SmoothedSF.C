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

#include "SmoothedSF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(SmoothedSF, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionForceModel,
        SmoothedSF,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::SmoothedSF::SmoothedSF
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
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceTensionForceModels::SmoothedSF::correct()
{
    
    //Smoothed fields
    volScalarField sigmaK_sm = this->interface_.sigmaK();
    volScalarField alpha1_sm = this->mixture_.alpha1();

    //Smoothing operations
    sigmaK_sm = fvc::average( fvc::interpolate(sigmaK_sm) );
    alpha1_sm = fvc::average( fvc::interpolate(alpha1_sm) );


    Fstffv = fvc::interpolate( sigmaK_sm )*fvc::snGrad(alpha1_sm);

}

bool Foam::surfaceTensionForceModels::SmoothedSF::
read(const dictionary& surfaceTensionForceProperties)
{
    surfaceTensionForceModel::read(surfaceTensionForceProperties);

    return true;
}


// ************************************************************************* //
