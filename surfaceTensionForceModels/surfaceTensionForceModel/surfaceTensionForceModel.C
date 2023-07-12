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

#include "surfaceTensionForceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceTensionForceModel, 0);
    defineRunTimeSelectionTable(surfaceTensionForceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModel::surfaceTensionForceModel
(
    const word& name,
    const dictionary& surfaceTensionForceProperties,
    const incompressibleThreePhaseMixture& mixture,
    const threePhaseInterfaceProperties& interface
)
:
    name_(name),
    surfaceTensionForceProperties_(surfaceTensionForceProperties),
    mixture_(mixture),
    interface_(interface),
    Fstffv
    (
        IOobject
        (
            "SurfaceTensionForce",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.alpha1().mesh(),
        dimensionedScalar( "dummy", dimensionSet(1,-2,-2,0,0,0,0), 0 )
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::surfaceTensionForceModel::pcap() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "DummyPC",
                mixture_.alpha1().time().timeName(),
                mixture_.alpha1().mesh()
            ),
            mixture_.alpha1().mesh(),
            dimensionedScalar("0", dimensionSet(1, -1, -2, 0, 0), 0)
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::surfaceTensionForceModel::phi_c(const surfaceScalarField& rAUf_) const
{
    return 
    tmp<surfaceScalarField>( this->Fstff() * rAUf_ * mixture_.alpha1().mesh().magSf() );
}


bool Foam::surfaceTensionForceModel::read
(
    const dictionary& surfaceTensionForceProperties
)
{
    surfaceTensionForceProperties_ = surfaceTensionForceProperties;
    
    return true;
}

// ************************************************************************* //
