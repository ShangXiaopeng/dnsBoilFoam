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

#include "temperatureDependentBrackbill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(temperatureDependentBrackbill, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel, temperatureDependentBrackbill, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::temperatureDependentBrackbill::temperatureDependentBrackbill
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

    TName_(surfaceTensionForceProperties.lookupOrDefault<word>("T", "T")),
    sigma_(Function1<scalar>::New("sigma", surfaceTensionForceProperties))
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceTensionForceModels::temperatureDependentBrackbill::correct()
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionSet(1, 0, -2, 0, 0)
        )
    );
    volScalarField& sigma = tsigma.ref();

    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    sigma.field() = sigma_->value(T.field());

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaBf, patchi)
    {
        sigmaBf[patchi] = sigma_->value(TBf[patchi]);
    }

    const scalar CSK = 0.5;
    const label nSmoothingAlpha(5);
//    volScalarField alpha1s(this->mixture_.alpha1());
    volScalarField alpha1s(max(this->mixture_.alpha1(), 0.0));

    for (label jj = 1; jj<=nSmoothingAlpha; jj++)
    {
	alpha1s = 
            CSK*(fvc::average(fvc::interpolate(alpha1s))) 
	  + (1.0 - CSK)*alpha1s;
    }

    volVectorField gradAlpha = fvc::grad(alpha1s);
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + this->interface_.deltaN()));
    this->interface_.correctContactAngle(nHatfv.boundaryFieldRef());

    volScalarField K = -fvc::div(nHatfv & mesh_.Sf());

//    Fstffv = fvc::interpolate(tsigma*K)*fvc::snGrad(this->mixture_.alpha1());
    Fstffv = fvc::interpolate(tsigma*K)*fvc::snGrad(alpha1s);
//    Fstffv = fvc::interpolate(tsigma*this->interface_.K())*fvc::snGrad(max(this->mixture_.alpha1(), 0.0));
//							Info << "k = " << min(K*this->interface_.nearInterface()) <<endl;
//							Info << "k = " << max(K*this->interface_.nearInterface()) <<endl;

}

bool Foam::surfaceTensionForceModels::temperatureDependentBrackbill::
read(const dictionary& surfaceTensionForceProperties)
{
    surfaceTensionForceModel::read(surfaceTensionForceProperties);

    return true;
}

/*
bool Foam::surfaceTensionForceModels::temperatureDependentBrackbill::writeData(Ostream& os) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_ << token::END_STATEMENT << nl;
        return os.good();
    }
    else
    {
        return false;
    }
}
*/
// ************************************************************************* //
