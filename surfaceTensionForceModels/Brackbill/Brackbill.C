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

#include "Brackbill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(Brackbill, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel, Brackbill, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::Brackbill::Brackbill
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

void Foam::surfaceTensionForceModels::Brackbill::correct()
{

    const scalar CSK = 0.5;
    const label nSmoothingAlpha(5);
    volScalarField alpha1s(this->mixture_.alpha1());

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
//    Fstffv = fvc::interpolate(this->interface_.sigma()*K)*fvc::snGrad(this->mixture_.alpha1());
//    Fstffv = fvc::interpolate(this->interface_.sigma()*K)*fvc::snGrad(alpha1s);

////////////////////////////////////////////////////////////////////////
    volScalarField w = 
	Foam::sqrt
	(
	      max( this->mixture_.alpha1()*(1.0 - this->mixture_.alpha1()), 0.) 
	    + 1.0E-6
	);

    volScalarField Coeff = 2.0*Foam::sqrt
	(
	    max
	    (
		this->mixture_.alpha1()*(1.0 - this->mixture_.alpha1()), 0.
	    )
	);

    //pass 1
    volScalarField Ks_star = 
        fvc::average(fvc::interpolate(K*w))
	/fvc::average(fvc::interpolate(w));

    volScalarField Ks = Coeff*K + (1.0 - Coeff)*Ks_star;

    //pass 2
    Ks_star = 
        fvc::average(fvc::interpolate(Ks*w))
	/fvc::average(fvc::interpolate(w));

    Ks = Coeff*K + (1.0 - Coeff)*Ks_star;

    surfaceScalarField Kf = fvc::interpolate(w*Ks)/fvc::interpolate(w); 
//    Fstffv = fvc::interpolate(this->interface_.sigma())*Kf*fvc::snGrad(alpha1s);
    Fstffv = fvc::interpolate(this->interface_.sigma())*Kf*fvc::snGrad(this->mixture_.alpha1());
////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------------------------//
/*
//    volVectorField ns(gradAlpha/(mag(gradAlpha) + this->interface_.deltaN()));
//    volScalarField K = fvc::div(ns);
    volScalarField K = fvc::div(nHatfv & mesh_.Sf());

    volScalarField w = 
	Foam::sqrt
	(
	      max( alpha1s*(1.0 - alpha1s), 0.) 
	    + 1.0E-6
	);

    volScalarField Coeff = 2.0*Foam::sqrt
	(
	    max
	    (
		alpha1s*(1.0 - alpha1s), 0.
	    )
	);

    //pass 1
    volScalarField Ks_star = 
        fvc::average(fvc::interpolate(K*w))
	/fvc::average(fvc::interpolate(w));

    volScalarField Ks = Coeff*K + (1.0 - Coeff)*Ks_star;

    //pass 2
    Ks_star = 
        fvc::average(fvc::interpolate(Ks*w))
	/fvc::average(fvc::interpolate(w));

    Ks = Coeff*K + (1.0 - Coeff)*Ks_star;

    surfaceScalarField Kf = fvc::interpolate(w*Ks)/fvc::interpolate(w); 


    const scalar Cpc = 0.5;
    volScalarField alpha1_pc =
        1.0/(1.0 - Cpc)*
	(
	      min
	      (
		 max(this->mixture_.alpha1(), Cpc/2.0), (1.0 - Cpc/2.0)
//		 max(alpha, 0.), (1.0 - Cpc/2.0)
	      ) 
	    - Cpc/2.0
	);
    surfaceScalarField deltasf = fvc::snGrad(alpha1_pc);

//    surfaceScalarField deltasf = fvc::snGrad(alpha1s);

    Fstffv = fvc::interpolate(this->interface_.sigma()*K)*deltasf;
*/
//-----------------------------------------------------------------------------------------------//


//    Fstffv = fvc::interpolate(this->interface_.sigmaK())*fvc::snGrad(this->mixture_.alpha1());

//							Info << "k = " << min(this->interface_.sigmaK()*this->interface_.nearInterface()) <<endl;
//							Info << "k = " << max(this->interface_.sigmaK()*this->interface_.nearInterface()) <<endl;

}

bool Foam::surfaceTensionForceModels::Brackbill::
read(const dictionary& surfaceTensionForceProperties)
{
    surfaceTensionForceModel::read(surfaceTensionForceProperties);

    return true;
}


// ************************************************************************* //
