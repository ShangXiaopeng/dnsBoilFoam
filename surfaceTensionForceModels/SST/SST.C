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

#include "SST.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(SST, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel, SST, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::SST::SST
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
    Cpc( surfaceTensionForceProperties.lookupOrDefault<scalar>("Cpc", 0.5) ),
    phi_c_thresholdFactor
    (
        surfaceTensionForceProperties.lookupOrDefault<scalar>
        ("ThresholdFactor", 0.01) 
    ),
    Cfc_
    (
        surfaceTensionForceProperties.lookupOrDefault<scalar>
        ("Cfc", 0.1)
    ),
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
    //Sharp surface force (capillary forces on cell faces)
    fcf
    (
        IOobject
        (
            "fcf",
            mixture.alpha1().time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        ("fcf0", dimMass/(dimLength*dimLength*dimTime*dimTime), 0)
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceTensionForceModels::SST::correct()
{


    //Step 1: smoothing the phase fraction field (2 passes)
    scalar CSK = 0.5;
    //pass 1
    volScalarField alpha1s =
          CSK*
	  (
	      fvc::average(fvc::interpolate(this->mixture_.alpha1()))
	  ) 
        + (1.0 - CSK) * this->mixture_.alpha1();
        
    //pass 2
    alpha1s = 
        CSK*
	(
	    fvc::average(fvc::interpolate(alpha1s))
	) 
	+ (1.0 - CSK)*alpha1s;


    const volScalarField& alpha(alpha1s);


    //Step 2: initialize interface curvature Kappa
    const volVectorField gradAlpha = fvc::grad(alpha1s);
    
    const dimensionedScalar 
    deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);
    
    const volVectorField ns(gradAlpha/(mag(gradAlpha) + deltaN));

    volScalarField K = fvc::div(ns);

    //Step 3: smooth curvature field (2 passes)
    volScalarField w = 
	Foam::sqrt
	(
	      max( alpha*(1.0 - alpha), 0.) 
	    + 1.0E-6
	);

    volScalarField Coeff = 2.0*Foam::sqrt
	(
	    max
	    (
		alpha*(1.0 - alpha), 0.
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

    //Step 4: compute smoothed curvature on faces
    surfaceScalarField Kf = fvc::interpolate(w*Ks)/fvc::interpolate(w); 

    //Step 5: compute interface delta function from sharpened interface field
    volScalarField alpha1_pc =
        1.0/(1.0 - Cpc)*
	(
	      min
	      (
//		 max(alpha, Cpc/2.0), (1.0 - Cpc/2.0)
		 max(alpha, 0.), (1.0 - Cpc/2.0)
	      ) 
	    - Cpc/2.0
	);

    surfaceScalarField deltasf = fvc::snGrad(alpha1_pc);

    //Step 6: compute surface tension force on faces

    fcf = -fvc::interpolate(this->interface_.sigma())*Kf*deltasf;

    Fstffv = fcf;
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::surfaceTensionForceModels::SST::
phi_c(const surfaceScalarField& rAUf_) const
{
    const surfaceScalarField phi_c_i( Fstffv * rAUf_ * mesh_.magSf() );

    //Apply limiting (filtering)
    const dimensionedScalar
    dummyFlux("dummyFlux", dimensionSet(0,3,-1,0,0,0,0), 1.0);

    dimensionedScalar phi_c_thresh
    ( 
          phi_c_thresholdFactor*dummyFlux*gSum( mag(phi_c_i.field()) ) 
        / gSum( SMALL + pos( mag(phi_c_i.field()) - SMALL ) ) 
    );

    //Return filtered phi_c
    return tmp<surfaceScalarField>
    ( 
	phi_c_i - max( min(phi_c_i, phi_c_thresh), -phi_c_thresh )
    );
}

bool Foam::surfaceTensionForceModels::SST::
read(const dictionary& surfaceTensionForceProperties)
{
    surfaceTensionForceModel::read(surfaceTensionForceProperties);

    return true;
}


// ************************************************************************* //
