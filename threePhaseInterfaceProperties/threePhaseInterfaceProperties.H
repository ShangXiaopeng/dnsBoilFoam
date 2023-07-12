/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::threePhaseInterfaceProperties

Description
    Properties to aid interFoam :
    1. Correct the alpha boundary condition for dynamic contact angle.
    2. Calculate interface curvature.

SourceFiles
    threePhaseInterfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef threePhaseInterfaceProperties_H
#define threePhaseInterfaceProperties_H

#include "incompressibleThreePhaseMixture.H"
#include "surfaceFields.H"
#include "surfaceTensionForceModel.H"

#include "RiddersRoot.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//class surfaceTensionForceModel;

/*---------------------------------------------------------------------------*\
              Class threePhaseInterfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class threePhaseInterfaceProperties
{
    // Private data

        const incompressibleThreePhaseMixture& mixture_;

        //- Compression coefficient
        scalar cAlpha_;
        scalar cAlpha2_;
        scalar cAlpha3_;

//        dimensionedScalar theta0_;

        //- Surface tension 1-2
        dimensionedScalar sigma12_;

        //- Surface tension 1-3
        dimensionedScalar sigma13_;

        //- f inverse
        dimensionedScalar fInverse_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        surfaceScalarField nHatf_;
        surfaceScalarField nHatfSmooth_;
        volScalarField K_;

autoPtr<surfaceTensionForceModel> stf_;

//surfaceScalarField surfaceTension_;///////////////////////////////////////////////////////////////////////////////

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        threePhaseInterfaceProperties(const threePhaseInterfaceProperties&);
        void operator=(const threePhaseInterfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle.
        //  Calculated from the component of U parallel to the wall
//        void correctContactAngle
//        (
//            surfaceVectorField::Boundary& nHat
//        ) const;

        //- Re-calculate the interface curvature
        void calculateK();


//-----------------------------------------------------------------------------------------//


public:

    //eb - Function object class, inverse of the Hoffman function
    class InverseHoffmanFunction : public FuncBase
    {

        // Private data

            //- Equilibrium contact angle (advancing or receding), as parameter
            scalar thetaE_;

    public:

        // Constructors

            //- Construct from data
            InverseHoffmanFunction(const scalar& thetaE)
            :
                thetaE_(thetaE)
            {}

        // operator()
            scalar operator()(scalar fHI) const
            {
                return
                  pow(fHI,0.706)
                -   pow((1+1.31*pow(fHI,0.99)),0.706)
                   *(1.0/5.16)*atanh((1-cos(thetaE_))/2);
            }
    };


    //eb - Evaluate the value of the Hoffman function
    scalar HoffmanFunction
    (
        const scalar& x
    ) const;

//-----------------------------------------------------------------------------------------//


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field alpha and IOdictionary
        threePhaseInterfaceProperties
        (
            const incompressibleThreePhaseMixture& mixture
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        scalar cAlpha2() const
        {
            return cAlpha2_;
        }

        scalar cAlpha3() const
        {
            return cAlpha3_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        surfaceScalarField& nHatfSmooth()
        {
            return nHatfSmooth_;
        }

        const volScalarField& K() const
        {
            return K_;
        }

        tmp<volScalarField> sigma() const
        {

            volScalarField limitedAlpha2(max(mixture_.alpha2(), scalar(0)));
            volScalarField limitedAlpha3(max(mixture_.alpha3(), scalar(0)));

            return
                (limitedAlpha2*sigma12_ + limitedAlpha3*sigma13_)
               /(limitedAlpha2 + limitedAlpha3 + small);

/*
            volScalarField limitedAlpha1(max(mixture_.alpha1(), scalar(0)));
            volScalarField limitedAlpha2(max(mixture_.alpha2(), scalar(0)));
            volScalarField limitedAlpha3(max(mixture_.alpha3(), scalar(0)));

	    volScalarField rho
	    (
		mixture_.rho1()*limitedAlpha1 + mixture_.rho2()*limitedAlpha2 + mixture_.rho3()*limitedAlpha3
	    );

	    dimensionedScalar rhog
	    (
		(mixture_.rho2() + mixture_.rho3())/2.0
	    );

            return
                0.5*(sigma12_ + sigma13_)*2.0*rho
               /(mixture_.rho1() + rhog);
*/
        }

        tmp<volScalarField> sigmaK() const
        {
            return sigma()*K_;
        }

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        void correct()
        {
            calculateK();
        }
/*
        tmp<surfaceScalarField> 
        phi_c(const surfaceScalarField& rAUf) const
	{
	    return stf_->phi_c(rAUf);
	}
*/
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat
//    const surfaceVectorField::Boundary& gradAlphaf
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
