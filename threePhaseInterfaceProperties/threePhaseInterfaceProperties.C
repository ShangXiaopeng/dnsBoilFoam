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

\*---------------------------------------------------------------------------*/

#include "threePhaseInterfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

//#include "surfaceTensionForceModel.H"
//#include "fvcSmooth.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::threePhaseInterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::threePhaseInterfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb
//    const surfaceVectorField::Boundary& gradAlphaf
) const
{

//////////////////////////////////////////////////////////////////////////////////////////////

    const volScalarField::Boundary& alpha1 =
        mixture_.alpha1().boundaryField();
    const volScalarField::Boundary& alpha2 =
        mixture_.alpha2().boundaryField();
    const volScalarField::Boundary& alpha3 =
        mixture_.alpha3().boundaryField();
    const volVectorField::Boundary& U =
        mixture_.U().boundaryField();

    const fvMesh& mesh = mixture_.U().mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& a2cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha2[patchi]);

            const alphaContactAngleFvPatchScalarField& a3cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha3[patchi]);

            scalarField twoPhaseAlpha2(max(a2cap, scalar(0)));
            scalarField twoPhaseAlpha3(max(a3cap, scalar(0)));

            scalarField sumTwoPhaseAlpha
            (
                twoPhaseAlpha2 + twoPhaseAlpha3 + small
            );

            twoPhaseAlpha2 /= sumTwoPhaseAlpha;
            twoPhaseAlpha3 /= sumTwoPhaseAlpha;

            fvsPatchVectorField& nHatp = nHatb[patchi];

            scalarField theta
            (
                convertToRad
              * (
                   twoPhaseAlpha2*(180 - a2cap.theta(U[patchi], nHatp))
                 + twoPhaseAlpha3*(180 - a3cap.theta(U[patchi], nHatp))
                )
            );

            vectorField nf(boundary[patchi].nf());


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*

    const volScalarField::Boundary& alpha1 =
        mixture_.alpha1().boundaryField();
    const volScalarField::Boundary& alpha2 =
        mixture_.alpha2().boundaryField();
    const volScalarField::Boundary& alpha3 =
        mixture_.alpha3().boundaryField();
    const volVectorField::Boundary& U =
        mixture_.U().boundaryField();

    const fvMesh& mesh = mixture_.U().mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& a2cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha2[patchi]);

            const alphaContactAngleFvPatchScalarField& a3cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha3[patchi]);

            scalarField twoPhaseAlpha2(max(a2cap, scalar(0)));
            scalarField twoPhaseAlpha3(max(a3cap, scalar(0)));

            scalarField sumTwoPhaseAlpha
            (
                twoPhaseAlpha2 + twoPhaseAlpha3 + small
            );


            twoPhaseAlpha2 /= sumTwoPhaseAlpha;
            twoPhaseAlpha3 /= sumTwoPhaseAlpha;

            fvsPatchVectorField& nHatp = nHatb[patchi];

            vectorField nf(boundary[patchi].nf());

//------------------------------------------------------------------------------//

	    const scalar ks0(4.5276e10);
	    const scalar h(6.626e-34);
	    const scalar Vm(3e-29);
	    const scalar l(0.5e-9);
	    const scalar n(Foam::pow(l, -2));
	    const scalar kB(1.38e-23);
	    const scalarField sigma(sigma12_.value()*twoPhaseAlpha2 + sigma13_.value()*twoPhaseAlpha3);

	    const dimensionedScalar nu
	    (
		"nu",
		dimensionSet(1, 1, -3, -1, 0, 0, 0),
		mixture_.subDict(mixture_.phase1Name()).lookup("nu")
	    );

	    scalarField mu
	    (
		mixture_.mu()().boundaryField()[patchi]

	    );
//----------------------------------------------------------------------------//

//	    const fvPatchVectorField& Uwall = U[patchi];
//	    vectorField Up(Uwall.patchInternalField() - Uwall);
//
//	    scalarField Ucl
//	    (
//		(-Up & nHatp)/sqrt(1.0 - (nHatp & nf))
//	    );




	    const fvPatchVectorField& Uwall = U[patchi];
	    vectorField Up(Uwall.patchInternalField() - Uwall);
	    Up -= (nf & Up)*nf;

	    vectorField nWall(nHatp - (nf & nHatp)*nf);
	    nWall /= (mag(nWall) + small);

	    scalarField Ucl((twoPhaseAlpha2 + twoPhaseAlpha3)*(-nWall & Up));

//---------------------------------------------------------------------------//

	    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
	    const fvPatchScalarField& Twall = T.boundaryField()[patchi];
	    scalarField Tp(Twall.patchInternalField());

            scalarField theta0
            (
                convertToRad
              * (
                   twoPhaseAlpha2*(180 - a2cap.theta(U[patchi], nHatp))
                 + twoPhaseAlpha3*(180 - a3cap.theta(U[patchi], nHatp))

//                   twoPhaseAlpha2*(180 - theta0_.value())
//                 + twoPhaseAlpha3*(180 - theta0_.value())
                )
            );

	    scalarField K1
	    (
		2.*ks0*l*h/mu/Vm
	    );

	    const scalarField K2(sigma/2./n/kB);

	    scalarField theta
	    (
	        Foam::acos
		(
		    min
		    (
			max(cos(theta0) - Foam::asinh(Ucl/K1)*Tp/(K2 + small), -1.0),
//			1.0
			0.9998*pos(twoPhaseAlpha2 + twoPhaseAlpha3) + 1.0*neg0(twoPhaseAlpha2 + twoPhaseAlpha3)
		    )
		)
	    );
//theta.write();
//	    theta = max(min(theta, 1.5707), -1.5707)

Info << "min(theta) = " << min(theta) <<endl;
Info << "max(theta) = " << max(theta) <<endl;


*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    const volScalarField::Boundary& alpha1 =
        mixture_.alpha1().boundaryField();
    const volScalarField::Boundary& alpha2 =
        mixture_.alpha2().boundaryField();
    const volScalarField::Boundary& alpha3 =
        mixture_.alpha3().boundaryField();
    const volVectorField::Boundary& U =
        mixture_.U().boundaryField();

    const fvMesh& mesh = mixture_.U().mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& a2cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha2[patchi]);

            const alphaContactAngleFvPatchScalarField& a3cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha3[patchi]);

            scalarField twoPhaseAlpha2(max(a2cap, scalar(0)));
            scalarField twoPhaseAlpha3(max(a3cap, scalar(0)));

            scalarField sumTwoPhaseAlpha
            (
                twoPhaseAlpha2 + twoPhaseAlpha3 + small
            );


            twoPhaseAlpha2 /= sumTwoPhaseAlpha;
            twoPhaseAlpha3 /= sumTwoPhaseAlpha;

            fvsPatchVectorField& nHatp = nHatb[patchi];

            vectorField nf(boundary[patchi].nf());

//------------------------------------------------------------------------------//

	    const scalar sigma(sigma12_.value()/2 + sigma13_.value()/2);

	    const dimensionedScalar nu
	    (
		"nu",
		dimensionSet(1, 1, -3, -1, 0, 0, 0),
		mixture_.subDict(mixture_.phase1Name()).lookup("nu")
	    );

	    scalarField mu
	    (
		mixture_.mu()().boundaryField()[patchi]

	    );
//----------------------------------------------------------------------------//

	    const fvPatchVectorField& Uwall = U[patchi];
	    vectorField Up(Uwall.patchInternalField() - Uwall);
	    Up -= (nf & Up)*nf;

	    vectorField nWall(nHatp - (nf & nHatp)*nf);
	    nWall /= (mag(nWall) + small);

	    scalarField Ucl(-nWall & Up);

	    scalarField Ca( mu*Ucl/sigma );

	    scalarField x
	    (
		max
		(
		    Ca + fInverse_.value(),
		    0.0
		)
	    );

//---------------------------------------------------------------------------//
            scalarField theta
	    (
		Foam::acos
		(
		    min
		    (
			max
			(
			    1.0 - 2.0*Foam::tanh
			    (
				5.16*Foam::pow
				(
				    x/(1.0 + 1.31*Foam::pow(x, 0.99)), 0.706
				)
			    ),
			    -1.0
			),
			1.0
		    )
		)
	    );


Info << "min(theta) = " << min(theta) <<endl;
Info << "max(theta) = " << max(theta) <<endl;

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
	theta = max
	(
	    min
	    (
		theta, 3.054326
	    ),
	    0.087266
	);
*/

            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatp & nf);

            scalarField b1(cos(theta));

            scalarField b2(nHatp.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;

            nHatp /= (mag(nHatp) + deltaN_.value());


//////////////////////////////////////////////////////////////////////////////////////
/*
            alphaContactAngleFvPatchScalarField& a1cap =
//                refCast<const alphaContactAngleFvPatchScalarField>
//                (alpha1[patchi]);
	    const_cast<alphaContactAngleFvPatchScalarField&>
	    (
		refCast<const alphaContactAngleFvPatchScalarField>
		(
		    alpha1[patchi]
		)
	    );

	    const volVectorField gradAlpha(fvc::grad(mixture_.alpha1(), "nHat"));
	    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
            a1cap.gradient() = (nf & nHatp)*mag(gradAlphaf.boundaryField()[patchi]);
            a1cap.evaluate();
*/
//////////////////////////////////////////////////////////////////////////////////////////
        }
    }
}


void Foam::threePhaseInterfaceProperties::calculateK()
{
    const volScalarField& alpha1 = mixture_.alpha1();

    const fvMesh& mesh = alpha1.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha1));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryFieldRef());
//    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    // volVectorField nHat = gradAlpha/(mag(gradAlpha) + deltaN_);
    // nHat.boundaryField() = nHatfv.boundaryField();
    // K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);

//-----------------------------------------------------//
    const scalar CSK = 0.5;
    const label nSmoothingAlpha(5);
    volScalarField alpha1s(alpha1);

    for (label jj = 1; jj<=nSmoothingAlpha; jj++)
    {
	alpha1s = 
            CSK*(fvc::average(fvc::interpolate(alpha1s))) 
	  + (1.0 - CSK)*alpha1s;
    }

    volVectorField gradAlphaSmooth = fvc::grad(alpha1s);
    surfaceVectorField gradAlphafSmooth(fvc::interpolate(gradAlphaSmooth));
    surfaceVectorField nHatfvSmooth(gradAlphafSmooth/(mag(gradAlphafSmooth) + deltaN_));
    correctContactAngle(nHatfvSmooth.boundaryFieldRef());
    nHatfSmooth_ = nHatfvSmooth & Sf;
//--------------------------------------------------------//

//    nHat_.boundaryField() = nHatfv.boundaryField();
    stf_->correct();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threePhaseInterfaceProperties::threePhaseInterfaceProperties
(
    const incompressibleThreePhaseMixture& mixture
)
:
    mixture_(mixture),
    cAlpha_
    (
        readScalar
        (
            mixture.U().mesh().solverDict
            (
                mixture_.alpha1().name()
            ).lookup("cAlpha")
        )
    ),
    cAlpha2_
    (
        readScalar
        (
            mixture.U().mesh().solverDict
            (
                mixture_.alpha1().name()
            ).lookup("cAlpha2")
        )
    ),
    cAlpha3_
    (
        readScalar
        (
            mixture.U().mesh().solverDict
            (
                mixture_.alpha1().name()
            ).lookup("cAlpha3")
        )
    ),

//    theta0_("theta0", dimensionSet(0, 0, 0, 0, 0), mixture),

    sigma12_("sigma12", dimensionSet(1, 0, -2, 0, 0), mixture),
    sigma13_("sigma13", dimensionSet(1, 0, -2, 0, 0), mixture),

    fInverse_("fInverse", dimensionSet(0, 0, 0, 0, 0), mixture),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mixture.U().mesh().V()), 1.0/3.0)
    ),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh()
        ),
        mixture.alpha1().mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    nHatfSmooth_
    (
        IOobject
        (
            "nHatfSmooth",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh()
        ),
        mixture.alpha1().mesh(),
        dimensionedScalar("nHatfSmooth", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh()
        ),
        mixture.alpha1().mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{

    stf_ = 
	surfaceTensionForceModel::New
        (
	    mixture.subDict("surfaceTensionModel").lookup("model"),
	    mixture.subDict("surfaceTensionModel"),
	    mixture,
	    *this
        );

    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::threePhaseInterfaceProperties::surfaceTensionForce() const
{
//    return fvc::interpolate(sigmaK())*fvc::snGrad(mixture_.alpha1());
    return stf_->Fstff();
}


Foam::tmp<Foam::volScalarField>
Foam::threePhaseInterfaceProperties::nearInterface() const
{
    return max
    (
        pos0(mixture_.alpha1() - 0.01)*pos0(0.99 - mixture_.alpha1()),
        pos0(mixture_.alpha2() - 0.01)*pos0(0.99 - mixture_.alpha2())
    );
}

Foam::scalar Foam::threePhaseInterfaceProperties::HoffmanFunction
(
    const scalar& x
) const
{
    return acos(1 - 2*tanh(5.16*pow(x/(1+1.31*pow(x, 0.99)), 0.706)));
}
// ************************************************************************* //
