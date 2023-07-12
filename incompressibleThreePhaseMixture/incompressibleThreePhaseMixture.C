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

#include "incompressibleThreePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleThreePhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1_*rho1_ + alpha2_*rho2_ + alpha3_*rho3_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseMixture::incompressibleThreePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase3Name_(wordList(lookup("phases"))[2]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.db()
        ),
        U.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    nuModel3_
    (
        viscosityModel::New
        (
            "nu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),
    rho3_("rho", dimDensity, nuModel3_->viscosityProperties()),

    MW1_("MW", dimensionSet(0,0,0,0,0,0,0), subDict(phase1Name_)),
    MW2_("MW", dimensionSet(0,0,0,0,0,0,0), subDict(phase2Name_)),
    MW3_("MW", dimensionSet(0,0,0,0,0,0,0), subDict(phase3Name_)),

    a_("a", dimless, subDict(phase2Name_)),
    b_("b", dimless, subDict(phase2Name_)),
    c_("c", dimless, subDict(phase2Name_)),
    d_("d", dimless, subDict(phase2Name_)),

    Pc_("Pc", dimensionSet(1,-1,-2,0,0,0,0), subDict(phase2Name_)),
    Tc_("Tc", dimensionSet(0,0,0,1,0,0,0), subDict(phase2Name_)),

    hEvap_("hEvap", dimensionSet(0,2,-2,0,0,0,0), subDict(phase2Name_)),

    cp1_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), subDict(phase1Name_)),
    cp2_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), subDict(phase2Name_)),
    cp3_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), subDict(phase3Name_)),

    k1_("k", dimensionSet(1, 1, -3, -1, 0, 0, 0), subDict(phase1Name_)),
    k2_("k", dimensionSet(1, 1, -3, -1, 0, 0, 0), subDict(phase2Name_)),
    k3_("k", dimensionSet(1, 1, -3, -1, 0, 0, 0), subDict(phase3Name_)),

    LD_("LD", dimensionSet(0, 1, 0, 0, 0, 0, 0), subDict(phase2Name_))
{
    alpha3_ == 1.0 - alpha1_ - alpha2_;
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            alpha1_*rho1_*nuModel1_->nu()
          + alpha2_*rho2_*nuModel2_->nu()
          + alpha3_*rho3_*nuModel3_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "mu",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
          + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nu",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
              + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
            )/(alpha1f*rho1_ + alpha2f*rho2_ + alpha3f*rho3_)
        )
    );
}


bool Foam::incompressibleThreePhaseMixture::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(*this)
         && nuModel2_().read(*this)
         && nuModel3_().read(*this)
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
            nuModel3_->viscosityProperties().lookup("rho") >> rho3_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
