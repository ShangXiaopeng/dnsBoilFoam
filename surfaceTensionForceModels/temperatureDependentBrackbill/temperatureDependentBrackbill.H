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

Class
    Foam::surfaceTensionForceModel::temperatureDependentBrackbill

Description
    Default surface tension force model, from temperatureDependentBrackbill, J. U., Kothe, D. B.,
    & Zemach, C. (1992). A continuum method for modeling surface tension.
    Journal of computational physics, 100(2), 335-354.

SourceFiles
    temperatureDependentBrackbill.C

\*--------------------------------------------------------------------*/

#ifndef temperatureDependentBrackbill_H
#define temperatureDependentBrackbill_H

#include "surfaceTensionForceModel.H"
#include "fvCFD.H"
#include "threePhaseInterfaceProperties.H"

#include "Function1.H"
// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{

/*--------------------------------------------------------------------*\
                             class temperatureDependentBrackbill
\*--------------------------------------------------------------------*/

class temperatureDependentBrackbill
:
    public surfaceTensionForceModel
{

// Private Member Functions

    //- This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //- Surface tension forces on cell faces
    surfaceScalarField Fstffv;


word TName_;
autoPtr<Function1<scalar>> sigma_;


public:

    //- Runtime type information
    TypeName("temperatureDependentBrackbill");


    // Constructors

    //- construct from components
    temperatureDependentBrackbill
    (
        const word& name,
        const dictionary& surfaceTensionForceProperties,
        const incompressibleThreePhaseMixture& mixture,
	const threePhaseInterfaceProperties& interface
    );



    //- Destructor
    ~temperatureDependentBrackbill()
    {}


    // Member Functions

        //- Correct the surface tension model
        virtual void correct();
        
        //- Returns surface tension force
        virtual tmp<surfaceScalarField> Fstff() const
        {
            return Fstffv;
        }


        //- Read the transportProperties dictionary and update
        bool read(const dictionary& surfaceTensionForceProperties);

        //- Write in dictionary format
//        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionForceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
