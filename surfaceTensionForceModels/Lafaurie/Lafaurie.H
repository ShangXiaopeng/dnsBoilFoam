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
    Foam::surfaceTensionForceModel::Lafaurie

Description
    Sharp surface tension force model, from Raeini, A. Q., Blunt, M. J., & 
    Bijeljic, B. (2012). Modelling two-phase flow in porous media at the pore 
    scale using the volume-of-fluid method. Journal of Computational Physics, 
    231(17), 5653-5668.


SourceFiles
    Lafaurie.C

\*--------------------------------------------------------------------*/

#ifndef Lafaurie_H
#define Lafaurie_H

#include "surfaceTensionForceModel.H"
#include "fvCFD.H"
#include "threePhaseInterfaceProperties.H"
// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{

/*--------------------------------------------------------------------*\
                             Class Lafaurie
\*--------------------------------------------------------------------*/

class Lafaurie
:
    public surfaceTensionForceModel
{

// Private Member Functions

    // This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //- Surface tension forces on cell faces
    surfaceScalarField Fstffv;
    const label nSmoothingAlpha_;

public:

    //- Runtime type information
    TypeName("Lafaurie");

    // Constructors

    //- construct from components
    Lafaurie
    (
        const word& name,
        const dictionary& surfaceTensionForceProperties,
        const incompressibleThreePhaseMixture& mixture,
	const threePhaseInterfaceProperties& interface
    );


    //- Destructor
    ~Lafaurie()
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
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionForceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
