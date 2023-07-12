//labeledInterface = 0.0;
//surfaceScalarField Flag(fvc::interpolate(flagField + liquidField));

forAll(mesh.cells(), cellI)
{

    if(flagField[cellI] == 1.0)
    {
        const cell& curCell = mesh.cells()[cellI];
	const labelList& faces = mesh.cells()[cellI];
	const labelList& neiCells = mesh.cellCells()[cellI];

	const labelUList& owners = mesh.faceOwner();
	const labelUList& neighbs = mesh.faceNeighbour();

	surfaceScalarField alpha1f = fvc::interpolate(alpha1);
	const faceList& allFacesList = mesh.faces();

	forAll(neiCells, nCellI)
	{
	    const label& neiCell = neiCells[nCellI];

//	    if(flagField[neiCell] == 0 && alpha1[neiCell] > EvapThresh.value()) // One neighbour cell is liquid
	    if(liquidField[neiCell] == 1.0) // One neighbour cell is liquid
	    {
		const labelList& neiFaces = mesh.cells()[neiCell];

		forAll (neiFaces, neiFaceI)
		{
		    const label& neiFace = neiFaces[neiFaceI];
		    const label& sharedFaceIndex = findIndex(faces, neiFace);

		    if (sharedFaceIndex != -1 ) // the shared face between the interface cell and neighbour cell
		    {
//			surfaceScalarField Tf = fvc::interpolate(T);
//	        	const faceList& allFacesList = mesh.faces();

	        	const label& oppFaceLabel = curCell.opposingFaceLabel(faces[sharedFaceIndex], allFacesList);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------//
/*

			if (mesh.isInternalFace(oppFaceLabel) && (liquidField[owners[oppFaceLabel]] + liquidField[neighb[oppFaceLabel]] == 1.0))
			{
			    continue;
			}

			else
			{

	        	    Tgf[faces[sharedFaceIndex]] = 2.0*T[cellI] - Tf[oppFaceLabel];

			    labeledInterface[cellI] += 1.0;
			}
*/
//---------------------------------------------------------------//

			if (mesh.isInternalFace(oppFaceLabel))
			{
			    const label& own = owners[oppFaceLabel];
			    const label& nei = neighbs[oppFaceLabel];

			    if (liquidField[own] * liquidField[nei] == -1.0)
			    {
Info << "Surrounded by two liquid cells! Ignored! " << endl;
				continue;
			    }
			    else
			    {
				alpha1gf[faces[sharedFaceIndex]] = 2.0*alpha1[cellI] - alpha1f[oppFaceLabel];
//				labeledInterface[cellI] += 1.0;
			    }
			}
			else
			{
Info << "Surrounded by one liquid cell and one boundary cell! " << endl;
			    alpha1gf[faces[sharedFaceIndex]] = 2.0*alpha1[cellI] - alpha1f[oppFaceLabel];
//			    labeledInterface[cellI] += 1.0;
			}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		    }
		}
	    }
	}
    }
}

//Info << "stop 1" << endl;

gradAlpha1g = fvc::grad(alpha1gf);



/*
fTan = 
    sigmaT*fvc::interpolate(mag(fvc::grad(alpha1)))*
    (
	(fvc::interpolate(gradTg) & (mesh.Sf()/mesh.magSf()))
      + fvc::interpolate(gradTg & nHat)*nHatf/mesh.magSf()
    ); // nHat reverse direction

fTan2 = dimensionedScalar("dummy", dimensionSet(1, -2, -2, 0, 0, 0, 0), scalar(0.0));

forAll (interfaceCellFaces, index)
{
    label faceI = interfaceCellFaces[index];
    fTan2[faceI] = fTan[faceI];
}

fSigma = fvc::reconstruct
            (
                fTan2*mesh.magSf()
            )*flagField;




fTan = 
    sigmaT*fvc::interpolate(mag(fvc::grad(alpha1)))*
    (
	(fvc::interpolate(gradTg) & (mesh.Sf()/mesh.magSf()))
      + fvc::interpolate(gradTg & nHat)*nHatf/mesh.magSf()
//      + fvc::interpolate(gradTg & nHat)*(-mixture.nHatf())/mesh.magSf()
    )*pos(fvc::interpolate(flagField)); // nHat reverse direction



fSigma = fvc::reconstruct
            (
                fTan*mesh.magSf()
            );


*/
/*

fTan = 

    sigmaT*fvc::interpolate(mag(fvc::grad(alpha1)))*fvc::interpolate(gradTg - (gradTg & nHat)*nHat) & (mesh.Sf()/mesh.magSf());

fTan2 = dimensionedScalar("dummy", dimensionSet(1, -2, -2, 0, 0, 0, 0), scalar(0.0));

forAll (interfaceCellFaces, index)
{
    label faceI = interfaceCellFaces[index];
    fTan2[faceI] = fTan[faceI];
}

fSigmaR = fvc::reconstruct
            (
                fTan*mesh.magSf()
            );
*/

//fSigma = sigmaT*mag(fvc::grad(alpha1))*(gradTg - (gradTg & nHat)*nHat);




