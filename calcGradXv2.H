if(Pstream::parRun())
{

//Parallel
// Loop over processor patches
Info << " Sending stuff " << endl ;
forAll (mesh.boundaryMesh(), patchInd)
{
    const polyPatch& patch = mesh.boundaryMesh()[patchInd];
    if (typeid(patch) == typeid(processorPolyPatch))
    {
	const processorPolyPatch& procpatch = dynamic_cast<const processorPolyPatch&>(patch);
	//only for master chunk
	if (procpatch.myProcNo() > procpatch.neighbProcNo())
	{
	    //Make buffer
	    Field<scalar> mybuffer_alpha1(patch.size());
	    Field<scalar> mybuffer_flagField(patch.size());
	    Field<scalar> mybuffer_countField(patch.size());

	    const labelList& internalcells = patch.faceCells();
	    forAll(internalcells, ind)
	    {
		label curcell = internalcells[ind];
		mybuffer_alpha1[ind] = alpha1[curcell];
		mybuffer_flagField[ind] = flagField[curcell];
		mybuffer_countField[ind] = countField[curcell];
//		patchboundary[curcell] = procpatch.myProcNo(); //for debugging
	    }

	    //Send buffer to neighbor
	    OPstream tNP(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 0);
	    tNP << mybuffer_alpha1 << endl;
	    OPstream tNP2(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 1);
	    tNP2 << mybuffer_flagField << endl;
	    OPstream tNP3(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 2);
	    tNP3 << mybuffer_countField << endl;
	}
    }
}


Info << "Messages sent" << endl;
//Receive messages
forAll (mesh.boundaryMesh(), patchInd)
{
    const polyPatch& patch = mesh.boundaryMesh()[patchInd];
    if (typeid(patch) == typeid(processorPolyPatch))
    {
	const processorPolyPatch& procpatch = dynamic_cast<const processorPolyPatch&>(patch);
	if (procpatch.myProcNo() < procpatch.neighbProcNo())
	{ 
	    //only for slave chunk
	    //Make buffer
	    Field<scalar> yourbuffer_alpha1(patch.size());
	    IPstream fNP(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 0);
	    fNP >> yourbuffer_alpha1;

	    Field<scalar> yourbuffer_flagField(patch.size());
	    IPstream fNP2(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 1);
	    fNP2 >> yourbuffer_flagField;

	    Field<scalar> yourbuffer_countField(patch.size());
	    IPstream fNP3(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 2);
	    fNP3 >> yourbuffer_countField;

	    {
		const labelList& internalcells = patch.faceCells();
		forAll(internalcells, cellJ)
		{
		    label curcell = internalcells[cellJ];
		    alpha1[curcell] = yourbuffer_alpha1[cellJ];
		    flagField[curcell] = yourbuffer_flagField[cellJ];
		    countField[curcell] = yourbuffer_countField[cellJ];
		}
	    }
	}
    }
}

}



countField = 0.0;
gradXv = dimensionedScalar
(
    "gradXv",
    dimensionSet(0,-1,0,0,0,0,0),
    scalar(0.0)
);

forAll(mesh.cells(), cellI)
{

    if(flagField[cellI] == 1.0)
    {
        labelHashSet setNBCells(1);
//        labelList lCP = mesh.cellEdges(cellI);
        labelList lCP = mesh.cellEdges()[cellI];

        forAll(lCP, edgeI)
        {
//          labelList lPC = mesh.edgeCells(lCP[edgeI]);
            labelList lPC = mesh.edgeCells()[lCP[edgeI]];

            setNBCells.insert(lPC);
        }

        labelList NBCells = setNBCells.toc();

	scalar COS(-1.0);
	label marker = -2000;

        for(int j = 0; j < NBCells.size(); j++)
        {
	    if(flagField[NBCells[j]] == 0.0 && alpha1[NBCells[j]] < EvapThresh.value())
	    {
		vector d = mesh.C()[NBCells[j]] - mesh.C()[cellI];
	        vector nd = d/mag(d);

		if((nHat[cellI]&nd) > COS)
		{
		    gradXv[cellI] = (XvSat[cellI] - Xv[NBCells[j]])/mag(d)*(nd&nHat[cellI]);
		    COS = nd&nHat[cellI];
//		    countField[NBCells[j]] = 1.0;
		    Xv[cellI] = (XvSat[cellI] + Xv[NBCells[j]])/2.0;
		    rhogp[cellI] = (rho2.value() + rhogp[NBCells[j]])/2.0;
		    marker = NBCells[j];
		}
	    }
        }
/*
if(gradXv[cellI] < 0)
{
	Info << "cellI = " << cellI << endl;
	Info << "XvSat = " << XvSat[cellI] << endl;
	Info << "Xv = " << Xv[marker] << endl;
	Info << "alpha1[cellI] = " << alpha1[cellI] << endl;
	Info << "alpha1[marker] = " << alpha1[marker] << endl;
	Info << "alpha2[cellI] = " << alpha2[cellI] << endl;
	Info << "alpha2[marker] = " << alpha2[marker] << endl;
	Info << "T[cellI] = " << T[cellI] << endl;
		    countField[marker] = 1.0;


        labelHashSet setNBCells2(1);
//        labelList lCP = mesh.cellEdges(cellI);
        labelList lCP2 = mesh.cellEdges()[marker];

        forAll(lCP2, edgeI)
        {
//          labelList lPC = mesh.edgeCells(lCP[edgeI]);
            labelList lPC2 = mesh.edgeCells()[lCP2[edgeI]];

            setNBCells2.insert(lPC2);
        }

        labelList NBCells2 = setNBCells2.toc();
        for(int j = 0; j < NBCells2.size(); j++)
	{
	    if(flagField[NBCells2[j]] == 1.0)
	    Info << "T of NBCells = " << T[NBCells2[j]] << endl;
	}

}
*/

	countField[marker] = 1.0;
	if(COS < 0)
	{
	    Info << "The angle between the interface normal and the closest cell larger than 90 degree!!!";
	}
    }
}




















