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
        labelList lCP = mesh.cellEdges()[cellI];

        forAll(lCP, edgeI)
        {
            labelList lPC = mesh.edgeCells()[lCP[edgeI]];
            setNBCells.insert(lPC);
        }

        labelList NBCells = setNBCells.toc();

//	label marker = -2000;

        for(int j = 0; j < NBCells.size(); j++)
        {
	    if(flagField[NBCells[j]] == 0.0 && alpha1[NBCells[j]] < EvapThresh.value())
	    {
		vector d = mesh.C()[NBCells[j]] - mesh.C()[cellI];
	        vector nd = d/mag(d);

		if((nHat[cellI]&nd) > 0.0 && pos(XvSat[cellI] - Xv[NBCells[j]]))
		{
		    gradXv[NBCells[j]] = (gradXv[NBCells[j]]*countField[NBCells[j]] + (XvSat[cellI] - Xv[NBCells[j]])/mag(d)*(nd&nHat[cellI]))/(countField[NBCells[j]] + 1.0);
		    countField[NBCells[j]] += 1.0;

//		    marker = NBCells[j];
		}
	    }
        }

//	if(marker < 0)
//	{
//	    Info << "The angle between the interface normal and the closest cell larger than 90 degree!!!";
//	}
    }
}


forAll(mesh.cells(), cellI)
{

    if(flagField[cellI] == 1.0)
    {
        labelHashSet setNBCells(1);
        labelList lCP = mesh.cellEdges()[cellI];

        forAll(lCP, edgeI)
        {
            labelList lPC = mesh.edgeCells()[lCP[edgeI]];
            setNBCells.insert(lPC);
        }

        labelList NBCells = setNBCells.toc();

	label marker = -2000;

        for(int j = 0; j < NBCells.size(); j++)
        {
	    if(flagField[NBCells[j]] == 0.0 && alpha1[NBCells[j]] < EvapThresh.value())
	    {
		vector d = mesh.C()[NBCells[j]] - mesh.C()[cellI];
	        vector nd = d/mag(d);

		if((nHat[cellI]&nd) > 0.0 && pos(XvSat[cellI] - Xv[NBCells[j]]))
		{
		    gradXv[cellI] = (gradXv[cellI]*countField[cellI] + gradXv[NBCells[j]])/(countField[cellI] + 1.0);
		    Xv[cellI] = (Xv[cellI]*countField[cellI] + Xv[NBCells[j]])/(countField[cellI] + 1.0);
		    rhogp[cellI] = (rhogp[cellI]*countField[cellI] + rhogp[NBCells[j]])/(countField[cellI] + 1.0);
		    countField[cellI] +=1.0;

		    marker = NBCells[j];
		}
	    }
        }

//        Xv[cellI] = (XvSat[cellI] + Xv[cellI])/2.0;
        Xv[cellI] = XvSat[cellI];
	rhogp[cellI] = (rho2.value() + rhogp[cellI])/2.0;



	if(marker < 0)
	{
	    Info << "The angle between the interface normal and the closest cell larger than 90 degree!!!";
	}
    }
}


















