
MeshGraph InterfaceMeshGraph(mesh, alpha1);

std::vector<MeshGraph::CellFacePair> EvapIntCellFacePairs;

// Find internal interface cell pairs using graph traversal
InterfaceMeshGraph.Reset();

InterfaceMeshGraph.GetInterfaceCellFacePairs
(
    EvapIntCellFacePairs, 
    EvapThresh.value()
);

// Reset interface field, then interpolate
flagField = 0.0;
liquidField = -1.0;
// Loop through evap cells:
for
(
    std::vector<MeshGraph::CellFacePair>::iterator it =
        EvapIntCellFacePairs.begin();
    it != EvapIntCellFacePairs.end();
    it++
)
{

    flagField[(*it).c1] = 1.0;  
//    flagField[(*it).c2] = 1.0;

}


if(Pstream::parRun())
{
//    #include "adjacentCellsPar.H"


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

	    const labelList& internalcells = patch.faceCells();
	    forAll(internalcells, ind)
	    {
		label curcell = internalcells[ind];
		mybuffer_alpha1[ind] = alpha1[curcell];
		mybuffer_flagField[ind] = flagField[curcell];
//		patchboundary[curcell] = procpatch.myProcNo(); //for debugging
	    }

	    //Send buffer to neighbor
	    OPstream tNP(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 0);
	    tNP << mybuffer_alpha1 << endl;
	    OPstream tNP2(UPstream::commsTypes::blocking, procpatch.neighbProcNo(), patch.size()*8, 1);
	    tNP2 << mybuffer_flagField << endl;
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
	    {
		const labelList& internalcells = patch.faceCells();
		forAll(internalcells, cellJ)
		{
		    label curcell = internalcells[cellJ];
		    alpha1[curcell] = yourbuffer_alpha1[cellJ];
		    flagField[curcell] = yourbuffer_flagField[cellJ];
		}
	    }
	}
    }
}


}


{

    forAll(flagField, cellI)
    {
        if(flagField[cellI] == 0 && alpha1[cellI] > EvapThresh.value())
        {
            labelHashSet setNBCells(1);
            labelList lCP = mesh.cellEdges()[cellI];

            forAll(lCP, edgeI)
            {
	        labelList lPC = mesh.edgeCells()[lCP[edgeI]];
	        setNBCells.insert(lPC);
            }

            labelList NBCells = setNBCells.toc();

            for(int j = 0; j < NBCells.size(); j++)
            {
	        if(alpha1[NBCells[j]] < EvapThresh.value())
	        {
	            flagField[cellI] = 1.0;
	            break;
	        }
	
            }
        }
    }


    forAll(liquidField, cellI)
    {
        if(flagField[cellI] == 0 && alpha1[cellI] > EvapThresh.value())
        {
	    liquidField[cellI] = 1.0;
        }
    }

}





