{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

//--------------------------------------------------------------//
/*
    dimensionedScalar maxU
    (
	"maxU",
	dimensionSet(0,1,-1,0,0,0,0),
	gMax
	(
	    (fvc::interpolate(flagField)*mag(phi/mesh.magSf()))()
	)
    );
*/

    dimensionedScalar maxU
    (
	"maxU",
	dimensionSet(0,1,-1,0,0,0,0),
	gMax
	(
	    mag(phi/mesh.magSf())()
	)
    );

//--------------------------------------------------------------//
    surfaceScalarField phir
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
//	mixture.cAlpha()*mag(phi/mesh.magSf())*mixture.nHatfSmooth()
	mixture.cAlpha()*mag(phi/mesh.magSf())*mixture.nHatf()
    );

///////////////////////////////////////////////////////////////////////////////////////
/*
    dimensionedScalar averU
    (
	"averU",
	dimensionSet(0,1,-1,0,0,0,0),
	gSum
	(
	    (fvc::interpolate(flagField)*mag(phi/mesh.magSf()))()
	)/
	gSum(fvc::interpolate(flagField)())
    );
*/


    dimensionedScalar averU
    (
	"averU",
	dimensionSet(0,1,-1,0,0,0,0),
	gAverage
	(
	    mag(phi/mesh.magSf())()
	)
    );


//----------------------------------------------------------------------//

/*
volScalarField averU
(
    IOobject
    (
        "averU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("averU", dimensionSet(0,1,-1,0,0,0,0), gAverage(mag(phi/mesh.magSf())()))
);


    const volScalarField::Boundary& alpha1bf =
        alpha1.boundaryField();
    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchI)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1bf[patchI]))
	{
//	    const polyPatch& cPatch = boundary[patchI];
	    const polyPatch& cPatch = mesh.boundaryMesh()[patchI];
	    forAll(cPatch, faceI)
	    {
		label faceCellI = cPatch.faceCells()[faceI];
		averU[faceCellI] = 0.0;

//		labelList neighCells = mesh.cellCells()[faceCellI];
//		forAll(neighCells, nCellI)
//		{
//		    averU[nCellI] = 0.0;
//		}
	    }
	}
    }
*/

//----------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////

//*--------------------------------------------------------------------*//
/*
    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1[patchi]))
        {

            const fvPatchVectorField& nHatp =
//                refCast<fvPatchScalarField>
                nHat.boundaryField()[patchi];

vectorField& nHatpIn = nHatp.patchInternalField()();


	    forAll(nHatpIn, cellI)
	    {
		nHatpIn[cellI] = vector(0.0, 0.0, 0.0);
	    }
	}
    }
*/
//*--------------------------------------------------------------------*//


    surfaceScalarField phir2
    (
        IOobject
        (
            "phir2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	min
	(
            mixture.cAlpha2()*
	    (
//	        fvc::interpolate(pos0(alpha1 - EvapThresh.value())*neg0(flagField - 0.9)*averU)
	        pos(fvc::interpolate(flagField + liquidField))*averU
	    )
	    +
            mixture.cAlpha3()*
	    (
	        mag(phi/mesh.magSf())*pos(fvc::interpolate(flagField + liquidField))
//	        mag(phi/mesh.magSf())
	    )
	    ,
	    maxU
	)
	*nHatf
//	*mixture.nHatf()
    );


/*
    surfaceScalarField phir2
    (
        IOobject
        (
            "phir2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.cAlpha()*
	(
	    mag(fvc::interpolate(flagField))*averU + mag(phi/mesh.magSf())
	)*mixture.nHatf()
    );
*/


////////////////////////////////////////////////////////////////////////////////////////

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {

        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );

        // Create the complete convection flux for alpha1
        surfaceScalarField alphaPhi1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )

          + fvc::flux
            (
                -fvc::flux(-phir, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
        );
//--------------------------------------------------------------//
/*
        // Create the bounded (upwind) flux for alpha1
        surfaceScalarField alphaPhi1BD
        (
            upwind<scalar>(mesh, phi).flux(alpha1)
        );

        // Calculate the flux correction for alpha1
        alphaPhi1 -= alphaPhi1BD;

        // Calculate the limiter for alpha1
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha1,
                alphaPhi1BD,
                alphaPhi1,
	    	(-mDotEvap/rho2 + divU)(),
	    	(-mDotEvap*(alpha2 + alpha3)/rho1)(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha1,
                alphaPhi1BD,
                alphaPhi1,
	    	(-mDotEvap/rho2 + divU)(),
	    	(-mDotEvap*(alpha2 + alpha3)/rho1)(),
                oneField(),
                zeroField()
            );
        }

        alphaPhi1 = alphaPhi1BD + lambda*alphaPhi1;
        // Solve for alpha1
        solve(fvm::ddt(alpha1) + fvc::div(alphaPhi1));
*/
//--------------------------------------------------------------------//

	MULES::explicitSolve
	(
	    geometricOneField(),
	    alpha1,
	    phi,
	    alphaPhi1,
//	    (-mDotEvap*(1.0/rho2 - 1.0/rho1) + divU)(),
//	    ( - mDotEvap/rho1)(),
	    (-mDotEvap/rho2 + divU)(),
	    (-mDotEvap*(alpha2 + alpha3)/rho1)(),
//	    zeroField(),
//	    ( - mDotEvap/rho1)(),
            oneField(),
            zeroField()
	);

//-------------------------------------------------------------------------//

        // Create the complete flux for alpha2
        surfaceScalarField alphaPhi2
        (
            fvc::flux
            (
                phi,
                alpha2,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir2, alpha1, alpharScheme),
                max(alpha2, 0.0),
                alpharScheme
            )
//------------------------------------------------------//
/*
          + fvc::flux
            (
                -fvc::flux(phir2, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
*/
//-------------------------------------------------------//

        );

//////////////////////////////////////////////////////////////////////////////////////////////

        // Create the bounded (upwind) flux for alpha2
        surfaceScalarField alphaPhi2BD
        (
            upwind<scalar>(mesh, phi).flux(alpha2)
        );

        // Calculate the flux correction for alpha2
        alphaPhi2 -= alphaPhi2BD;


        // Further limit the limiter for alpha2
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha2,
                alphaPhi2BD,
                alphaPhi2,
                (-mDotEvap*(1.0/rho2 - 1.0/rho1) + divU)(),
                (mDotEvap/rho2)(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha2,
                alphaPhi2BD,
                alphaPhi2,
                (-mDotEvap*(1.0/rho2 - 1.0/rho1) + divU)(),
                (mDotEvap/rho2)(),
                oneField(),
                zeroField()
            );
        }

        // Construct the limited fluxes
        alphaPhi2 = alphaPhi2BD + lambda*alphaPhi2;

//////////////////////////////////////////////////////////////////////////////////////////////



        // Create the diffusion coefficients for alpha2<->alpha3
//        volScalarField Dc23(D23*max(alpha3, scalar(0))*pos0(alpha2));
//        volScalarField Dc32(D23*max(alpha2, scalar(0))*pos0(alpha3));
//        volScalarField Dc(D23*pos(alpha2)*pos(alpha3));
//        volScalarField Dc(D23*neg(alpha1 - EvapThresh.value()));
/*
        volScalarField Dcf
	(
	    D23*
	    min
	    (
		neg0(alpha1 - EvapThresh.value()) + flagField,
		1.0
	    )
	);
*/

/*
    tensor I
    (
	1, 0, 0,
	0, 1, 0,
	0, 0, 1
    );

volTensorField Dc = 
    	    D23*
	    min
	    (
		neg0(alpha1 - EvapThresh.value()) + flagField,
		1.0
	    )*I;
*/

/*
Dc23 = Dc23*	    
    min
    (
	neg0(alpha1 - EvapThresh.value()) + flagField,
	1.0
    );

Dc32 = Dc32*	    
    min
    (
	neg0(alpha1 - EvapThresh.value()) + flagField,
	1.0
    );
*/

    surfaceScalarField Dcf
    (
        IOobject
        (
            "Dcf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	D23*neg0(fvc::interpolate(flagField + liquidField))
    );


        // Add the diffusive flux for alpha3->alpha2
//        alphaPhi2 += 0.5*D23*mesh.magSf()*fvc::snGrad(alpha2)*pos(fvc::interpolate(flagField + liquidField))*neg(fvc::interpolate(flagField + liquidField - 1.0)); //
//	alphaPhi2 -= fvc::interpolate(Dc32)*mesh.magSf()*fvc::snGrad(alpha1);

        // Solve for alpha2
        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2)
          + fvc::div(alphaPhi2)
//          - fvm::laplacian(Dc23 + Dc32, alpha2)
          - fvm::laplacian(Dcf, alpha2)
	  - fvm::Sp(-mDotEvap*(1.0/rho2 - 1.0/rho1), alpha2)
	  - fvm::Sp(divU, alpha2)
	  - (mDotEvap/rho2)
        );
        alpha2Eqn.solve();


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// solve alpha3 //////////////////////////////////////////////////////////////////////////

/*
        // Create the complete flux for alpha3
        surfaceScalarField alphaPhi3
        (
            fvc::flux
            (
                phi,
                alpha3,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir2, alpha1, alpharScheme),
                max(alpha3, 0.0),
                alpharScheme
            )
        );

        // Create the bounded (upwind) flux for alpha3
        surfaceScalarField alphaPhi3BD
        (
            upwind<scalar>(mesh, phi).flux(alpha3)
        );

        // Calculate the flux correction for alpha3
        alphaPhi3 -= alphaPhi3BD;


        // Further limit the limiter for alpha3
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha3,
                alphaPhi3BD,
                alphaPhi3,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha3,
                alphaPhi3BD,
                alphaPhi3,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        // Construct the limited fluxes
        alphaPhi3 = alphaPhi3BD + lambda*alphaPhi3;

        // Solve for alpha3
        fvScalarMatrix alpha3Eqn
        (
            fvm::ddt(alpha3)
          + fvc::div(alphaPhi3)
          - fvm::laplacian(Dc, alpha3)
//	  - fvm::Sp(-mDotEvap*(1.0/rho2 - 1.0/rho1), alpha2)
//	  - fvm::Sp(divU, alpha2)
//	  - (mDotEvap/rho2)
        );
        alpha3Eqn.solve();

alpha1 = alpha1/(alpha1 + alpha2 + alpha3);
alpha2 = alpha2/(alpha1 + alpha2 + alpha3);
//alpha3 = alpha1/(alpha1 + alpha2 + alpha3);

*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////// end of solve alpha3 ////////////////////////////////////////////////////////////////////////

        // Construct the complete mass flux
        rhoPhi =
              alphaPhi1*(rho1 - rho3)
            + (alphaPhi2 + alpha2Eqn.flux())*(rho2 - rho3)
            + phi*rho3;

        alpha3 = max(1.0 - max(alpha1, 0.0) - max(alpha2, 0.0), 0.0);
alpha1Phi = alphaPhi1;
alpha2Phi = alphaPhi2 + alpha2Eqn.flux();
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "Vapor phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << endl;

    Info<< "Air phase volume fraction = "
        << alpha3.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha3.name() << ") = " << min(alpha3).value()
        << "  Max(" << alpha3.name() << ") = " << max(alpha3).value()
        << endl;
}
