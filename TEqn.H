volScalarField k = alpha1*mixture.k1() + alpha2*mixture.k2() + (1.0 - alpha1 - alpha2)*mixture.k3();
volScalarField rhoCp = alpha1*rho1*mixture.cp1() + alpha2*rho2*mixture.cp2() + (1.0 - alpha1 - alpha2)*rho3*mixture.cp3();
rhoCp.oldTime() = alpha1.oldTime()*rho1*mixture.cp1() + alpha2.oldTime()*rho2*mixture.cp2() + (1.0 - alpha1.oldTime() - alpha2.oldTime())*rho3*mixture.cp3();

//surfaceScalarField alpha2Phi = (rhoPhi - phi*mixture.rho3() - alpha1Phi*(mixture.rho1() - mixture.rho3()))/(mixture.rho2() - mixture.rho3());

surfaceScalarField rhoCpPhi = alpha1Phi*(mixture.rho1()*mixture.cp1() - mixture.rho3()*mixture.cp3()) + alpha2Phi*(mixture.rho2()*mixture.cp2() - mixture.rho3()*mixture.cp3()) + phi*mixture.rho3()*mixture.cp3();





//------------------------------------------------------------//
/*
volScalarField k = 1.0/
	(
	    alpha1/mixture.k1() + alpha2/mixture.k2() + (1.0 - alpha1 - alpha2)/mixture.k3()
	);
volScalarField rhoCp = alpha1*rho1*mixture.cp1() + alpha2*rho2*mixture.cp2() + (1.0 - alpha1 - alpha2)*rho3*mixture.cp3();
*/
//----------------------------------------------------------------//


// energy equation
fvScalarMatrix TEqn
(
//	fvm::ddt(T) + fvm::div(phi, T) == fvm::laplacian(k/rhoCp, T) + hDotEvap/rhoCp + fvm::Sp(fvc::div(phi), T)

	fvm::ddt(rhoCp, T) + fvm::div(rhoCpPhi, T) == fvm::laplacian(k, T) + hDotEvap + fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)

//	fvm::ddt(rhoCp, T) + fvm::div(rhoCpPhi, T) == fvm::laplacian(k, T) + hDotEvap
);


TEqn.relax();

TEqn.solve();
//T.correctBoundaryConditions();
							Info << "T = " << min(T) <<endl;

							Info << "T = " << max(T*flagField) <<endl;

