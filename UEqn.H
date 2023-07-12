    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      + turbulence->divDevRhoReff(rho, U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
//                    mixture.Fstff()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }

fSigma = fvc::reconstruct
            (
                mixture.surfaceTensionForce()* mesh.magSf()
            );
//							Info << "sigma = " << min(mixture.surfaceTensionForce()) <<endl;

//							Info << "sigma = " << max(mixture.surfaceTensionForce()) <<endl;
