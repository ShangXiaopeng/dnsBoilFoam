
#include "interIdentification.H"

scalar cutoff = 0.01;

rhogp =	(rho3*alpha3 + rho2*alpha2 + 0.5*SMALL*(rho2+rho3))/
	(alpha2 + alpha3 +SMALL)*neg(alpha1 + cutoff - 1.0)*pos0(alpha2);

Xv =	(alpha2*rho2)/
	(
	    alpha2*rho2 + alpha3*rho3 + 0.5*SMALL*(rho2+rho3)
	)*neg(alpha1 - 1.0 + cutoff)*pos0(alpha2);

// volScalarField pSat = exp(11.495 - mixture.Tc()/617.65*5141.36/T)*101325*mixture.Pc()/21.05e5;


volScalarField tau(1.0-T/mixture.Tc());

volScalarField pSat = 
exp
(
    (mixture.a()*tau + mixture.b()*pow(tau, 1.5) + mixture.c()*pow(tau, 3.0) + mixture.d()*pow(tau, 6.0))
//    (mixture.a()*tau + mixture.b()*pow(tau, 1.5) + mixture.c()*pow(tau, 2.5) + mixture.d()*pow(tau, 5.0))
   /(1.0 - tau)
)*mixture.Pc();


volScalarField& pStat(p_rgh);

volScalarField XvSat =  
(
    (pSat/pStat)*mixture.MW2()/
    (
	(pSat/pStat)*mixture.MW2() + (1.0 - (pSat/pStat))*mixture.MW3()
    )
);
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
forAll(alpha1, cellI)
{
    if (alpha1[cellI] >= EvapThresh.value())
    {
	Xv[cellI] = XvSat[cellI];
    }
}


gradXv = fvc::grad(Xv)*flagField;
//gradXv = (fvc::grad(XvSat*alpha1 + (1.0 - alpha1)*Xv))*flagField;
*/

#include "calcGradXv.H"
Info << "gradXv*nHat = " << min(gradXv) << endl;
Info << "gradXv*nHat = " << max(gradXv) << endl;


Info << "rhogp = " << max(rhogp*flagField) << endl;
Info << "Xv = " << max(Xv) << endl;
Info << "p_rgh = " << min(p_rgh) << endl;


////////////////////////////////////////////////////////////////////////////////////////////////////////////

gradAlpha1 = mag(fvc::grad(alpha1));
gradAlpha1f = fvc::interpolate(gradAlpha1);
gradAlpha1 =fvc::surfaceSum(gradAlpha1f);

/*
dimensionedScalar N("N", dimensionSet(0,0,0,0,0,0,0), scalar(0.));
N = fvc::domainIntegrate(mag(gradAlpha1)).value()/
max
(
    fvc::domainIntegrate(flagField*mag(gradAlpha1)).value(), 
    scalar(1e-12)
);
Info << "N = " << N.value() << endl;
*/

mDotEvap = 
max
(
    min(
            D23*rhogp/(1.0 - Xv + SMALL)
//             D23*rhogp/(1.0 - XvSat + SMALL)
          * gradXv*flagField
//            *mag(fvc::grad(alpha1)),
          * gradAlpha1,
            alpha1*rho1/runTime.deltaT()
       ),
    mDotZero
);
mDotEvap = mDotZero;

hDotEvap = -mixture.hEvap()*mDotEvap;

//Info << "mDotEvap = " << min(mDotEvap) << endl;
Info << "mDotEvap = " << max(mDotEvap) << endl;
Info << "hDotEvap = " << min(hDotEvap) << endl;



