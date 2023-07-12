# dnsEvapFoam
a solver for direct numerical simulation of vaporization problems
Including liquid, vapor and inert gas
The vaporization rate determined directly by physical variables in the solution
A sharp-interface algebraic algorithm used to resolve the liquid-gas interface
Please feel free to contact Dr. Shang (xshang002@e.ntu.edu.sg) for any questions
You are welcome to cite the paper "X. Shang, X. Zhang, T.B. Nguyen, T. Tran, Direct numerical simulation of evaporating droplets based on a sharp-interface algebraic VOF approach, International Journal of Heat and Mass Transfer 184: 122282"

Make/files
MeshGraph/MeshGraph.C
incompressibleThreePhaseMixture/incompressibleThreePhaseMixture.C
threePhaseInterfaceProperties/threePhaseInterfaceProperties.C

surfaceTensionForceModels/surfaceTensionForceModel/surfaceTensionForceModel.C
surfaceTensionForceModels/surfaceTensionForceModel/newSurfaceTensionForceModel.C
surfaceTensionForceModels/Brackbill/Brackbill.C
surfaceTensionForceModels/SST/SST.C
surfaceTensionForceModels/Lafaurie/Lafaurie.C
surfaceTensionForceModels/SmoothedSF/SmoothedSF.C
surfaceTensionForceModels/temperatureDependentBrackbill/temperatureDependentBrackbill.C

immiscibleIncompressibleThreePhaseMixture/immiscibleIncompressibleThreePhaseMixture.C
dnsEvapFoam.C

EXE = $(FOAM_USER_APPBIN)/dnsEvapFoam

Make/options
EXE_INC = \
    -I. \
    -I.. \
    -I../../VoF \
    -I$(LIB_SRC)/transportModels/twoPhaseMixture/lnInclude \
    -IimmiscibleIncompressibleThreePhaseMixture \
    -IincompressibleThreePhaseMixture \
    -IthreePhaseInterfaceProperties \
    -IMeshGraph \
    -IRiddersRoot \
    -I$(LIB_SRC)/transportModels/interfaceProperties/lnInclude \
    -I$(LIB_SRC)/transportModels/twoPhaseProperties/alphaContactAngle/alphaContactAngle \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/dynamicFvMesh/lnInclude \
    -I$(LIB_SRC)/transportModels \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/sampling/lnInclude \
-IsurfaceTensionForceModels/surfaceTensionForceModel \
-IsurfaceTensionForceModels/Brackbill \
-IsurfaceTensionForceModels/SST \
-IsurfaceTensionForceModels/Lafaurie \
-IsurfaceTensionForceModels/SmoothedSF \
-IsurfaceTensionForceModels/temperatureDependentBrackbill

EXE_LIBS = \
    -ltwoPhaseMixture \
    -ltwoPhaseProperties \
    -lincompressibleTransportModels \
    -lturbulenceModels \
    -lincompressibleTurbulenceModels \
    -lfiniteVolume \
    -ldynamicFvMesh \
    -lmeshTools \
    -lfvOptions \
    -lsampling
