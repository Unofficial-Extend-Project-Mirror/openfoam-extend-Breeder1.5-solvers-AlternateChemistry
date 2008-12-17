//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/createFields.H 4044 2008-11-01T22:28:04.889474Z bgschaid  $ 
Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<hCombustionThermo> thermo
(
    hCombustionThermo::New(mesh)
);

combustionMixture& composition = thermo->composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo->lookup("inertSpecie"));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo->rho()
);


volScalarField& p = thermo->p();
//const volScalarField& psi = thermo->psi();
const volScalarField& T = thermo->T();
volScalarField& h = thermo->h();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

dimensionedScalar pMin
(
    mesh.solutionDict().subDict("SIMPLE").lookup("pMin")
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::RASModel> turbulence
(
    compressible::RASModel::New
    (
        rho,
        U,
        phi,
        thermo()
    )
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

// Info<< "Creating field DpDt\n" << endl;
// volScalarField DpDt = 
//     fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

Info << "Constructing chemical mechanism" << endl;
autoPtr<alternateChemistryModel> chemistry=alternateChemistryModel::New
                                             (
                                                 thermo(),
                                                 rho,
                                                 true
                                             );

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

for(label i=0; i<Y.size(); i++)
{
    fields.add(Y[i]);
}
fields.add(h);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
);
