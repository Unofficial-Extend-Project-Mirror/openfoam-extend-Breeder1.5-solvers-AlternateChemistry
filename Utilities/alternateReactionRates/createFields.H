//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Utilities/alternateReactionRates/createFields.H 4044 2008-11-01T22:28:04.889474Z bgschaid  $ 
    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<hCombustionThermo> thermo
    (
        hCombustionThermo::New(mesh)
    );

    combustionMixture& composition = thermo->composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermo->rho()
    );

    PtrList<volScalarField> RR(Y.size());

    forAll(RR,yI) {
        RR.set(
            yI,
            new volScalarField(
                IOobject
                (
                    "ReactionRate"+Y[yI].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("source",dimensionSet(1,-3,-1,0,0,0,0),0)
            )
        );
    }

    volScalarField RRsum
    (
        IOobject
        (
            "ReactionRateSum",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0)
    );

Info << "Constructing chemical mechanism" << endl;
autoPtr<alternateChemistryModel> chemistry=alternateChemistryModel::New
                                             (
                                                 thermo(),
                                                 rho,
                                                 steady
                                             );
