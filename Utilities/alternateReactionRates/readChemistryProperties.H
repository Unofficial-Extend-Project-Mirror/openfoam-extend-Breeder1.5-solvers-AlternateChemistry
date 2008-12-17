//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Utilities/alternateReactionRates/readChemistryProperties.H 3877 2008-10-05T19:22:47.837696Z bgschaid  $ 
Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch turbulentReaction(chemistryProperties.lookup("turbulentReaction"));

dimensionedScalar Cmix("Cmix", dimless, 1.0);

if (turbulentReaction)
{
    chemistryProperties.lookup("Cmix") >> Cmix;
}

scalar steadyChemistryDeltaT=readScalar(chemistryProperties.lookup("steadyChemistryDeltaT"));
