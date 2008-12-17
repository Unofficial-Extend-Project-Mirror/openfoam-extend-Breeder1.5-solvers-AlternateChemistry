//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Transient/alternateReactingFoam/readChemistryProperties.H 3602 2008-08-21T16:17:37.921551Z bgschaid  $ 
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
