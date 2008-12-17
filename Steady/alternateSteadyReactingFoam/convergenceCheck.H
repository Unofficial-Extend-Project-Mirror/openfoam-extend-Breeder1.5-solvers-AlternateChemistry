//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/convergenceCheck.H 3673 2008-09-05T13:46:35.475871Z bgschaid  $ 
// check convergence

if (maxResidual < convergenceCriterion)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;
    runTime.writeAndEnd();
    Info<< "latestTime = " << runTime.timeName() << endl;
}

