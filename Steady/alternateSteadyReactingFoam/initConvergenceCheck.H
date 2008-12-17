//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/initConvergenceCheck.H 3673 2008-09-05T13:46:35.475871Z bgschaid  $ 
// initialize values for convergence checks

    scalar eqnResidual = 1, maxResidual = 0;
    scalar convergenceCriterion = 0;

    simple.readIfPresent("convergence", convergenceCriterion);

