//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/UEqn.H 3673 2008-09-05T13:46:35.475871Z bgschaid  $ 
    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::Sp(fvc::div(phi), U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
    );

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
