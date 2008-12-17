//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Transient/alternateReactingFoam/UEqn.H 3602 2008-08-21T16:17:37.921551Z bgschaid  $ 
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
