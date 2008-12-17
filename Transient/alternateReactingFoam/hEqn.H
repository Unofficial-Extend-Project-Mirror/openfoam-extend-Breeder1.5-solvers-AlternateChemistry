//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Transient/alternateReactingFoam/hEqn.H 3626 2008-08-28T14:34:55.168701Z bgschaid  $ 
{
    solve
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        DpDt
    );

    thermo->correct();

    chemistry().calcDQ(dQ);
}
