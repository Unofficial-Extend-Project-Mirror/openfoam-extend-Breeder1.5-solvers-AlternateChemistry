//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/hEqn.H 4036 2008-10-30T15:07:24.059047Z bgschaid  $ 
{
    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p, "div(U,p)"))
      - p*fvc::div(phi/fvc::interpolate(rho))
    );

    hEqn.relax();

    if(fixedCells.size()>0) {
        Info << "Fixing enthalpy in " << fixedCells.size() << " cells" << endl;

        scalarField hFix=dynamic_cast<basicThermo&>(thermo()).h(fixedTemp,fixedCells);
        hEqn.setValues(fixedCells,hFix);
    }

    eqnResidual = hEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    thermo->correct();

    chemistry().calcDQ(dQ);
}
