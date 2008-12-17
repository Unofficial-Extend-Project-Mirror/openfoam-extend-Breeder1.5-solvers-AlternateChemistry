//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Transient/alternateReactingFoam/chemistry.H 4040 2008-10-31T12:22:18.094531Z bgschaid  $ 
{
    Info << "Solving chemistry" << endl;

    chemistry().solve
    (
        runTime.value() - runTime.deltaT().value(),
        runTime.deltaT().value()
    );

    // turbulent time scale
    if (turbulentReaction)
    {
        volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
        volScalarField tc = chemistry().tc();
	Info << "ChemicalTimeScale tc in [ " << min(tc).value() << " ; " << max(tc).value() << "] s" << endl;

        // Chalmers PaSR model
        kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
    }
    else
    {
        kappa = 1.0;
    }
}
