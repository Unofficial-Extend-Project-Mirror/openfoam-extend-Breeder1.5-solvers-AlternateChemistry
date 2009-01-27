//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/chemistry.H 4040 2008-10-31T12:22:18.094531Z bgschaid  $ 
{
    Info << "Solving chemistry" << endl;

    scalar chemTime=chemistry().solve
    (
        runTime.value() - steadyChemistryDeltaT*runTime.deltaT().value(),
        steadyChemistryDeltaT*runTime.deltaT().value()
    );

    Info << " Characteristic time of chemistry: " << chemTime << endl;
        
    // turbulent time scale
    if (turbulentReaction)
    {
        volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
        volScalarField tc = chemistry().tc();
	volScalarField tf = chemistry().tf();	
	
        // Chalmers PaSR model
        //kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
	kappa = (tf + tc)/(tf + tc + tk);
    }
    else
    {
        kappa = 1.0;
    }
}
