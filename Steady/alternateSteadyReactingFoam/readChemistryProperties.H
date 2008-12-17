//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/readChemistryProperties.H 4040 2008-10-31T12:22:18.094531Z bgschaid  $ 
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

labelList fixedCells;
scalarField fixedTemp;

if(chemistryProperties.found("fixedTemperatur")) {
    List<dictionary> ignDict(chemistryProperties.lookup("fixedTemperatur"));
    fixedCells.setSize(ignDict.size());
    fixedTemp.setSize(ignDict.size());
    forAll(ignDict,i) {
        const dictionary &iDi=ignDict[i];
        const vector pos(iDi.lookup("position"));
        fixedTemp[i]=readScalar(iDi.lookup("temperature"));
        scalar minDist=GREAT;
        const volVectorField& centres = mesh.C();

        forAll(centres,celli) {
            scalar dist=mag(pos-centres[celli]);
            if(dist<minDist) {
                minDist=dist;
                fixedCells[i]=celli;
            }
        }
    }
}

//YEqn stabilisation parameters - to prevent steep changes in Y and T if one of the 
//checks finds a problem (positive result of the check) kappa is reduced.

// 5 - number of performed stabilisation loops
label recalcMax   = readLabel(chemistryProperties.lookup("recalcMax"));	

// 1.001  - maximum YTotal=sum(Yi) allowed (see 1st check)
scalar yOverThres = readScalar(chemistryProperties.lookup("YOvershootThreshold"));

// 1e-3 - if deltaY is smaller than this the 2nd check is not performed
scalar minAbsDeltaY = readScalar(chemistryProperties.lookup("minAbsDeltaY"));

// 0.1 - if the absolute change of Y is larger than this the 2nd check turns positive
scalar absDeltaYThres = readScalar(chemistryProperties.lookup("absDeltaYThres"));

// 0.5 - if the relative change of Y is larger than this the 2nd check turns positive
scalar relDeltaYThres = readScalar(chemistryProperties.lookup("relDeltaYThres"));

// 1e-6 - a threshold needed for the calculation of relDeltaY like in case Yold=0
//relDeltaY = absDeltaY/max(relDeltaYCalcThres,Yold);
//first it was set to SMALL but this resulted in huge relDeltaY!
scalar relDeltaYCalcThres = readScalar(chemistryProperties.lookup("relDeltaYCalcThres"));

// 0.25 - kappa is scaled by that value if 1st or 2nd check is positive
scalar kappaScaleDelta = readScalar(chemistryProperties.lookup("kappaScaleDelta"));

