//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/YEqn.H 4040 2008-10-31T12:22:18.094531Z bgschaid  $ 
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    List<scalarField> Yold(Y.size(),scalarField(kappa.size()));

    for(label i=0; i<Y.size(); i++)
    {
        Yold[i]=Y[i].internalField();
    }

    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    bool recalculate = true;
    label recalc_count= 0;

    while( (recalculate == 1) && recalc_count < recalcMax)
    {
        recalculate = false;
        recalc_count++;

        Yt = 0.0*Y[0];

        for(label i=0; i<Y.size(); i++)
        {
            if (Y[i].name() != inertSpecie)
            {
                volScalarField& Yi = Y[i];

                tmp<fvScalarMatrix> YEqn
                    (
                        - fvm::Sp(fvc::div(phi), Yi)
                        + mvConvection->fvmDiv(phi, Yi)
                        - fvm::laplacian(turbulence->muEff(), Yi)
                        ==
                        kappa*chemistry().RR(i)
                    );

                YEqn().relax(mesh.relaxationFactor("Yi"));

                eqnResidual = solve(YEqn() , mesh.solver("Yi")).initialResidual();
                maxResidual = max(eqnResidual, maxResidual);

                Yi.max(0.0);
                Yt += Yi;
            }
            else
            {
                inertIndex = i;
            }
        }

        Info << "Max summe Yi = " << max(Yt) << " Min summe Yi = " << min(Yt) <<endl;

        //Check Solution if sum of species mass fraction >1 => reduce kappa and recalculate YEqn.
        //This stabilizes the solution process at startup. This should not happen towards the end of the solution!

        //1st check: Yt
        scalar maxYt=max(Yt.internalField());
        reduce(maxYt,maxOp<scalar>());

        if (maxYt > yOverThres)
        {
            label cnt=0;

            recalculate = true;
            Info << "Yt ("<< maxYt <<")> " << yOverThres 
                << " recalculating YEqn! (#"<< recalc_count
                << ")==========================="<< endl;
            for(label celli=0; celli<kappa.size(); celli++)
            {
                if (Yt[celli]>yOverThres)
                {
                    cnt++;

                    Pout << "Yt = "<< Yt[celli] << " !!! for cell " << celli 
                        << ". Setting kappa from " <<	kappa[celli] 
                        << " to " << kappa[celli]*kappaScaleDelta << endl;

		    kappa[celli] = kappa[celli]*kappaScaleDelta;
                }
            }

            reduce(cnt,plusOp<label>());
            Info << " Fixed kappa in " << cnt << " cells due to Yt>" << yOverThres << endl;

        }

        //2nd check: deltaY
        {
            label cnt=0;

            for(label celli=0; celli<kappa.size(); celli++) //calculate difference
            {
                for(label i=0; i<Y.size(); i++)
                {   
                    scalar absDeltaY = fabs(Y[i][celli]-Yold[i][celli]);
                    scalar relDeltaY = absDeltaY/max(relDeltaYCalcThres,Yold[i][celli]);


                    if ( (absDeltaY>minAbsDeltaY) && 
                    ((absDeltaY>absDeltaYThres) || (relDeltaY>relDeltaYThres)) )
                    {       
                        recalculate = true;
                        cnt++;

                        Pout << "diffY to big!  = absDeltaY="<< absDeltaY 
                            << " relDeltaY="<< relDeltaY <<" for cell " 
                            << celli << ". Setting kappa from " << kappa[celli] 
                            << " to " << kappa[celli]*kappaScaleDelta << " count = "
                            << recalc_count <<endl;

			kappa[celli] = kappa[celli]*kappaScaleDelta;
                        break; //only once for every cell
                    }
                }
            }
            
            reduce(recalculate,orOp<bool>());
            reduce(cnt,plusOp<label>());
            Info << " Fixed kappa in " << cnt << " cells due to large deltaY" << endl;
        }

        if ( (recalculate) && (recalc_count<recalcMax))
        {
            for(label i=0; i<Y.size(); i++)
            {
                Y[i].internalField()=Yold[i];
            }
        }
    }

    if(recalc_count>1) {
        Info << "Recalculated chemistry " << recalc_count << " of " 
            << recalcMax << " possible times" << endl;
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}

