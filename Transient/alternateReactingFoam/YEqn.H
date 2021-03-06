//  ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Transient/alternateReactingFoam/YEqn.H 3602 2008-08-21T16:17:37.921551Z bgschaid  $ 
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
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    for(label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            solve
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                kappa*chemistry().RR(i),
                mesh.solver("Yi")
            );

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }
    
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
