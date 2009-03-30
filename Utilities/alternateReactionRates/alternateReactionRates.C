/*---------------------------------------------------------------------------*\
This file written by Institute of Energy Process Enineering and Chemical
	Engineering TU Freiberg  http://www.iec.tu-freiberg.de
and ICE Stroemungsfoschungs GmbH http://www.ice-sf.at
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
   thermoFieldsCompressible

Description

 ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Utilities/alternateReactionRates/alternateReactionRates.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "wallFvPatch.H"

#include "reactingMixture.H"
#include "chemistryReader.H"
#include "alternateChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "addRegionOption.H"

    argList::validOptions.insert("steady","");
    argList::validOptions.insert("writeTc","");

#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createNamedMesh.H"

#   include "readChemistryProperties.H"

    bool steady=args.options().found("steady");
    bool writeTc=args.options().found("writeTc");

#   include "createFields.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info << "Time = " << runTime.timeName() << endl;

        scalar chemTime=chemistry().solve
            (
                runTime.value() - steadyChemistryDeltaT*runTime.deltaT().value(),
                steadyChemistryDeltaT*runTime.deltaT().value()
            );
        Info << " Characteristic time of chemistry: " << chemTime << endl;

        for(label i=0; i<Y.size(); i++)
        {
            RR[i]=chemistry().RR(i);
            RRsum+=RR[i];
            RR[i].write();
        }
        RRsum.write();

        if(writeTc) {
            Info << "Writing tc" << endl;
            
            volScalarField tc(chemistry().tc());

            tc.write();
        }

        Info << endl;
    }

    Info << "\n Done" << endl;

    return(0);
}


// ************************************************************************* //
