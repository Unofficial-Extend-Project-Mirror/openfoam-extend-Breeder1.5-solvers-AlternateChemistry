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
    reactingFoam

Description
    Chemical reaction code.

 ICE Revision: $Id: /local/openfoam/Solvers/AlternateChemistry/Steady/alternateSteadyReactingFoam/alternateSteadyReactingFoam.C 4282 2008-12-16T23:01:02.470981Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "compressible/RASModel/RASModel.H"
#include "alternateChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

#   include "readChemistryProperties.H"
#   include "readEnvironmentalProperties.H"

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"
#       include "initConvergenceCheck.H"

        p.storePrevIter();
        rho.storePrevIter();

#       include "chemistry.H"

        // Pressure-velocity SIMPLE corrector
        {
#           include "UEqn.H"
#           include "YEqn.H"
#           include "hEqn.H"
#           include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
