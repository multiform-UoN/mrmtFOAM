/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiRateScalarTransportFoam

Description
    Multi-rate solver for advection-diffusion of a scalar field. Coupling is
    obtained using pimple.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiContinuumModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createTimeControls.H"
#   include "createMesh.H"

#   include "createFields.H"

    pimpleControl pimple(mesh);

    autoPtr<multiContinuumModels::multiContinuumModel> multiRate
    (
        multiContinuumModels::multiContinuumModel::New
        (
            multiRateProperties,
            c,
            beta,
            "multiRateMassTransfer"
        ).ptr()
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"
#   include "constrainDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
   
	    multiRate->update();
        runTime++;

	Info<< "Time = " << runTime.timeName() << nl << endl;

        //- solve the coupled multi-continuum equations
        //  until convergence criteria are met
        while (pimple.loop())
        {

            //- update concentration in the immobile regions
            multiRate->correct();

            //- more iterations are required for explicit fluxes
            //  arising, for example, from non-orthogonal grids
            while (pimple.correctNonOrthogonal())
            {
                //- solve scalar transport for c
                fvScalarMatrix cEqn
                (
                    fvm::ddt(beta,c)
                  + multiRate->source()
                  + fvm::div(phi, c)
                  - fvm::laplacian(Dc, c)
                );
                //- relax and solve the scalar transport equation
                cEqn.relax();
                cEqn.solve();
            }

        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
