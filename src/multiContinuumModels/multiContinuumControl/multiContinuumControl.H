/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::multiContinuumControl

Description
    Multi-continuum control class. As a PIMPLE, it provides time-loop control
    methods which exit the simulation once convergence criteria have been
    reached. Also provides Pimple-loop control methods which exit the iteration
    once corrector convergence criteria have been met.
    Furthermore, it stores the multi-continuum model.

Author
    Federico Municchi, University of Nottingham (2019)

\*---------------------------------------------------------------------------*/

#ifndef multiContinuumControl_H
#define multiContinuumControl_H

#include "solidNoLoopControl.H"
#include "pimpleLoop.H"
#include "singleRegionConvergenceControl.H"
#include "singleRegionCorrectorConvergenceControl.H"
#include "multiContinuumModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class multiContinuumControl Declaration
\*---------------------------------------------------------------------------*/

class multiContinuumControl
:
    public solidNoLoopControl,
    public pimpleLoop
{
private:

    //- Multi continuum model
    autoPtr<multiContinuumModels::multiContinuumModel> mcModel_;

public:

    // Static data members

        //- Run-time type information
        TypeName("multiContinuumControl");


    // Constructors

        //- Construct from a mesh, the field in the mobile region
        //  and the name of the algorithm
        multiContinuumControl
        (
            fvMesh& mesh,
            const volScalarField& psiM,
            const volScalarField& beta,
            const word& algorithmName="MultiContinuumCoupling"
        );


    //- Destructor
    virtual ~multiContinuumControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();


        // Evolution

            //- Pimple loop
            bool loop();

            //- Return source term for the psi equation
            tmp<fvScalarMatrix> immobileRegionsSource();

            //- Correct multi-continuum models
            void correctImmobileRegions();

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);

            //- Update multi-continuum model
            void updateMultiContinuumModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
