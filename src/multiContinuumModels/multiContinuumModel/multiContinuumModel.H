/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Library
    libmatrixmultiContinuumModels

Description
    Abstract class for multi-continuum models for porous media.

Developers
    Federico Municchi, University of Nottingham (2019)

\*---------------------------------------------------------------------------*/
#ifndef multiContinuumModel_H
#define multiContinuumModel_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

namespace multiContinuumModels
{
class multiContinuumModel
{

private:

    //- Reference to field in the mobile region
    const volScalarField& psiM_;

    //- Reference to capacity in the mobile zone
    const volScalarField& beta_;

    //- Adaptive time stepping
    Switch adaptiveTS_;

    //- Max Courant for adaptive time stepping
    scalar maxCo_;

    //- dictionary
    dictionary dict_;

    //- Adjust time step
    void adjustTimeStep();

protected:

    //- Access functions for derived classes
    const volScalarField& psiM() {return psiM_;};

    const volScalarField& beta() {return beta_;};

    dictionary& dict() {return dict_;};

    //- Returns the maximum current Courant
    virtual scalar maxCurrentCo() = 0;

public:

    //- Runtime type information
    TypeName("multiContinuumModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiContinuumModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& psi,
                const volScalarField& beta
            ),
            (dict, psi, beta)
        );



    // Constructors

        multiContinuumModel
        (
            const dictionary& dict,
            const volScalarField& psi,
            const volScalarField& beta
        );


    //- Destructor
    virtual ~multiContinuumModel();


    // Selectors

        static autoPtr<multiContinuumModel> New
        (
            const dictionary& dict,
            const volScalarField& psi,
            const volScalarField& beta,
            const word& multiContinuumModelType
        );


    // Member Functions

    //- Return scalar source matrix
    virtual tmp<fvScalarMatrix> source();

    //- Correct mass transfer
    virtual void correct();

    //- Update for new time step
    virtual void update();

    //- Read fracturedPorousProperties dictionary
    virtual bool read() const;


    inline const Switch& adaptiveTS() {return adaptiveTS_;};

    inline const scalar& maxCo() {return maxCo_;};
};

}
}

#endif
