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

Developers
    Federico Municchi, University of Nottingham (2019)

\*---------------------------------------------------------------------------*/
#include "multiContinuumModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace multiContinuumModels
    {
        defineTypeNameAndDebug(multiContinuumModel, 0);
        defineRunTimeSelectionTable(multiContinuumModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::multiContinuumModels::multiContinuumModel::multiContinuumModel
(
  const dictionary& dict,
  const volScalarField& psi,
  const volScalarField& beta
)
:
    psiM_(psi),
    beta_(beta),
    adaptiveTS_(dict.lookupOrDefault<bool>("adaptiveTimeStepping",false)),
    maxCo_(1.0)
{

    if(adaptiveTS())
    {
        maxCo_ = readScalar(dict.lookup("maxCo"));
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiContinuumModels::multiContinuumModel::~multiContinuumModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<fvScalarMatrix>
Foam::multiContinuumModels::multiContinuumModel::source()
{
    tmp<fvScalarMatrix> tSource
    (
        new fvScalarMatrix
        (
            psiM_,
            psiM_.dimensions()/dimTime*dimVol
        )
    );

    return tSource;
}

void
Foam::multiContinuumModels::multiContinuumModel::correct()
{}

void
Foam::multiContinuumModels::multiContinuumModel::update()
{
    if(adaptiveTS())
    {
        adjustTimeStep();
    }
}

bool
Foam::multiContinuumModels::multiContinuumModel::read() const
{
    return false;
}

void
Foam::multiContinuumModels::multiContinuumModel::adjustTimeStep()
{
    scalar deltaT(psiM().mesh().time().deltaT().value());

    scalar currCo(maxCurrentCo());

    Info<<"\nCurrent reaction Co: " << currCo << " maximum reaction Co: "
        << maxCo() << endl;

    if(currCo > maxCo() + small)
    {

        scalar newDeltaT(min(deltaT/currCo*maxCo(),deltaT));

        if(adaptiveTS())
        {
            Info << "Setting new time step to " << newDeltaT << " s" << endl;

            Time& runTime(const_cast<Time&>(psiM().mesh().time()));
            runTime.setDeltaT(newDeltaT);
        }
        else
        {
            WarningInFunction
                << "\nMaximum reaction Courant number is Co = " << currCo
                << " larger than maxCo = "<< maxCo() <<" !" << endl
                << "Try reducing your time step below " << newDeltaT
                << " s or set : \n" << endl
                << "adaptiveTimeStepping true; \n" << endl
                << "in " << dict().name() << endl;
        }
    }

}

// ************************************************************************* //
