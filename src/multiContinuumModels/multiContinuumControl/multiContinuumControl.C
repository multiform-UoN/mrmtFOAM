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

\*---------------------------------------------------------------------------*/

#include "multiContinuumControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiContinuumControl, 0);
}

using namespace multiContinuumModels;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiContinuumControl::multiContinuumControl
(
    fvMesh& mesh,
    const volScalarField& psiM,
    const volScalarField& beta,
    const word& algorithmName
)
:
    solidNoLoopControl(mesh, algorithmName, *this),
    pimpleLoop(static_cast<solutionControl&>(*this))
{
    read();

    printResidualControls();

    if (nCorrPimple_ > 1)
    {
        printCorrResidualControls(nCorrPimple_);
    }
    else
    {
        Info<< nl << algorithmName << ": Operating solver in PISO mode" << nl
            << endl;
    }

    Info<< nl << "Reading multiContinuumProperties" <<endl;
    IOdictionary multiContinuumProperties
    (
        IOobject
        (
            "multiContinuumProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    mcModel_.set
    (
        multiContinuumModel::New
        (
            multiContinuumProperties,
            psiM,
            beta,
            word(multiContinuumProperties.lookup("Type"))
        ).ptr()

    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiContinuumControl::~multiContinuumControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::multiContinuumControl::read()
{
    if (!solidNoLoopControl::read() || !pimpleLoop::read())
    {
        return false;
    }

    nCorrPimple_ = dict().lookupOrDefault<label>("nOuterCorrectors", 1);

    return true;
}


bool Foam::multiContinuumControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        mesh().data::remove("finalIteration");

        return false;
    }

    storePrevIterFields();

    if (finalPimpleIter())
    {
        mesh().data::add("finalIteration", true);
    }

    return true;
}


bool Foam::multiContinuumControl::run(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.run();
}


bool Foam::multiContinuumControl::loop(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.loop();
}

tmp<fvScalarMatrix> Foam::multiContinuumControl::immobileRegionsSource()
{
    return mcModel_->source();
}

void Foam::multiContinuumControl::correctImmobileRegions()
{
    mcModel_->correct();
}

void Foam::multiContinuumControl::updateMultiContinuumModel()
{
    mcModel_->update();
}

// ************************************************************************* //
