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
#include "FirstOrderRegions.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    namespace immobileRegions
    {
        defineTypeNameAndDebug(FirstOrderRegions, 0);
        addToRunTimeSelectionTable
        (
            immobileRegion,
            FirstOrderRegions,
            dictionary
        );
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immobileRegions::FirstOrderRegions::FirstOrderRegions
(
    word  name,
    dictionary dict,
    const volScalarField& psiM
)
:
    immobileRegion(name,dict,psiM)
{
    computeCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immobileRegions::FirstOrderRegions::~FirstOrderRegions()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immobileRegions::FirstOrderRegions::computeCoeffs()
{

    scalarField alphas("alphaCoeffs",dict(),nOfTerms());
    scalarField betas("betaCoeffs",dict(),nOfTerms());

    alphaCoeffs_ = alphas;
    betaCoeffs_ = betas;

    betaCoeffs_/=gSum(betaCoeffs_);

    if(dict().found("debug"))
    {
        scalar sumBetas(0);

        Info<<"\nDebug from region " << name() << endl;
        for(int term=0;term<alphaCoeffs_.size();term++)
        {
            Info<<"alpha_" << term << ": "<<alphaCoeffs_[term]<<endl;
            Info<<"beta_" << term << ": "<<betaCoeffs_[term]<<endl;
            sumBetas+=betaCoeffs_[term];
        }
        Info<<"Sum of betas: " << sumBetas << endl;

    }
}
