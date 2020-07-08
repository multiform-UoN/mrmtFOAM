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
#include "Cylinders.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

#define numberOfZeros 31831

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    namespace immobileRegions
    {
        defineTypeNameAndDebug(Cylinders, 0);
        addToRunTimeSelectionTable
        (
            immobileRegion,
            Cylinders,
            dictionary
        );
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immobileRegions::Cylinders::Cylinders
(
    word  name,
    dictionary dict,
    const volScalarField& psiM
)
:
    immobileRegion(name,dict,psiM),
    zerosJ0_(numberOfZeros)
{
    #include "Jzeros.H"

    computeCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immobileRegions::Cylinders::~Cylinders()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immobileRegions::Cylinders::computeCoeffs()
{

    if (alphaCoeffs_.size()>zerosJ0_.size())
    {
        FatalErrorInFunction
            <<"You can not have more than " << numberOfZeros
            <<" expansion terms for cylindrical immobile regions!" << endl
            <<"I simply do not have enough zeros of J0 "
            <<"to fullfill your request" << endl
            <<"Why would you need so many expansion terms anyway?"
            << abort(FatalError);
    }

    Info << nl << "Compute series expansion terms for region "
         << name() <<endl;

    //- Compute series expansion
    for(int term=0;term<alphaCoeffs_.size();term++)
    {
        alphaCoeffs_[term] = sqr(zerosJ0_[term]);
        betaCoeffs_[term] = scalar(4.)/sqr(zerosJ0_[term]);
    }

    if(rescaleBetas())
	{
        betaCoeffs_/=gSum(betaCoeffs_);
    }

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
