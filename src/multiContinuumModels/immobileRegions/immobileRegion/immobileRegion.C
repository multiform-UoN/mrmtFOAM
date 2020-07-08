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
#include "immobileRegion.H"
#include "subCycle.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    namespace immobileRegions
    {
        defineTypeNameAndDebug(immobileRegion, 0);
        defineRunTimeSelectionTable(immobileRegion, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immobileRegions::immobileRegion::immobileRegion
(
    word  name,
    dictionary dict,
    const volScalarField& psiM
)
:
    name_(name),
    nOfTerms_(readLabel(dict.lookup("numberOfTermsInExpansion"))),
    dict_(dict),
    psiImms_(nOfTerms_),
    psiImm_
    (
        IOobject
        (
            psiM.name() + "." + name,
            psiM.mesh().time().timeName(),
            psiM.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        psiM.mesh()
    ),
    beta_
    (
        IOobject
        (
            "beta." + name,
            psiM.mesh().time().timeName(),
            psiM.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        psiM.mesh()
    ),
    omega_
    (
        IOobject
        (
            "omega." + name,
            psiM.mesh().time().timeName(),
            psiM.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        psiM.mesh()
    ),
    rescaleBetas_(dict.lookupOrDefault<bool>("rescaleBetas",false)),
    alphaCoeffs_(nOfTerms_),
    betaCoeffs_(nOfTerms_)
{
    //- Initalise immobile fields
    Info<<nl<<"Initialising immobile fields for region " <<name<<endl;
    forAll(psiImms_,term)
    {
        psiImms_.set
        (
            term,
            new volScalarField
            (
                IOobject
                (
                    psiM.name() + "_" + name + "_term" + std::to_string(term),
                    psiM.mesh().time().timeName(),
                    psiM.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                psiImm_
            )
        );

    //    psiImm_[term].storeOldTime();
    }


    //- Compute coefficients of the series
    computeCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immobileRegions::immobileRegion::~immobileRegion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::immobileRegions::immobileRegion::correct
(
    const volScalarField& psiM
)
{

    //- Reset total fields
    psiImm_*=scalar(0);
    scalar sumBetaI(0);

    dimensionedScalar deltaT(omega().mesh().time().deltaT());

    //- Solve equations for immobile regions
    Info<<"Solving equations for immobile region : "<< name()<<endl;

    solverPerformance::debug=0;

    forAll(psiImms_,term)
    {

        volScalarField& psiI = psiImms_[term];
        dimensionedScalar alphaI
        (
            "alphaI",
            dimless,
            alphaCoeffs_[term]
        );

        dimensionedScalar betaI
        (
            "betaI",
            dimless,
            betaCoeffs_[term]
        );


        solveImmobileRegion(psiM,psiI,alphaI,deltaT);

        //- Add to total
        psiImm_+=betaI*psiI;

        sumBetaI += betaI.value();
    }

    psiImm_/=sumBetaI;

    solverPerformance::debug=1;

}

Foam::scalar
Foam::immobileRegions::immobileRegion::currentMaxCo()
{
    //- Get the effective max reaction rate
    scalar maxOmega(gMax(omega()));

    //- Get maximum from all processors
    reduce(maxOmega,maxOp<scalar>());

    //- Calculate max reaction rate for each term in the expansion
    scalarField maxK(maxOmega*alphaCoeffs_);

    //- Get current deltaT
    scalar deltaT(omega().mesh().time().deltaT().value());

    return deltaT*max(maxK);

}

void
Foam::immobileRegions::immobileRegion::solveImmobileRegion
(
    const volScalarField& psiMobile,
    volScalarField&       psiImmobile,
    dimensionedScalar&    alphaI,
    dimensionedScalar     deltaT
)
{

    //- backward discretization
    // scalar coefft(3.0/2.0);
    // scalar coefft0(-2.0);
    // scalar coefft00(0.5);
    //
    // fvScalarMatrix ddtPsiI
    // (
    //         fvm::Sp
    //         (
    //             coefft/deltaT + omega_()*alphaI,
    //             psiImmobile
    //         )
    //      + (coefft0/deltaT)*psiImmobile.oldTime()
    //      + (coefft00/deltaT)*psiImmobile.oldTime().oldTime()
    // );

    solve
    (
        fvm::ddt(psiImmobile)
      + fvm::Sp(omega()*alphaI,psiImmobile)
      ==
        omega()*alphaI*psiMobile
    );

}

scalar
Foam::immobileRegions::immobileRegion::truncBetaError() const
{
    if(rescaleBetas_)
    {
        return scalar(0);
    }
    else
    {
        return scalar(1.0) - gSum(betaCoeffs());
    }
}
