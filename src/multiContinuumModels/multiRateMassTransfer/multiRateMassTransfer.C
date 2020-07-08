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
#include "multiRateMassTransfer.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    namespace multiContinuumModels
    {
        defineTypeNameAndDebug(multiRateMassTransfer, 0);
        addToRunTimeSelectionTable
        (
            multiContinuumModel,
            multiRateMassTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiContinuumModels::multiRateMassTransfer::multiRateMassTransfer
(
    const dictionary& dict,
    const volScalarField& C,
    const volScalarField& beta
)
:
    multiContinuumModel(dict, C, beta),
    psiAve_
    (
        IOobject
        (
            C.name() + ".total",
            C.mesh().time().timeName(),
            C.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C
    )
{
    //- Read immobile regions
    PtrList<entry> regions(dict.lookup("immobileRegions"));
    immobileRegions_.setSize(regions.size());

    forAll(immobileRegions_,immId)
    {
        immobileRegions_.set
        (
            immId,
            immobileRegion::New
            (
                regions[immId].keyword(),
                regions[immId].dict(),
                psiM()
            ).ptr()
        );
    }

    //- Check that input is consistent
    checkConsistency();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiContinuumModels::multiRateMassTransfer::~multiRateMassTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::multiContinuumModels::multiRateMassTransfer::source()
{
    Foam::tmp<Foam::fvScalarMatrix> tfvsm
    (
        new fvScalarMatrix
        (
            psiM(),
            psiM().dimensions()/dimTime*dimVol
        )
    );


    const volScalarField& betaM = beta();

    fvScalarMatrix& fvsm = tfvsm.ref();

    //- Add contributions to the source term
    forAll(immobileRegions_,immId)
    {

        const scalarField& alphaCoeffs
        (
            immobileRegions_[immId].alphaCoeffs()
        );

        const scalarField& betaCoeffs
        (
            immobileRegions_[immId].betaCoeffs()
        );

        const PtrList<volScalarField>& psiImms
        (
            immobileRegions_[immId].psiImms()
        );


	    const volScalarField& betaI = immobileRegions_[immId].beta();
        const volScalarField& omegaI = immobileRegions_[immId].omega();

        //- Add contributions from multi-rate to matrix
        forAll(psiImms, rateTerm)
        {
            fvsm += omegaI*betaI*
            (
                  fvm::Sp
                  (
                      alphaCoeffs[rateTerm]*betaCoeffs[rateTerm],
                      psiM()
                  )
                - (
                      psiImms[rateTerm]()
                      *alphaCoeffs[rateTerm]
                      *betaCoeffs[rateTerm]
                  )
            );

        }

        //- Contribution from truncation errors (equilibrium modes)
        fvsm += betaI*immobileRegions_[immId].truncBetaError()
                * fvm::ddt(psiM());

    }

    fvsm *= betaM;

    return tfvsm;

}

void
Foam::multiContinuumModels::multiRateMassTransfer::correct()
{
    psiAve_=psiM();

    volScalarField betaTot(beta());


    //- Correct psi in the immobile regions
    forAll(immobileRegions_,immId)
    {
	    const volScalarField& betaI =
			immobileRegions_[immId].beta();

        const volScalarField& psiImm
        (
            immobileRegions_[immId].psiImm()
        );

        immobileRegions_[immId].correct(psiM());
        psiAve_+=betaI*psiImm;
		betaTot+=betaI;
    }

    psiAve_/=betaTot;
}

Foam::scalar
Foam::multiContinuumModels::multiRateMassTransfer::maxCurrentCo()
{
    scalar currentCo(0.);

    forAll(immobileRegions_,immId)
    {
        if
        (
            immobileRegions_[immId].currentMaxCo() > currentCo
        )
        {
            currentCo = immobileRegions_[immId].currentMaxCo();
        }
    }

    return currentCo;
}

void
Foam::multiContinuumModels::multiRateMassTransfer::checkConsistency()
{
    //- Check that sum of the capcities do not exceed 1
    scalar totBeta(1.0);

    scalar maxBetam(gMax(beta()));
    reduce(maxBetam,maxOp<scalar>());

    forAll(immobileRegions_,immId)
    {
        scalar maxBetaId(gMax(immobileRegions_[immId].beta()));
        reduce(maxBetaId,maxOp<scalar>());
        totBeta += maxBetaId;
    }

    totBeta *= maxBetam;

    if(totBeta > 1.0)
    {
        FatalErrorInFunction
        <<"You entered a maximum total capacity equal to " << totBeta
        <<", which is larger than one!" << endl
        <<"Check your input carefully and remember that  "
        <<"capacities in the immobile regions are defined realtive " 
        <<"to the mobile regions." << endl
        <<"Therefore your capacities should satisfy: " << endl
        <<"beta_m * ( 1.0 + beta_1 + beta_2 + ... ) <= 1"
        << abort(FatalError);        
    }
}