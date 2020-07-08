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
    libmatriximmobileRegions

Description
    Base class for immobile regions.

Developers
    Federico Municchi, University of Nottingham (2019)

\*---------------------------------------------------------------------------*/
#ifndef immobileRegion_H
#define immobileRegion_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

namespace immobileRegions
{
class immobileRegion
{

private:

    //- Name (external) of the immobile region
    word    name_;

    //- Number of terms in the series expansion
    label nOfTerms_;

    //- dictionary
    dictionary dict_;

    //- Scalar fields in the immobile regions
    PtrList<volScalarField> psiImms_;

	//- Avrage scalar field
	volScalarField psiImm_;

	//- Region capacity
	volScalarField beta_;

    //- Region exchange frequency
	volScalarField omega_;

    //- Rescale betas or not?
    Switch rescaleBetas_;

    //- Solve for one immobile region
    void solveImmobileRegion
    (
        const volScalarField& psiMobile,
        volScalarField&       psiImmobile,
        dimensionedScalar&    alphaI,
        dimensionedScalar     deltaT
    );

protected:

    //- Alpha coefficient for each term in the series
    scalarField alphaCoeffs_;

    //- Beta coefficient for each term in the series
    scalarField betaCoeffs_;

    dictionary& dict() {return dict_;};

    virtual void computeCoeffs() {};

public:

    //- Runtime type information
    TypeName("immobileRegion");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            immobileRegion,
            dictionary,
            (
                word  name,
                dictionary dict,
                const volScalarField& psiM
            ),
            (name, dict, psiM)
        );



    // Constructors

        immobileRegion
        (
            word  name,
            dictionary dict,
            const volScalarField& psiM
        );


    //- Destructor
    virtual ~immobileRegion();


    // Selectors

        static autoPtr<immobileRegion> New
        (
            word  name,
            dictionary dict,
            const volScalarField& psiM
        );


    // Member Functions

    //- Correct mass transfer
    void correct(const volScalarField& psiM) ;

    //- Access functions for derived classes
    inline const scalarField& alphaCoeffs() const
    {
        return alphaCoeffs_;
    };

    inline const scalarField& betaCoeffs() const
    {
        return betaCoeffs_;
    };

    inline const PtrList<volScalarField>&
    psiImms() const {return psiImms_;};

    inline const word name() const  {return name_;};

    inline label nOfTerms() {return nOfTerms_;};

    inline const volScalarField&  psiImm() {return psiImm_;}

    inline const volScalarField&  beta() {return beta_;};

    inline const volScalarField&  omega() {return omega_;};

    inline const Switch& rescaleBetas() {return rescaleBetas_;};

    //- Return maximum reactive Courant correspoding to the actual TS
    scalar currentMaxCo();

    //- Return truncation error on betas
    scalar truncBetaError() const;

};

}
}

#endif
