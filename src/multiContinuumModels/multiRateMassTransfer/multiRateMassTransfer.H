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
    libmatrixmultiRateMassTransferMassTransfer

Description
    Multi-rate mass transfer model as from Haggerty & Gorelick (1995). In this
    model, immobile regions are modeled as additional scalar fields coupled to
    the mobile regions by mean of first order reaction.

Reference
    Haggerty, R., & Gorelick, S. M. (1995). Multiple-Rate Mass Transfer for
    Modeling Diffusion and Surface Reactions in Media with Pore-Scale
    Heterogeneity. Water Resources Research, 31(10), 2383–2400.
    https://doi.org/10.1029/95WR10583

Developers
    Federico Municchi, University of Nottingham (2019)

\*---------------------------------------------------------------------------*/
#ifndef multiRateMassTransfer_H
#define multiRateMassTransfer_H

#include "multiContinuumModel.H"
#include "immobileRegion.H"

namespace Foam
{

using namespace immobileRegions;

namespace multiContinuumModels
{

class multiRateMassTransfer
:
public multiContinuumModel
{
private:

    //- Immobile regions
    PtrList<immobileRegion> immobileRegions_;

    //- Total field value (including immo)
    volScalarField psiAve_;

    //- Check consistency of input parameters
    void checkConsistency();

protected:

    virtual scalar maxCurrentCo();

public:

    TypeName("multiRateMassTransfer");

    multiRateMassTransfer
    (
        const dictionary& dict,
        const volScalarField& psiM,
        const volScalarField& beta
    );

    ~multiRateMassTransfer();

    //- Return scalar source matrix
    tmp<fvScalarMatrix> source();

    //- Correct mass transfer
    void correct();

    //- Read fracturedPorousProperties dictionary
    bool read() const {return false;};

};
}
}
#endif
