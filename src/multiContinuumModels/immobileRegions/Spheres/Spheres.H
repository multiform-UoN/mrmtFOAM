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
    Class for spherical immobile regions for multi-rate models.

Developers
    Federico Municchi, University of Nottingham (2019)

\*---------------------------------------------------------------------------*/
#ifndef immobileSpheres_H
#define immobileSpheres_H

#include "fvCFD.H"
#include "immobileRegion.H"
#include "autoPtr.H"
#include "dictionary.H"

namespace Foam
{

namespace immobileRegions
{
class Spheres
:
public immobileRegion
{

private:


protected:

    //- Compute coefficients of the series
    void computeCoeffs();

public:

    //- Runtime type information
    TypeName("Spheres");


    // Constructors

    Spheres
    (
        word  name,
        dictionary dict,
        const volScalarField& psiM
    );


    //- Destructor
    virtual ~Spheres();


};

}
}

#endif
