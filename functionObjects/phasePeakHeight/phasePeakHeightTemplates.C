/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "phasePeakHeight.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::phasePeakHeight::foundObject
(
    const word& name,
    const bool verbose
) const
{
    if (fvMeshFunctionObject::foundObject<Type>(name))
    {
        return true;
    }
    else
    {
        if (debug || verbose)
        {
            Warning
                << "    functionObjects::" << type() << " " << this->name()
                << " cannot find required object " << name << " of type "
                << Type::typeName << endl;
        }

        return false;
    }
}

template<class Type>
Foam::scalar Foam::functionObjects::phasePeakHeight::findPeakHeight
(
    const word& fieldName,
    const scalar& fieldCrit_
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    scalar height = 0;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);

        // loop over the cells and find the peak height of the phase
        forAll(field, cellI)
        {
            vector pos=mesh_.C()[cellI];
            if((field[cellI]>fieldCrit_) && (pos[2]>height))
            {
                height = pos[2]; 
            }
        }
        
        Info<<"Max height of "<<fieldName<<": "<<height<<endl;
    }
    return height;
}



// ************************************************************************* //
