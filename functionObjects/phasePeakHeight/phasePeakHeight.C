/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "volFields.H"

#include "memInfo.H"
#include "OFstreamMod.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

#define NUMBER_OF_COLUMNS 10

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(phasePeakHeight, 0);
    addToRunTimeSelectionTable(functionObject, phasePeakHeight, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::phasePeakHeight::printInfo() const
{
    Info<<"***************************************************************"<<nl;
    Info<<"Dictionary parameters"<<nl;
    Info<< "field names: " << fieldNames_ << nl;
    Info<< "phase crit: " << phaseCrit_ << nl;
        
    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;
    Info<<"***************************************************************"<<nl;
}

bool Foam::functionObjects::phasePeakHeight::cellInsideTheBox(point& pos) const
{
    bool res = true;
    
    for(int i=0; i<3; i++)
    {
        res = res  &&  ( pos[i]>=minPoint_[i] && pos[i]<=maxPoint_[i] );
    }
    
    return res;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::phasePeakHeight::phasePeakHeight
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_),
    height_(0.0)
{
    if (debug)
    {
        Info<< "functionObjects::phasePeakHeight Constructor"<<nl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::phasePeakHeight::~phasePeakHeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::phasePeakHeight::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("phasePeakHeight::read")
              << "There is no fields parameter in phasePeakHeight dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<scalar>("phaseCrit", phaseCrit_) ){
        SeriousErrorIn("phasePeakHeight::read")
              << "There is no phaseCrit parameter in phasePeakHeight dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<point>("minPoint", minPoint_) ){
        SeriousErrorIn("phasePeakHeight::read")
              << "There is no minPoint parameter in phasePeakHeight dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("maxPoint", maxPoint_) ){
        SeriousErrorIn("phasePeakHeight::read")
              << "There is no maxPoint parameter in phasePeakHeight dictionary"
              << exit(FatalError);
    }
    
    fieldSet_.read(dict);
    
    Info<<"END READ"<<nl<<nl;
    
    printInfo();
    
    Info<<nl<<nl;
    
    return true;
}


bool Foam::functionObjects::phasePeakHeight::execute()
{
    if(debug)
    {
        Info<<"Integrate fields within the box boundaries"<<nl;
        printInfo();
    }

    scalar height = 0;
    
    //for (const word& fieldName : fieldSet_.selection())
    for (const word& fieldName : fieldNames_)
    {
        height = findPeakHeight<scalar>(fieldName, phaseCrit_);
    }

    height_ = height;

    write();
    
    return true;
}

bool Foam::functionObjects::phasePeakHeight::write()
{
    Info<<"Currently there is no write into a file"<<nl;
    string fileName = "phasePeakHeight.csv";
    std::ifstream file(fileName);
    std::fstream fout;
    fout.open(fileName, std::ios::out | std::ios::app);
    if (file.peek() == std::ifstream::traits_type::eof())
    {
        fout << "Time[s],totSufArea,averNuRate" << "\n";
    }
    fout << mesh_.time().timeName() << ", "
        << height_ << "\n";
    return true;
}


// ************************************************************************* //
