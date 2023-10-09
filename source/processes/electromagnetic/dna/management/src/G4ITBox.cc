//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITBox.hh"

G4ITBox::G4ITBox() 
= default;

G4ITBox::~G4ITBox()
{
    if( fNbIT != 0 )
    {
        G4IT * aIT = fpFirstIT;
        G4IT * nextIT;

        while( aIT != nullptr )
        {
            nextIT = aIT->GetNext();
            delete aIT;
            aIT = nextIT;
        }
    }

    if(fpPreviousBox != nullptr)    fpPreviousBox->SetNextBox(fpNextBox) ;
    if(fpNextBox != nullptr)        fpNextBox->SetPreviousBox(fpPreviousBox);
}

const G4ITBox & G4ITBox::operator=(const G4ITBox &right)
{
    fNbIT = right.fNbIT;
    fpFirstIT = right.fpFirstIT;
    fpLastIT = right.fpLastIT;
    fpPreviousBox = nullptr;
    fpNextBox = nullptr;
    return *this;
}

void G4ITBox::Push( G4IT * aIT )
{
    if( fNbIT == 0 )
    {
        aIT->SetPrevious( nullptr );
        fpFirstIT = aIT;
    }
    else
    {
        fpLastIT->SetNext( aIT );
        aIT->SetPrevious( fpLastIT );
    }
    fpLastIT = aIT;
    fNbIT++;
    aIT->SetITBox(this);
}

void G4ITBox::Extract( G4IT * aStackedIT )
{
    if( aStackedIT == fpFirstIT )
    {
        fpFirstIT = aStackedIT->GetNext();
    }
    else  if( aStackedIT == fpLastIT )
    {
        fpLastIT = aStackedIT->GetPrevious();

    }

    if( aStackedIT->GetNext() != nullptr)
        aStackedIT->GetNext()->SetPrevious(aStackedIT->GetPrevious());
    if( aStackedIT->GetPrevious() != nullptr)
        aStackedIT->GetPrevious()->SetNext(aStackedIT->GetNext());

    aStackedIT->SetNext(nullptr);
    aStackedIT->SetPrevious(nullptr);
    aStackedIT->SetITBox(nullptr);
    fNbIT--;
}

G4IT* G4ITBox::FindIT(const G4Track& track)
{
    if( fNbIT == 0 ) return nullptr;

    G4IT * temp = fpLastIT;
    G4bool find = false;

    while(!find && temp != nullptr)
    {
        if(temp-> GetTrack() == &track)
        {
            find = true;
            break;
        }
        temp = temp->GetPrevious();
    }

    return temp;
}

const G4IT* G4ITBox::FindIT(const G4Track& track) const
{
    if( fNbIT == 0 ) return nullptr;

    const G4IT * temp = fpLastIT;
    G4bool find = false;

    while(!find && temp != nullptr)
    {
        if(temp-> GetTrack() == &track)
        {
            find = true;
            break;
        }
        temp = temp->GetPrevious();
    }

    return temp;
}

void G4ITBox::TransferTo(G4ITBox * aStack)
{
    G4IT * ITToTransfer = fpFirstIT;
    while(fNbIT != 0)
    {
        Extract(ITToTransfer);
        aStack->Push(ITToTransfer);
        ITToTransfer = ITToTransfer->GetNext();
    }
}
