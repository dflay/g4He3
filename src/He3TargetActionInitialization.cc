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
/// \file He3TargetActionInitialization.cc
/// \brief Implementation of the He3TargetActionInitialization class

#include "He3TargetActionInitialization.hh"
#include "He3TargetPrimaryGeneratorAction.hh"
#include "He3TargetRunAction.hh"
#include "He3TargetEventAction.hh"
#include "He3TargetSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3TargetActionInitialization::He3TargetActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3TargetActionInitialization::~He3TargetActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3TargetActionInitialization::BuildForMaster() const
{
  He3TargetRunAction* runAction = new He3TargetRunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3TargetActionInitialization::Build() const
{
  SetUserAction(new He3TargetPrimaryGeneratorAction);

  He3TargetRunAction* runAction = new He3TargetRunAction;
  SetUserAction(runAction);
  
  He3TargetEventAction* eventAction = new He3TargetEventAction(runAction);
  SetUserAction(eventAction);
  
  SetUserAction(new He3TargetSteppingAction(eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
