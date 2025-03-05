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
/// \file hadronic/Hadr00/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
//
// Author:      V.Ivanchenko 30/01/01
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"

#include "HistoManagerMessenger.hh"

#include "G4HadronicProcessStore.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4StableIsotopes.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4IonTable.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fAnalysisManager = 0;
  fHistoName = "hadr00";

  fNeutron = G4Neutron::Neutron();
  fMessenger = new HistoManagerMessenger(this);
  fVerbose = 1;

  fParticleName = "proton";
  fElementName = "G4_Bi";
  fMaterialName = "G4_BGO";
  fIonZ = 6;
  fIonA = 12;
  fIonCharge = 0;


  fTargetMaterial = 0;

  fMinKinEnergy = 1 * GeV;
  fMaxKinEnergy = 1 * TeV;
  fMinMomentum = 1 * GeV;
  fMaxMomentum = 1 * TeV;

  fBinsE = 60;
  fBinsP = 60;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::BeginOfRun()
{
  G4double p1 = std::log10(fMinMomentum / GeV);
  G4double p2 = std::log10(fMaxMomentum / GeV);
  G4double e1 = std::log10(fMinKinEnergy / MeV);
  G4double e2 = std::log10(fMaxKinEnergy / MeV);

  // G4cout<<"e1= "<<e1<<" e2= "<<e2<<" p1= "<<p1<<" p2= "<<p2<<G4endl;

  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->OpenFile(fHistoName + ".root");
  fAnalysisManager->SetFirstHistoId(1);

  fAnalysisManager->CreateH1("h1", "Inelastic interaction length ;log10(E_k/MeV); Length(cm)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h2", "Elastic interaction length ;log10(E_k/MeV); Length(cm)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h3", "Total Hadronic interaction length ;log10(E_k/MeV); Length(cm)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h4", "Inelastic cross section per volume ;log10(E_k/MeV); Cross Section(barn)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h5", "Elastic cross section per volume ;log10(E_k/MeV); Cross Section(barn)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h6", "Total Hadronic cross section per volume ;log10(E_k/MeV); Cross Section(barn)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h7", "Inelastic cross section per volume ;log10(E_kn /MeV); Cross Section(barn)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h8", "Elastic cross section per volume ;log10(E_kn/MeV); Cross Section(barn)", fBinsE, e1, e2);
  fAnalysisManager->CreateH1("h9", "Total Hadronic cross section per volume ;log10(E_kn/MeV); Cross Section(barn)", fBinsE, e1, e2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::EndOfRun()
{
  if (fVerbose > 0) {
    G4cout << "HistoManager: End of run actions are started" << G4endl;
  }

  const G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName);
  const G4Element* elm =   (mat->GetNumberOfElements() == 1) ? mat->GetElement(0) : nullptr;
  if (elm) {
    G4cout << "This material is a single element: " << elm->GetName() << G4endl;
  } else {
      G4cout << "This material is a compound." << G4endl;
  }
  const G4ParticleDefinition* particle;
  if (fParticleName == "ion")
  {
    particle = G4IonTable::GetIonTable()->FindIon(fIonZ, fIonA, 0);
  }
  else
  {
    particle = G4ParticleTable::GetParticleTable()->FindParticle(fParticleName);
  }

  G4int A = 0;
  G4int Z = 0;
  if (particle) 
  {
    if (particle->GetParticleType() == "nucleus") {
      A = particle->GetAtomicMass();
      Z = particle->GetAtomicNumber();
    } 
    else if (particle->GetParticleName() == "proton") {
      A = 1;
      Z = 1;
    } 
    else if (particle->GetParticleName() == "neutron") {
        A = 1;
        Z = 0;
    } 
    else if (particle->GetParticleName() == "alpha") {
        A = 4;
        Z = 2;
    }
    else if (particle->GetParticleName() == "He3") {
        A = 3;
        Z = 2;
    }
    G4cout << "Particle: " << particle->GetParticleName()  << ", A = " << A << ", Z = " << Z << G4endl; 
  } 
  else G4cout << "Error: Particle not found!" << G4endl;


  G4int prec = G4cout.precision();
  G4cout.precision(6);

  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  G4double mass = particle->GetPDGMass();

  // Build histograms

  G4double p1 = std::log10(fMinMomentum / GeV);
  G4double p2 = std::log10(fMaxMomentum / GeV);
  G4double e1 = std::log10(fMinKinEnergy / MeV);
  G4double e2 = std::log10(fMaxKinEnergy / MeV);
  G4double de = (e2 - e1) / G4double(fBinsE);
  G4double dp = (p2 - p1) / G4double(fBinsP);

  G4double x = e1 - de * 0.5;
  G4double e, p, xel, xinel, xtot;
  G4double  kn_xel, kn_xinel, kn_xtot;
  G4int i;

  G4double coeff = 1.0;
  G4double density = 1.0;
  G4double Nb_atom  = 0;
  G4double Mass_atom = 0;
  G4double Nbatom_PerVolume = 0;
  G4double n_BGO=0;


  density = fTargetMaterial->GetDensity();
  coeff   = fTargetMaterial->GetDensity() * cm2 / g;
  Nbatom_PerVolume  = fTargetMaterial->GetTotNbOfAtomsPerVolume();

  const G4ElementVector* elementVector = fTargetMaterial->GetElementVector();
  const G4int* atomsVector = fTargetMaterial->GetAtomsVector(); // 每种元素的原子个数
  const G4double* atomDensities = fTargetMaterial->GetVecNbOfAtomsPerVolume();
  const G4double* fractionMass = fTargetMaterial->GetFractionVector();

  for (size_t i = 0; i < fTargetMaterial->GetNumberOfElements(); ++i) {
      G4double atomicMass = (*elementVector)[i]->GetAtomicMassAmu();  // amu 单位
      G4double numAtoms = (atomsVector) ? atomsVector[i] : -1;  // 获取原子个数
      Mass_atom += numAtoms * atomicMass;
      Nb_atom   += numAtoms;
      G4cout << "Element: " << (*elementVector)[i]->GetName()  << 
                "  Atoms per Volume : " << atomDensities[i] << " 1/mm³" <<
                "  Number of Atoms : " << numAtoms << 
                "  atomicMass : "   << atomicMass << 
                "  fractionMass : " << fractionMass[i] << 
                G4endl;
  }
  n_BGO = (density * cm3/g )/( Mass_atom * mole / Nb_atom ) * (CLHEP::Avogadro * mole );
  G4cout << "Calculated Molecular Mass: " << Mass_atom << " g/mol" << " , Number of Atoms in Molecular " << Nb_atom << G4endl;
  G4cout <<  "Material Density : " << G4BestUnit(density, "Volumic Mass") << G4endl;
  G4cout <<  "Mean Path : " << G4BestUnit(1/coeff,"Length") << G4endl;
  G4cout <<  "Total Number of Atoms per Volume :" << G4BestUnit(1/Nbatom_PerVolume,"Volume") << G4endl;
  G4cout << "n_BGO: " << n_BGO  << " cm-3 " << G4endl;

  G4cout << "### Fill Cross Sections for " << fParticleName << " off " << fMaterialName << G4endl;
  G4cout << "-----------------------------------------------------------------------------------------------" << G4endl;
  G4cout << "    N     E_k(MeV)   Lambda_Inel(cm)   Lambda_El(cm)   Inelastic(b)   Elastic(b)   Total(b)    " << G4endl;
  G4cout << "-----------------------------------------------------------------------------------------------" << G4endl;
  for (i = 0; i < fBinsE; i++) {
    x += de;
    e = std::pow(10., x) * MeV;
    G4cout << std::setw(5) << i << std::setw(12) << e ;
    if (fTargetMaterial) {
      xinel = store->GetInelasticCrossSectionPerVolume(particle, e, fTargetMaterial);
      xel   = store->GetElasticCrossSectionPerVolume(particle, e, fTargetMaterial);
      xtot  =  (xinel + xel);

      kn_xinel = store->GetInelasticCrossSectionPerVolume(particle, e*A, fTargetMaterial);
      kn_xel   = store->GetElasticCrossSectionPerVolume(particle, e*A, fTargetMaterial);
      kn_xtot  =  (kn_xinel + kn_xel);

      fAnalysisManager->FillH1(1, x, 1/xinel/cm );
      fAnalysisManager->FillH1(2, x, 1/xel/cm);
      fAnalysisManager->FillH1(3, x, 1/xtot/cm);
      fAnalysisManager->FillH1(4, x, xinel/n_BGO*cm*1e24 );
      fAnalysisManager->FillH1(5, x, xel/n_BGO*cm*1e24 );
      fAnalysisManager->FillH1(6, x, xtot/n_BGO*cm*1e24 );
      fAnalysisManager->FillH1(7, x, kn_xinel/n_BGO*cm*1e24 );
      fAnalysisManager->FillH1(8, x, kn_xel/n_BGO*cm*1e24 );
      fAnalysisManager->FillH1(9, x, kn_xtot/n_BGO*cm*1e24 );
      G4cout 
      << std::setw(12) << G4BestUnit(1/xinel/cm, "Length") 
      << std::setw(12) << G4BestUnit(1/xel/cm, "Length") 
      << std::setw(12) << G4BestUnit(1/xtot/cm, "Length") 
      << std::setw(12) << G4BestUnit(xinel/n_BGO*cm3, "Surface") 
      << std::setw(12) << G4BestUnit(xel/n_BGO*cm3, "Surface") 
      << std::setw(12) << G4BestUnit(xtot/n_BGO*cm3, "Surface") 
      // << std::setw(12) << xinel/n_BGO  << " barn "
      // << std::setw(12) << xel/n_BGO    << " barn "
      // << std::setw(12) << xtot/n_BGO * 1e24  << " barn "
      << G4endl;
    }
  }
  G4cout <<  "End Energy Loop"  << G4endl;


  if (fVerbose > 0) {
    G4cout << "-------------------------------------------------------------" << G4endl;
  }
  G4cout.precision(prec);
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
  fAnalysisManager->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetVerbose(G4int val)
{
  fVerbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
