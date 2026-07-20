//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "DDG4/Factories.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "G4EmProcessSubType.hh"
#include "G4OpticalPhoton.hh"
#include "G4Poisson.hh"
#include "G4ProcessType.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include "detectorSegmentations/SCEPCal_MainSegmentation_k4geo.h"

namespace SCEPCal {
// PDE SIPM 3x3
inline const std::vector<std::pair<double, double>>& SiPM_PDE_table() {
  static const std::vector<std::pair<double, double>> tbl = {
    {281.741, 0.005957}, {283.227, 0.010426}, {283.227, 0.017128}, {283.227, 0.025319},
    {287.686, 0.031277}, {290.658, 0.038723}, {293.631, 0.044681}, {296.603, 0.050638},
    {301.062, 0.055106}, {302.548, 0.061809}, {308.493, 0.067766}, {314.437, 0.073723},
    {320.382, 0.080426}, {324.841, 0.087872}, {330.786, 0.093830}, {333.758, 0.099787},
    {338.217, 0.105745}, {344.161, 0.110957}, {353.079, 0.116170}, {366.454, 0.119149},
    {375.372, 0.123617}, {382.803, 0.130319}, {387.261, 0.136277}, {393.206, 0.142979},
    {402.123, 0.149681}, {406.582, 0.154894}, {412.527, 0.160851}, {419.958, 0.165319},
    {425.902, 0.169787}, {431.847, 0.175000}, {442.251, 0.179468}, {454.140, 0.182447},
    {466.030, 0.183936}, {479.406, 0.182447}, {497.240, 0.180957}, {510.616, 0.179468},
    {525.478, 0.176489}, {535.881, 0.171277}, {544.798, 0.166809}, {555.202, 0.161596},
    {564.119, 0.155638}, {576.008, 0.150426}, {583.439, 0.145957}, {593.843, 0.140000},
    {602.760, 0.134787}, {613.163, 0.129574}, {623.567, 0.125106}, {635.456, 0.119894},
    {645.860, 0.113936}, {657.749, 0.107234}, {669.639, 0.101277}, {680.042, 0.096064},
    {693.418, 0.089362}, {705.308, 0.083404}, {717.197, 0.077447}, {729.087, 0.070745},
    {742.463, 0.064787}, {757.325, 0.059574}, {767.728, 0.054362}, {779.618, 0.049894},
    {794.480, 0.043191}, {809.342, 0.038723}, {827.176, 0.033511}, {842.038, 0.029787},
    {855.414, 0.025319}, {871.762, 0.022340}, {886.624, 0.019362}, {895.541, 0.017872}
  };
  return tbl;
}

// Interpolazione tra i punti della pde
inline double SiPM_PDE(double lambda_nm) {
  const auto& tbl = SiPM_PDE_table();
  if (lambda_nm <= tbl.front().first) return tbl.front().second;
  if (lambda_nm >= tbl.back().first)  return tbl.back().second;

  auto it = std::upper_bound(tbl.begin(), tbl.end(), lambda_nm,
      [](double val, const std::pair<double, double>& p) { return val < p.first; });
  auto it0 = std::prev(it);
  const double x0 = it0->first, y0 = it0->second;
  const double x1 = it->first,  y1 = it->second;
  const double t  = (lambda_nm - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}




class SCEPCal_MainSDAction {
public:
  std::size_t m_collectionID_scint{0};
  std::size_t m_collectionID_ceren{0};

  
  //double m_scintPhotoEleGeV{21372.59};   // PWO after pde
  double m_scintPhotoEleGeV{246606.84};   // BSO after pde
  double m_geomEff{0.09};  // (3mm x 3mm) / (10mm x 10mm) SiPM su cristallo
  //double m_geomEff{1};  // Nominal

  //double m_scintPhotoEleGeV{21372.59};
  // double m_scintPhotoEleGeV{1500000.0};   // BSO
  //double m_scintPhotoEleGeV{1500000.0};   // default value
  //double m_cherenBinomialCut{0.00147}; // default value


  

  // This function takes in input the ionizing energy deposit (in MeV) from a step in the main scepcal section
  // and returns the corresponding number of sipm fired cells (or photo-electrons).
  // It reproduces a scintillating light yield of 2000 p.e./GeV which is expected from test-beam data.
  // The smearing is done by a single poissoninan sampling. De facto one must have two samplings:
  // a poissonian sampling for light emission fluctuations and a Binomial sampling for light detection
  // fluctuations. Thanks to the Poissonian thinning theorem this is equivalent to a single poissonian
  // sampling reproducing on average the expected light yield.
  int SmearSsignal(G4double edep) const {
    return G4Poisson((edep / CLHEP::GeV) * m_scintPhotoEleGeV * m_geomEff); // MeV->GeV
  }

  // This function is a Bernoully trial to decide wether a Cerenkov photon is detected or not
  // in the scepcal. Cerenkov photons are created by Geant4 inside crystals, already taking
  // into account their poissonian emission fluctuations, therefore we only need to add a
  // Binomial sampling that reproduces the expected light yield of 100 p.e./GeV
  //bool SmearCsignal() const { return G4UniformRand() < m_cherenBinomialCut; } //originale
  bool SmearCsignal(double lambda_nm) const { return G4UniformRand() < SiPM_PDE(lambda_nm) * m_geomEff; } 

};

} // namespace SCEPCal

namespace dd4hep {
namespace sim {
  using namespace SCEPCal;

  template <>
  void Geant4SensitiveAction<SCEPCal_MainSDAction>::defineCollections() {
    m_collectionID = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainEdep");
    m_userData.m_collectionID_scint = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainScounts");
    m_userData.m_collectionID_ceren = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainCcounts");
  }

  template <>
  Geant4SensitiveAction<SCEPCal_MainSDAction>::Geant4SensitiveAction(Geant4Context* ctxt, const std::string& nam,
                                                                     DetElement det, Detector& lcdd_ref)
      : Geant4Sensitive(ctxt, nam, det, lcdd_ref), m_collectionID(0) {
    initialize();
    defineCollections();
    InstanceCount::increment(this);
    //declareProperty("scintPhotoElePerGeV", m_userData.m_scintPhotoEleGeV = 21372.59); // al posto di 130000
    //declareProperty("scintPhotoElePerGeV", m_userData.m_scintPhotoEleGeV = 1500000); 
    declareProperty("scintPhotoElePerGeV", m_userData.m_scintPhotoEleGeV = 246606.84); // al posto di 1500000
    //declareProperty("cherenBinomialCut", m_userData.m_cherenBinomialCut = 0.00147);
    //declareProperty("scintPhotoElePerGeV", m_userData.m_scintPhotoEleGeV = 130000);
    declareProperty("geomEfficiency", m_userData.m_geomEff = 0.09);
  }

  template <>
  bool Geant4SensitiveAction<SCEPCal_MainSDAction>::process(const G4Step* step, G4TouchableHistory* /*hist*/) {
    G4StepPoint* thePrePoint = step->GetPreStepPoint();
    G4TouchableHandle thePreStepTouchable = thePrePoint->GetTouchableHandle();

    auto vID = thePreStepTouchable->GetCopyNumber(0);
    G4Track* track = step->GetTrack();
    // get phi from grandparent
    auto cellID2 = thePreStepTouchable->GetCopyNumber(2);

    dd4hep::Segmentation* _geoSeg = &m_segmentation;
    auto segmentation = dynamic_cast<dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo*>(_geoSeg->segmentation());
    G4double edep = step->GetTotalEnergyDeposit();

    auto newOrExistingHitIn = [&](std::size_t id, bool isCherenkov = false) -> Geant4Calorimeter::Hit* {
      Geant4HitCollection* coll = collection(id);

      auto cellID_phi = segmentation->setPhi(vID, cellID2);
      auto cIDwithChannel = segmentation->setCellID(cellID_phi, isCherenkov);

      auto* hit = coll->findByKey<Geant4Calorimeter::Hit>(cIDwithChannel);

      if (!hit) {
        DDSegmentation::Vector3D pos =
            segmentation->position(cellID_phi); // always use scintillation channel for position
        Position global(pos.x(), pos.y(), pos.z());
        hit = new Geant4Calorimeter::Hit(global / dd4hep::mm);

        hit->cellID = cIDwithChannel;
        coll->add(cIDwithChannel, hit);
      }

      return hit;
    };

    // Cerenkov hits
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      auto procSubType = track->GetCreatorProcess()->GetProcessSubType();
      bool isCerenkov = (procSubType == G4EmProcessSubType::fCerenkov);
      if (!isCerenkov) {
        track->SetTrackStatus(fStopAndKill);
        return true;
      } else if (track->GetCurrentStepNumber() == 1) {
<<<<<<< HEAD
        // create a hit only if the photon is detected
        if (m_userData.SmearCsignal()) {
          auto* hitC = newOrExistingHitIn(m_userData.m_collectionID_ceren, true);
          // The EDM4HEP converter divides this entry by CLHEP::GeV to save energies in GeV unit.
          // Here, we are saving C counts so we multiply by CLHEP::GeV to have the correct number in the output file.
=======
        auto* hitC = newOrExistingHitIn(m_userData.m_collectionID_ceren, true);

        // lungh. d'onda fotone Cherenkov dall'energia cinetica
        const double photonEnergy_eV = track->GetKineticEnergy() / CLHEP::eV;
        const double lambda_nm       = 1239.84198 / photonEnergy_eV; // hc [eV*nm]



        // The EDM4HEP converter divides this entry by CLHEP::GeV to save energies in GeV unit.
        // Here, we are saving C counts so we multiply by CLHEP::GeV to have the correct number in the output file.
        // prossima riga da commentare se vuoi solo il numero di fotoni Cherenkov di Geant4
        if (m_userData.SmearCsignal(lambda_nm))
>>>>>>> 309c3bac (Change crystal and PDE response function)
          hitC->energyDeposit += 1 * CLHEP::GeV;
        } // if SmearCsignal

        // counted -> kill || not detected -> kill
        track->SetTrackStatus(fStopAndKill);
      } else // keep it until it passes the reflection
        return true;
    } // if optical photon

    if (edep <= 0.)
      return true; // do not save hits with no energy deposit

    // Scintillation signal for this step.  NOTE: this early return also gates the edep collection
    // below, so SCEPCal_MainEdep is deliberately NOT a complete record of the ionizing energy
    // deposited in the crystals.  It stores only the deposits that produced a scintillation signal,
    // so that every edep hit -- and the MC-truth contribution attached to it -- has a matching S hit
    // in the same cell.
    auto Scount = m_userData.SmearSsignal(edep);
    if (Scount <= 0)
      return true;

    // edep hits (truth): only for steps with a scintillation signal -- see the note above
    auto* hitedep = newOrExistingHitIn(m_collectionID, false);
    HitContribution contrib = Geant4Calorimeter::Hit::extractContribution(step);
    hitedep->energyDeposit += contrib.deposit;
    hitedep->truth.emplace_back(contrib);
    Geant4StepHandler h(step);
    mark(h.track);

    // scintillation hits
    auto* hitS = newOrExistingHitIn(m_userData.m_collectionID_scint, false);
    // The EDM4HEP converter divides this entry by CLHEP::GeV to save energies in GeV unit.
    // Here, we are saving S counts so we multiply by CLHEP::GeV to have the correct number in the output file.
    hitS->energyDeposit += Scount * CLHEP::GeV;

    return true;
  } // process
} // namespace sim
} // namespace dd4hep

namespace dd4hep {
namespace sim {
  typedef Geant4SensitiveAction<SCEPCal_MainSDAction> SCEPCal_MainSDAction;
}
} // namespace dd4hep

DECLARE_GEANT4SENSITIVE(SCEPCal_MainSDAction)
