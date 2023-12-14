#include <string>

#include <TFile.h>
#include <TChain.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TGeoGlobalMagField.h>
#include <vector>

#include <Generators/GeneratorFactory.h>
#include "FairPrimaryGenerator.h"
#include "FairGenerator.h"
#include "FairBoxGenerator.h"
#include <FairLogger.h>
#include <SimConfig/SimConfig.h>
#include <Generators/GeneratorFromFile.h>

#include "SimulationDataFormat/MCTrack.h"
#include "MathUtils/Cartesian.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/TrackParametrizationWithError.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <Generators/GeneratorFactory.h>
#include "FairPrimaryGenerator.h"
#include "FairGenerator.h"
#include "FairBoxGenerator.h"
#include <FairLogger.h>
#include <SimConfig/SimConfig.h>
#include <Generators/GeneratorFromFile.h>

#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "Framework/DataTypes.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "ITSBase/GeometryTGeo.h"
#include "TStopwatch.h" //for bechmarking

/// get mass from TParticlePDG
double getMass(int input_pdg){
  double mass = 0;
  if(TDatabasePDG::Instance()){
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(input_pdg);
    if(particle)  mass = particle->Mass();
    else      std::cout << "===> particle mass equal to 0" << std::endl;
  }
  return mass;
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
template <typename McParticleType>
o2::track::TrackParCov convertMCParticleToO2Track(McParticleType& particle)
{
  // FIXME: this is a fundamentally important piece of code.
  // It could be placed in a utility file instead of here.
  auto pdgInfo = TDatabasePDG::Instance()->GetParticle(particle.GetPdgCode());
  int charge = 0;
  if (pdgInfo != nullptr) {
    charge = pdgInfo->Charge() / 3;
  }
  std::array<float, 5> params;
  std::array<float, 15> covm = {0.};
  float s, c, x;
  o2::math_utils::sincos(particle.GetPhi(), s, c);
  o2::math_utils::rotateZInv(particle.Vx(), particle.Vy(), x, params[0], s, c);
  params[1] = particle.Vz();
  params[2] = 0.; // since alpha = phi
  auto theta = 2. * std::atan(std::exp(-particle.GetEta()));
  params[3] = 1. / std::tan(theta);
  params[4] = charge / particle.GetPt();

  // Return TrackParCov
  return o2::track::TrackParCov(x, particle.GetPhi(), params, covm);
}

/// function to calculate track length of this track up to a certain radius
/// \param track the input track
/// \param radius the radius of the layer you're calculating the length to
/// \param magneticField the magnetic field to use when propagating
float trackLength(o2::track::TrackParCov track, float radius, float magneticField)
{
  // don't make use of the track parametrization
  float length = -100;

  o2::math_utils::CircleXYf_t trcCircle;
  float sna, csa;
  track.getCircleParams(magneticField, trcCircle, sna, csa);

  // distance between circle centers (one circle is at origin -> easy)
  float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

  // condition of circles touching - if not satisfied returned length will be -100
  if (centerDistance < trcCircle.rC + radius && centerDistance > fabs(trcCircle.rC - radius)) {
    length = 0.0f;

    // base radical direction
    float ux = trcCircle.xC / centerDistance;
    float uy = trcCircle.yC / centerDistance;
    // calculate perpendicular vector (normalized) for +/- displacement
    float vx = -uy;
    float vy = +ux;
    // calculate coordinate for radical line
    float radical = (centerDistance * centerDistance - trcCircle.rC * trcCircle.rC + radius * radius) / (2.0f * centerDistance);
    // calculate absolute displacement from center-to-center axis
    float displace = (0.5f / centerDistance) * TMath::Sqrt(
                                                 (-centerDistance + trcCircle.rC - radius) *
                                                 (-centerDistance - trcCircle.rC + radius) *
                                                 (-centerDistance + trcCircle.rC + radius) *
                                                 (centerDistance + trcCircle.rC + radius));

    // possible intercept points of track and TOF layer in 2D plane
    float point1[2] = {radical * ux + displace * vx, radical * uy + displace * vy};
    float point2[2] = {radical * ux - displace * vx, radical * uy - displace * vy};

    // decide on correct intercept point
    std::array<float, 3> mom;
    track.getPxPyPzGlo(mom);
    float scalarProduct1 = point1[0] * mom[0] + point1[1] * mom[1];
    float scalarProduct2 = point2[0] * mom[0] + point2[1] * mom[1];

    // get start point
    std::array<float, 3> startPoint;
    track.getXYZGlo(startPoint);

    if(std::hypot(startPoint[0], startPoint[1])>radius) return -100.0f;
    
    float cosAngle = -1000, modulus = -1000;

    if (scalarProduct1 > scalarProduct2) {
      modulus = std::hypot(point1[0] - trcCircle.xC, point1[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
      cosAngle = (point1[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point1[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    } else {
      modulus = std::hypot(point2[0] - trcCircle.xC, point2[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
      cosAngle = (point2[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point2[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    }
    cosAngle /= modulus;
    length = trcCircle.rC * TMath::ACos(cosAngle);
    length *= sqrt(1.0f + track.getTgl() * track.getTgl());
  }
  return length;
}



void runacceptancetest(TString lDataPath = "./", TString lOutFile =  "output.root")
{
  
  TStopwatch lWatchTracking, lWatchAnalysis, lWatchSmearing;
  
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  std::cout << "\e[1;31m      --- ALICE 3 Acceptance tester ---        \e[0;00m" << std::endl;
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  
  // master control of layer setup
  //Double_t layers[] = {0.5, 1.2, 2.5, 3.75, 7, 12, 20, 30, 45, 60, 80};
  
  // 10 layers (hypothetical descoped OT)
  //Double_t layers[] = {0.5f,1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 50.0f, 80.0f};
  
  // 13 layers (hypothetical descoped OT)
  // Double_t layers[] = {0.5f,1.2f, 2.5f, 3.75f, 5.0f, 7.0f, 9.0f, 15.0f, 20.0f, 30.0f, 45.0f, 65.0f, 80.0f};
  
  // 12 layers
  Double_t layers[] = {0.5f,1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 38.0f, 51.0f, 68.0f};
  
  //Double_t layers[] = {0.5, 1.2, 2.5, 3.75, 5.0, 6.5, 8.0, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 22.5, 25.0, 30.0, 35.0, 45.0, 60.0, 80.0};
  
  Int_t nlayers = sizeof(layers)/sizeof(Double_t);
  int minHits = 6;
  cout<<"N layers = "<<nlayers<<endl;
  
  cout<<"Main configuration: "<<endl;
  cout<<"Data path for analysis.......: "<<lDataPath.Data()<<endl;
  cout<<"Data path for LUTs...........: "<<lOutFile.Data()<<endl;
  
  //Open GRP
  const auto grp = o2::parameters::GRPObject::loadFrom(Form("%so2sim_grp.root",lDataPath.Data()));
  if (!grp) {
    LOG(FATAL) << "Cannot run w/o GRP object";
  }
  
  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  if (!field) {
    LOG(FATAL) << "Failed to load mag field";
  }
  double origD[3] = {0., 0., 0.};
  
  //Operational parameters
  const Double_t lMagneticField = field->GetBz(0,0,0);
  cout<<"Magnetic field auto-detected to be "<<lMagneticField<<endl;
  
  TChain mcTree("o2sim");
  mcTree.AddFile(Form("%so2sim_Kine.root", lDataPath.Data()));
  
  o2::dataformats::MCEventHeader* mcHead = nullptr;
  mcTree.SetBranchStatus("*", 0); //disable all branches
  mcTree.SetBranchStatus("MCTrack*", 1);
  mcTree.SetBranchStatus("MCEventHeader.*", 1);
  mcTree.SetBranchAddress("MCEventHeader.", &mcHead);
  
  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree.SetBranchAddress("MCTrack", &mcArr);

  Double_t pTbins[] = {0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f};
  Int_t nPtBins = sizeof(pTbins) / sizeof(Double_t) - 1;
  
  TFile *fOutput = new TFile(lOutFile.Data(), "RECREATE");
  //TH1D *hGeneratedXiCC = new TH1D("hGeneratedXiCC", "", 200, 0, 20);
  //TH1D *hReconstructedXiCC = new TH1D("hReconstructedXiCC", "", 200, 0, 20);
  TH1D *hGeneratedXiCC = new TH1D("hGeneratedXiCC", "", nPtBins, pTbins);
  TH1D *hReconstructedXiCC = new TH1D("hReconstructedXiCC", "", nPtBins, pTbins);
  TH1D *hLayerArrangement = new TH1D("hLayerArrangement", "", nlayers-1, layers);
  
  //store min hits in bin 1 of layer arrangement
  hLayerArrangement->SetBinContent(1, minHits);
  
  cout<<"Number of events in total: "<<mcTree.GetEntries()<<endl;
  for (uint32_t iEvent{0}; iEvent < mcTree.GetEntries(); ++iEvent) {
    
    //std::cout << "*************** Event " << iEvent << " ("<<mcArr->size()<<" particles) ****" << std::endl;
    mcTree.GetEvent(iEvent);
    Double_t pTxicc = 0.0f;
    for (Long_t iii=0; iii< mcArr->size(); iii++ ){
      auto part = mcArr->at(iii);
      // identify if XiCC as desired
      if(part.GetPdgCode()==4422 && TMath::Abs(part.GetEta())<1.5){
        pTxicc = std::hypot(part.Px(),part.Py());
        
        // inner loop: find daughters and check if they travel through certain layers
        o2::MCTrack firstXiCCDau = mcArr->at(part.getFirstDaughterTrackId());
        o2::MCTrack lastXiCCDau = mcArr->at(part.getLastDaughterTrackId());
        
        if(firstXiCCDau.GetPdgCode() != 4232 || lastXiCCDau.GetPdgCode() != 211){
          cout<< "PANIC ! Something very bad happened!"<<endl;
          return;
        }
        
        o2::MCTrack firstXiCDau = mcArr->at(firstXiCCDau.getFirstDaughterTrackId());
        o2::MCTrack secondXiCDau = mcArr->at(firstXiCCDau.getFirstDaughterTrackId()+1);
        o2::MCTrack lastXiCDau = mcArr->at(firstXiCCDau.getLastDaughterTrackId());
        
        o2::MCTrack firstXiDau = mcArr->at(firstXiCDau.getFirstDaughterTrackId());
        o2::MCTrack lastXiDau = mcArr->at(firstXiCDau.getLastDaughterTrackId());
        
        // protect against other decay channels
        if(firstXiDau.getLastDaughterTrackId()!=firstXiDau.getFirstDaughterTrackId()+1){
          continue;
        }
        
        o2::MCTrack firstLambdaDau = mcArr->at(firstXiDau.getFirstDaughterTrackId());
        o2::MCTrack lastLambdaDau = mcArr->at(firstXiDau.getLastDaughterTrackId());
        
        // protect against other decay channels (2)
        if(firstLambdaDau.GetPdgCode()!=2212 || lastLambdaDau.GetPdgCode()!=-211) continue;
        
        // at this stage, I already know what's what and can decide if stuff reached certain layers or not
        o2::track::TrackParCov daughterTrack[6];
        daughterTrack[0] = convertMCParticleToO2Track(lastXiCCDau);
        daughterTrack[1] = convertMCParticleToO2Track(secondXiCDau);
        daughterTrack[2] = convertMCParticleToO2Track(lastXiCDau);
        daughterTrack[3] = convertMCParticleToO2Track(lastXiDau);
        daughterTrack[4] = convertMCParticleToO2Track(lastLambdaDau);
        daughterTrack[5] = convertMCParticleToO2Track(firstLambdaDau);
        
        // only desired decay
        hGeneratedXiCC->Fill(pTxicc);
        
        bool masterDetect = true;
        int daughterHits[6]={0,0,0,0,0,0};
        for(Int_t id=0;id<6;id++){
          for(Int_t il=0;il<nlayers;il++){
            daughterHits[id] += static_cast<int>(trackLength(daughterTrack[id], layers[il], lMagneticField) > 0.0f);
          }
          if(daughterHits[id]<minHits) masterDetect = false;
        }
        if(masterDetect) hReconstructedXiCC->Fill(pTxicc);
      }
    }
  }
  
  fOutput->Write();
}

