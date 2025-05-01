#include "multiColStudy.h"
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdOut.h>
#include <calobase/RawTowerv1.h>
#include <phool/recoConsts.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <ffaobjects/EventHeader.h>
#include <jetbase/Jet.h>
#include <jetbase/JetContainerv1.h>
using namespace std;
static const float radius_EM = 93.5;
static const float minz_EM = -130.23;
static const float maxz_EM = 130.23;

static const float radius_IH = 127.503;
static const float minz_IH = -170.299;
static const float maxz_IH = 170.299;

static const float radius_OH = 225.87;
static const float minz_OH = -301.683;
static const float maxz_OH = 301.683;
float get_emcal_mineta_zcorrected(float zvertex) {
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

//____________________________________________________________________________..
multiColStudy::multiColStudy(const std::string &filename, const std::string &name, const int debug, const int issim):
  SubsysReco(name)
{
  _issim = issim;
  _name = name;
  _debug = debug;
  _filename = filename;
  _f = new TFile(_filename.c_str(), "RECREATE");
  _tree = new TTree("tree","a persevering date tree");
}

//____________________________________________________________________________..
multiColStudy::~multiColStudy()
{

}

//____________________________________________________________________________..
int multiColStudy::Init(PHCompositeNode *topNode)
{
  
  if(_debug > 1) cout << "Begin init: " << endl;
  _tree->Branch("njet",&_njet,"njet/I");
  _tree->Branch("jet_e",_jet_e,"jet_e[njet]/D");
  _tree->Branch("jet_et",_jet_et,"jet_et[njet]/D");
  _tree->Branch("jet_eta",_jet_eta,"jet_eta[njet]/D");
  _tree->Branch("jet_phi",_jet_phi,"jet_phi[njet]/D");
  if(!_issim) _tree->Branch("trigvec",&_trigvec,"trigvec/l");
  if(_issim)
    {
      _tree->Branch("tnjet",&_tnjet,"tnjet/I");
      _tree->Branch("tjet_e",_tjet_e,"tjet_e[tnjet]/D");
      _tree->Branch("tjet_et",_tjet_et,"tjet_et[tnjet]/D");
      _tree->Branch("tjet_eta",_tjet_eta,"tjet_eta[tnjet]/D");
      _tree->Branch("tjet_phi",_tjet_phi,"tjet_phi[tnjet]/D");
      _tree->Branch("tdphilead",&_tdphilead,"tdphilead/D");
      _tree->Branch("tisdijet",&_tisdijet,"tisdijet/I");
    }
  _tree->Branch("nzvtx",&_nzvtx,"nzvtx/I");
  if(_issim) _tree->Branch("tnzvtx",&_tnzvtx,"tnzvtx/I");
  _tree->Branch("rzvtx",_rzvtx,"rzvtx[nzvtx]/D");
  if(_issim) _tree->Branch("tzvtx",_tzvtx,"tzvtx[tnzvtx]/D");
  _tree->Branch("hitsgrone",&_hitsgrone,"hitsgrone/I");
  _tree->Branch("towgrone",_towgrone,"towgrone[hitsgrone][5]/D");
  _tree->Branch("mbdq",_mbdq,"mbdq[2][64]/D");
  _tree->Branch("frcem",_frcem,"frcem[njet]/D");
  _tree->Branch("frcoh",_frcoh,"frcoh[njet]/D");
  _tree->Branch("dphilead",&_dphilead,"dphilead/D");
  _tree->Branch("isdijet",&_isdijet,"isdijet/I");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int multiColStudy::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..

void print_debug(float jet_eta, float jet_phi, float tower_eta, float tower_phi, float dphi, float deta)
{
  cout << "printing debug info for dphi deta:" << endl;
  cout << "jet eta/phi: " << jet_eta << " " << jet_phi << endl;
  cout << "tower eta/phi: " << tower_eta << " " << tower_phi << endl;
  cout << "deta dphi:" << deta << " " << dphi << endl;
}

int multiColStudy::process_event(PHCompositeNode *topNode)
{

  

  MbdVertexMap* mbdvtxmap = findNode::getClass<MbdVertexMapv1>(topNode, "MbdVertexMap");

  float zvtx = NAN;

  if(!_issim)
    {
      Gl1Packetv2* gl1 = findNode::getClass<Gl1Packetv2>(topNode, "GL1Packet");
      if(!gl1)
	{
	  cout << "No trigger info!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
      _trigvec = gl1->getScaledVector();
    }
  _nzvtx = 0;
  if(mbdvtxmap)
    {
      for(auto iter = mbdvtxmap->begin(); iter != mbdvtxmap->end(); ++iter)
        {
          MbdVertex* mbdvtx = iter->second;
          if(mbdvtx)
	    {
	      _rzvtx[_nzvtx] = mbdvtx->get_z();
	      ++_nzvtx;
	      if(_nzvtx > _maxzvtx) break;
	    }
        }
    }
  zvtx = _rzvtx[0];  
  if(std::isnan(zvtx))
    {
      if(_debug > 1) cout << "no good zvtx!" << endl;
      goto badz;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }


  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainerv4>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainerv4>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainerv4>(topNode, "TOWERINFO_CALIB_HCALOUT");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r04");//"AntiKt_unsubtracted_r04");
  if(!jets) jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_unsubtracted_r04");
  JetContainer* truthjets = findNode::getClass<JetContainerv1>(topNode,"AntiKt_Truth_r04");
  MbdPmtContainer * mbdtow = findNode::getClass<MbdPmtContainer>(topNode,"MbdPmtContainer");
  if(!mbdtow) mbdtow = findNode::getClass<MbdPmtContainerV1>(topNode,"MbdPmtContainer");
  if(!mbdtow) mbdtow = findNode::getClass<MbdPmtSimContainerV1>(topNode,"MbdPmtContainer");
  if(_debug > 2) cout << towers[0] << " " << towers[1] << " " << towers[2] << endl;
  
  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT");
  
  //TowerInfoContainer* emcrt = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  int sectormb=128;
  if(mbdtow)
    {
      for(int i=0; i<sectormb; ++i)
	{
	  MbdPmtHit* mbdhit = mbdtow->get_pmt(i);
	  _mbdq[(i<64?0:1)][(i<64?i:i-64)] = mbdhit->get_q();
	}
    }
  else
    {
      cout << "No MBD info!" << endl;
    }
  _hitsgrone = 0;
  int nchan = 1536;
  for(int h=0; h<_ncalotype; ++h)
    {
      if(towers[h])
	{
	  if(_debug > 1) cout << "got towers " << h << endl;
	  for(int i=0; i<nchan; ++i)
	    {
	      TowerInfo* tower = towers[h]->get_tower_at_channel(i);
	      if(!tower->get_isGood()) continue;
	      int key = towers[h]->encode_key(i);
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(h<2?RawTowerDefs::CalorimeterId::HCALIN:RawTowerDefs::CalorimeterId::HCALOUT, towers[h]->getTowerEtaBin(key), towers[h]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[h<2?1:2]->get_tower_geometry(geomkey);
	      float radius = h<1?93.5:tower_geom->get_center_radius();
	      float ihEta = tower_geom->get_eta();
	      float emZ = radius/(tan(2*atan(exp(-ihEta))));
	      float newz = h<1?(emZ - _rzvtx[0]):(tower_geom->get_center_z() - _rzvtx[0]);
	      float newTheta = atan2(radius,newz);
	      float towerEta = -log(tan(0.5*newTheta));
	      float towerE = tower->get_energy();///cosh(towerEta);
	      if(towerE > 1)
		{
		  if(_debug > 2) cout << "hit greater than 1" << endl;
		  _towgrone[_hitsgrone][0] = towerE;
		  _towgrone[_hitsgrone][1] = towerEta;
		  _towgrone[_hitsgrone][2] = tower_geom->get_phi();
		  _towgrone[_hitsgrone][3] = tower_geom->get_eta();
		  _towgrone[_hitsgrone][4] = h;
		  ++_hitsgrone;
		} 
	    }
	}
      else if(_debug > 1)
	{
	  cout << "no towers " << h << endl;
	}
    }
  
  float maxJetE = 0;
  float maxJetPhi = 0;
  float subJetE = 0;
  float subJetPhi = 0;
  _njet = 0;

  if(jets)
    {
      int tocheck = jets->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *jet = jets->get_jet(i);
          if(jet)
            {
	      if(_debug > 5) cout << "getting jet E/eta" << endl;
	      float testJetE = jet->get_e();
	      float testJetPhi = jet->get_phi();
	      if(_debug > 5) cout << "jet E/eta: " << testJetE  << " " << jet->get_eta() << endl;
	      if(testJetE < 7) continue;
	      if(_debug > 3) cout << "got a candidate jet" << endl;
	      _jet_eta[_njet] = jet->get_eta();
	      if(check_bad_jet_eta(_jet_eta[_njet],_rzvtx[0],0.4)) continue;
	      _frcem[_njet] = 0;
	      _frcoh[_njet] = 0;
	      _jet_e[_njet] = testJetE;
	      _jet_et[_njet] = testJetE/cosh(_jet_eta[_njet]);
	      _jet_phi[_njet] = testJetPhi;
	      int ncomp = 0;
	      if(_jet_et[_njet] > subJetE && _jet_et[_njet] < maxJetE)
		{
		  subJetE = _jet_et[_njet];
		  subJetPhi = testJetPhi;
		}
	      if(_debug > 2) cout << "found a good jet!" << endl;
	      if(_jet_et[_njet] > maxJetE)
		{
		  subJetE = maxJetE;
		  subJetPhi = maxJetPhi;
		  maxJetE = _jet_et[_njet];
		  maxJetPhi = _jet_phi[_njet];
		}
	
	      if(_debug > 3) cout << "getting comp vec" << endl;
	      
	      for(auto comp: jet->get_comp_vec())
		{
		  ++ncomp;
		  unsigned int channel = comp.second;
		  TowerInfo* tower;
		  //cout << "type: " << comp.first << endl;
		  int towerType = -1;
		  if(comp.first == 5 || comp.first == 26)
		    {
		      towerType = 1;
		    }
		  else if(comp.first == 7 || comp.first == 27)
		    {
		      towerType = 2;
		    }
		  else if(comp.first == 13 || comp.first == 25 || comp.first == 28)
		    {
		      towerType = 0;
		    }
		  else
		    {
		      cout << "BAD TOWERTYPE!!" << endl;
		      continue;
		    }
		  tower = towers[towerType]->get_tower_at_channel(channel);
		  float towerE = tower->get_energy();
		  //int key = towers[towerType]->encode_key(channel);
		  //const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[towerType]->getTowerEtaBin(key), towers[towerType]->getTowerPhiBin(key));
		  if(_debug > 6) cout << "encoding tower geom" << endl;
		  //RawTowerGeom *tower_geom = geom[towerType==2?2:1]->get_tower_geometry(geomkey); //encode tower geometry
		  //float radius = towerType==0?93.5:tower_geom->get_center_radius();
		  //float emZ = radius/(tan(2*atan(exp(-tower_geom->get_eta()))));
		  //float newz = towerType==0?(emZ - zvtx):(tower_geom->get_center_z() - zvtx);
		  //float newTheta = atan2(radius,newz);
		  //float towerEta = -log(tan(0.5*newTheta));
		  //float towerPhi = tower_geom->get_phi();
		  if(towerType==0) _frcem[_njet] += towerE;
		  else if(towerType==2) _frcoh[_njet] += towerE;
		}
	      _frcoh[_njet] /= _jet_e[_njet];
	      _frcem[_njet] /= _jet_e[_njet];
	      ++_njet;
	      if(_njet > 9) break;
	    }
	  else
	    {
	      continue;
	    }
	}
      /*
      if(!_njet)
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      */
      _dphilead = abs(maxJetPhi-subJetPhi);
      if(_dphilead > M_PI) _dphilead = 2*M_PI - _dphilead;
      if(subJetE > 7) _isdijet = 1;
      else _isdijet = 0;
    }
  else
    {
      if(_debug > 0) cout << "no jets" << endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }

  _tnjet = 0;
  if(truthjets)
    {
      for(int i=0; i<truthjets->size(); ++i)
	{
	  Jet* jet = truthjets->get_jet(i);
	  if(_debug > 5) cout << "getting jet E/eta" << endl;
	  float testJetE = jet->get_e();
	  float testJetPhi = jet->get_phi();
	  if(_debug > 5) cout << "jet E/eta: " << testJetE  << " " << jet->get_eta() << endl;
	  if(testJetE < 4) continue;
	  if(_debug > 3) cout << "got a candidate jet" << endl;
	  _tjet_eta[_tnjet] = jet->get_eta();
	  if(check_bad_jet_eta(_tjet_eta[_tnjet],_rzvtx[0],0.4)) continue;
	  _tjet_e[_tnjet] = testJetE;
	  _tjet_et[_tnjet] = testJetE/cosh(_tjet_eta[_tnjet]);
	  _tjet_phi[_tnjet] = testJetPhi;
	  if(_tjet_et[_tnjet] > subJetE && _tjet_et[_tnjet] < maxJetE)
	    {
	      subJetE = _tjet_et[_tnjet];
	      subJetPhi = testJetPhi;
	    }
	  if(_debug > 2) cout << "found a good jet!" << endl;
	  if(_tjet_et[_tnjet] > maxJetE)
	    {
	      if(maxJetE > subJetE)
		{
		  subJetE = maxJetE;
		  subJetPhi = maxJetPhi;
		}
	      maxJetE = _tjet_et[_tnjet];
	    }
	}
      _tdphilead = abs(maxJetPhi-subJetPhi);
      if(_tdphilead > M_PI) _tdphilead = 2*M_PI - _tdphilead;
      if(subJetE > 4) _tisdijet = 1;
      else _tisdijet = 0;      
    }
  //if(maxJetE > 4)
  //{
  if(_debug > 0) cout << "filling jet tree" << endl;
  
 badz: _tree->Fill();
  //}
  
  if(_debug > 3) cout << "end event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int multiColStudy::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "multiColStudy::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int multiColStudy::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "multiColStudy::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int multiColStudy::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "multiColStudy::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  if(_debug > 1) cout << "ending run" << endl;
  _f->cd();
  _tree->Write();
  _f->Write();
  _f->Close();

  //delete jet_tree;
  //delete _f;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int multiColStudy::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "multiColStudy::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void multiColStudy::Print(const std::string &what) const
{
  std::cout << "multiColStudy::Print(const std::string &what) const Printing info for " << what << std::endl;
}
