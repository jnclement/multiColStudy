#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <caloreco/CaloTowerStatus.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbase/FastJetAlgo.h>
//#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/BeamBackgroundFilterAndQA.h>
#include <fstream>
#include <phool/recoConsts.h>
#include <TSystem.h>
#include <caloreco/CaloTowerCalib.h>
#include <frog/FROG.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <globalvertex/GlobalVertexReco.h>
#include <GlobalVertex.h>
#include <multicolstudy/multiColStudy.h>
#include <TruthJetInput.h>//#include <G4Setup_sPHENIX.C>
//#include <trigzvtxchecker/Trigzvtxchecker.h>
#include <MbdDigitization.h>
#include <MbdReco.h>
#include <Calo_Calib.C>
#include <calotrigger/CaloTriggerEmulator.h>
using namespace std;

R__LOAD_LIBRARY(libg4centrality.so)
//R__LOAD_LIBRARY(libFROG.so)
//R__LOAD_LIBRARY(libg4vertex.so)
//R__LOAD_LIBRARY(libglobalvertex.so);
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libmulticolstudy.so)
//R__LOAD_LIBRARY(libtrigzvtxchecker.so)
R__LOAD_LIBRARY(libcalotrigger.so)
//gSystem->Load("libg4detectors.so");
//gSystem->Load("libg4detectors.so");

bool file_exists(const char* filename)
{
  std::ifstream infile(filename);
  return infile.good();
}
int run_earlydata(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, string dir = ".", int issim = 0, int rn = 0)
{
  int verbosity = 0;//debug;
  string filename = dir+"/"+to_string(nproc)+"_multicol/events_"+tag+(tag==""?"":"_");
  if(!issim) filename += to_string(rn)+"_";
  filename += to_string(nproc)+"_";
  filename += to_string(nevt);
  filename += ".root";

  string tzfn = dir+"/"+to_string(nproc)+"_multicol/tz_"+tag+(tag==""?"":"_");
  if(!issim) tzfn += to_string(rn)+"_";
  tzfn += to_string(nproc)+"_";
  tzfn += to_string(nevt);
  tzfn += ".root";
  
  if(debug > 1) cout << "test1" << endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc =  recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG",issim?"MDC2":"ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP",issim?21:rn);
  
  se->Verbosity(verbosity);
  // just if we set some flags somewhere in this macro

  //Trigzvtxchecker* tz;
  //if(!issim) tz = new Trigzvtxchecker(tzfn, rn, nproc, debug, "tzvtx");
  //if(!issim) se->registerSubsystem(tz);

  Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  Fun4AllInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  Fun4AllInputManager *in_3 = new Fun4AllDstInputManager("DSTin3");
  Fun4AllInputManager *in_4 = new Fun4AllDstInputManager("DSTin4");
  Fun4AllInputManager *in_5 = new Fun4AllDstInputManager("DSTin5");
  cout << "get filenames" << endl;
  string line1, line2, line3, line4, line5;
  line1 = "./dsts/"+to_string(nproc)+"/calo_cluster_"+to_string(nproc)+".root";

  in_1->AddFile(line1);

  cout << "register managers" << endl;
  se->registerInputManager( in_1 );

  if(issim)
    {      
      line2 = "./dsts/"+to_string(nproc)+"/global_"+to_string(nproc)+".root";
      line3 = "./dsts/"+to_string(nproc)+"/mbd_epd_"+to_string(nproc)+".root";
      line4 = "./dsts/"+to_string(nproc)+"/truth_jet_"+to_string(nproc)+".root";
      line5 = "./dsts/"+to_string(nproc)+"/g4hits_"+to_string(nproc)+".root";
      in_2->AddFile(line2);
      in_3->AddFile(line3);
      in_4->AddFile(line4);
      in_5->AddFile(line5);
      //se->registerInputManager( in_2 );
      se->registerInputManager( in_3 );
      se->registerInputManager(in_4);
      se->registerInputManager(in_5);
    }
  std::cout << "status setters" << std::endl;

  CDBInterface::instance()->Verbosity(0);

  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR");
  te->setNSamples(12);
  te->setTriggerSample(6);
  te->setJetThreshold(8, 13, 17, 22);
  te->setPhotonThreshold(40, 56, 70, 78);
  // subrtraction delay of the post and pre sample                                                                     
  te->setTriggerDelay(5);
  te->SetIsData(true); // only for sim                                                                                
  te->Verbosity(0);
  te->setEmcalLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/emcal_ll1_lut_0.50tr_new.root");
  te->setHcalinLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/hcalin_ll1_lut_0.50tr_new.root");
  te->setHcaloutLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/hcalout_ll1_lut_0.50tr_new.root");
  //se->registerSubsystem(te);











  bool isSim = true;
  int data_sim_runnumber_thres = 1000;
  if (rc->get_uint64Flag("TIMESTAMP") > data_sim_runnumber_thres)
  {
    isSim = false;
  }
  std::cout << "Calo Calib uses runnumber " << rc->get_uint64Flag("TIMESTAMP") << std::endl;

  //////////////////////
  // Input geometry node
  std::cout << "Adding Geometry file" << std::endl;
  Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  ingeo->AddFile(geoLocation);
  se->registerInputManager(ingeo);

  //////////////////////////////
  // set statuses on raw towers
  std::cout << "status setters" << std::endl;
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  // MC Towers Status
  if(isSim) {
    // Uses threshold of 50% for towers be considered frequently bad.
    std::string calibName_hotMap = "CEMC_hotTowers_status";
    /* Systematic options (to be used as needed). */
    /* Uses threshold of 40% for towers be considered frequently bad. */
    // std::string calibName_hotMap = "CEMC_hotTowers_status_40";

    /* Uses threshold of 60% for towers be considered frequently bad. */
    // std::string calibName_hotMap = "CEMC_hotTowers_status_60";

    std::string calibdir = CDBInterface::instance()->getUrl(calibName_hotMap);
    statusEMC->set_directURL_hotMap(calibdir);
  }
  se->registerSubsystem(statusEMC);

  CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
  statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
  statusHCalIn->set_time_cut(2);
  se->registerSubsystem(statusHCalIn);

  CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
  statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
  statusHCALOUT->set_time_cut(2);
  se->registerSubsystem(statusHCALOUT);

  ////////////////////
  // Calibrate towers
  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  se->registerSubsystem(calibEMC);

  std::cout << "Calibrating OHcal" << std::endl;
  CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
  calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  se->registerSubsystem(calibOHCal);

  std::cout << "Calibrating IHcal" << std::endl;
  CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
  calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  se->registerSubsystem(calibIHCal);

  ////////////////
  // MC Calibration
  if (isSim)
  {
    std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
    if (MC_Calib.empty())
    {
      std::cout << "No MC calibration found :( )" << std::endl;
      gSystem->Exit(0);
    }
    CaloTowerCalib *calibEMC_MC = new CaloTowerCalib("CEMCCALIB_MC");
    calibEMC_MC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC_MC->set_inputNodePrefix("TOWERINFO_CALIB_");
    calibEMC_MC->set_outputNodePrefix("TOWERINFO_CALIB_");
    calibEMC_MC->set_directURL(MC_Calib);
    calibEMC_MC->set_doCalibOnly(true);
    se->registerSubsystem(calibEMC_MC);
  }










  

  JetReco *towerjetreco_noz = new JetReco("jetreco_noz");
  TowerJetInput* emtji_noz = new TowerJetInput(Jet::CEMC_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* ohtji_noz = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* ihtji_noz = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
  //towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  emtji_noz->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  ohtji_noz->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  ihtji_noz->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
      
  towerjetreco_noz->add_input(emtji_noz);
  towerjetreco_noz->add_input(ohtji_noz);
  towerjetreco_noz->add_input(ihtji_noz);
  towerjetreco_noz->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Tower_HIRecoSeedsRaw_r04_noz");
  towerjetreco_noz->set_algo_node("ANTIKT");
  towerjetreco_noz->set_input_node("TOWER");
  towerjetreco_noz->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco_noz);



  
  //auto mbddigi = new MbdDigitization();
  auto mbdreco = new MbdReco();
  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  
  //mbddigi->Verbosity(verbosity);
  //se->registerSubsystem(mbddigi);
  mbdreco->Verbosity(verbosity);
  se->registerSubsystem(mbdreco);
  
  gblvertex->Verbosity(verbosity);
  se->registerSubsystem(gblvertex);
  //Process_Calo_Calib();    
  se->Print("NODETREE");

  

  //TriggerRunInfoReco* tana = new TriggerRunInfoReco("tana");
  //se->registerSubsystem(tana);
  /*
  RetowerCEMC *rcemc = new RetowerCEMC();
  rcemc->set_towerinfo(true);
  rcemc->Verbosity(verbosity);
  se->registerSubsystem(rcemc);
  cout << "set up retower emcal" << endl;
  */
  /*
  JetReco *truthjetreco = new JetReco("tjetreco");
  TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
  tji->add_embedding_flag(0);  // changes depending on signal vs. embedded
  truthjetreco->add_input(tji);
  truthjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  truthjetreco->set_algo_node("ANTIKT");
  truthjetreco->set_input_node("TRUTH");
  if(issim) se->registerSubsystem(truthjetreco);
  */
  
  JetReco *towerjetreco = new JetReco("jetreco");
  TowerJetInput* emtji = new TowerJetInput(Jet::CEMC_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* ohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* ihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
  //towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  emtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  ohtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  ihtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
      
  towerjetreco->add_input(emtji);
  towerjetreco->add_input(ohtji);
  towerjetreco->add_input(ihtji);
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Tower_HIRecoSeedsRaw_r04");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);

  multiColStudy* mcs = new multiColStudy(filename, "mcs", debug, issim);
  se->registerSubsystem(mcs);

  
  cout << "test4" << endl;
  se->Print("NODETREE");
  cout << "run " << nevt << endl;
  se->skip(0);
  se->run(nevt);
  cout << "ran " << nevt << endl;
  cout << "Ran all events" << endl;
  se->Print("NODETREE");
  cout << "ending" << endl;
  se->End();
  cout << "printing timer" << endl;
  se->PrintTimer();
  cout << "Ended server" << endl;
  delete se;
  cout << "Deleted server" << endl;
  //cout << "wrote " << filename << endl;
  //gSystem->Exit(0);
  //cout << "Exited gSystem" << endl;
  return 0;

}
