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

int rnvals[42] = {49121,49125,49127,49128,49133,49136,49138,49140,49141,49143,49148,49157,49159,49165,49217,49218,49219,49224,49226,49227,49228,49229,49230,49233,49240,49241,49244,49245,49247,49248,49249,49250,49251,49254,49262,49263,49264,49265,49266,49267,49268,49270};

double lumivals[42] ={0.0567942,0.0516462,0.0642938,0.0354406,0.0143268,0.0412767,0.0238550,0.0116476,0.0245767,0.0052713,0.0169815,0.0210836,0.0154796,0.0065756,0.0309800,0.0160149,0.0832480,0.0123099,0.0130056,0.0192810,0.0217280,0.0598882,0.0129254,0.0044402,0.0409331,0.0264116,0.0432918,0.0188073,0.0003522,0.0666371,0.0388263,0.0269200,0.0107077,0.0131979,0.0228682,0.0243940,0.0414531,0.0592927,0.0187420,0.0446407,0.0159728,0.0176773};











int make_hists(string tag, vector<int> rns, vector<int> nfiles, int triggerbit = 18, int dodijetcut = 1, int jetcut = 1, int issim = 0)
{
  //define histograms
  string region = "";
  if(rns[0] < 49216) region = "RegionA";
  else if(rns[0] < 49238) region = "RegionB";
  else if(rns[0] < 49256) region = "RegionC";
  else region = "RegionD";
  TH2D* jet_eta_phi = new TH2D(("jet_eta_phi_"+region).c_str(),"",30,-1.5,1.5,64,-M_PI,M_PI);
  TH2D* jet_eta_phi_gr20 = new TH2D(("jet_eta_phi_gr20_"+region).c_str(),"",30,-1.5,1.5,64,-M_PI,M_PI);
  TH2D* tow_eta_phi = new TH2D(("tow_eta_phi_"+region).c_str(),"",30,-1.5,1.5,64,-M_PI,M_PI);
  TH2D* tow_eta_phi_gr5 = new TH2D(("tow_eta_phi_gr5_"+region).c_str(),"",30,-1.5,1.5,64,-M_PI,M_PI);
  TH2D* tow_deteta_phi = new TH2D(("tow_deteta_phi_"+region).c_str(),"",24,-1.1,1.1,64,-M_PI,M_PI);
  TH2D* tow_deteta_phi_gr5 = new TH2D(("tow_deteta_phi_gr5_"+region).c_str(),"",24,-1.1,1.1,64,-M_PI,M_PI);
  TH2D* jet_E_eta = new TH2D(("jet_E_eta_"+region).c_str(),"",60,0,60,30,-1.5,1.5);
  TH2D* tow_E_eta = new TH2D(("tow_E_eta_"+region).c_str(),"",30,0,15,30,-1.5,1.5);
  TH2D* tow_E_deteta_em = new TH2D(("tow_E_deteta_em_"+region).c_str(),"",30,0,15,30,-1.5,1.5);
  TH2D* tow_E_deteta_oh = new TH2D(("tow_E_deteta_oh_"+region).c_str(),"",30,0,15,30,-1.5,1.5);
  TH2D* tow_E_deteta = new TH2D(("tow_E_deteta_"+region).c_str(),"",30,0,15,24,-1.1,1.1);
  TH2D* jet_E_phi = new TH2D(("jet_E_phi_"+region).c_str(),"",60,0,60,64,-M_PI,M_PI);
  TH2D* tow_E_phi = new TH2D(("tow_E_phi_"+region).c_str(),"",30,0,15,64,-M_PI,M_PI);
  TH2D* jet_ET_eta = new TH2D(("jet_ET_eta_"+region).c_str(),"",60,0,60,30,-1.5,1.5);
  TH2D* tow_ET_eta = new TH2D(("tow_ET_eta_"+region).c_str(),"",30,0,15,30,-1.5,1.5);
  TH2D* tow_ET_deteta = new TH2D(("tow_ET_deteta_"+region).c_str(),"",30,0,15,24,-1.1,1.1);
  TH2D* jet_ET_phi = new TH2D(("jet_ET_phi_"+region).c_str(),"",60,0,60,64,-M_PI,M_PI);
  TH2D* tow_ET_phi = new TH2D(("tow_ET_phi_"+region).c_str(),"",30,0,15,64,-M_PI,M_PI);
  TH1D* zhist = new TH1D(("zhist_"+region).c_str(),"",400,-200,200);
  TH1D* zhist_gr20 = new TH1D(("zhist_gr20_"+region).c_str(),"",400,-200,200);
  TH1D* mbdn = new TH1D(("mbdn_"+region).c_str(),"",100,0,10);
  TH1D* mbds = new TH1D(("mbds_"+region).c_str(),"",100,0,10);
  TH1D* mbdt = new TH1D(("mbdt_"+region).c_str(),"",100,0,10);
  TH2D* jet_E_frcem = new TH2D(("jet_E_frcem_"+region).c_str(),"",60,0,60,24,-0.1,1.1);
  TH2D* frcem_frcoh = new TH2D(("frcem_frcoh_"+region).c_str(),"",24,-0.1,1.1,24,-0.1,1.1);
  TH2D* frcem_frcoh_gr20 = new TH2D(("frcem_frcoh_gr20_"+region).c_str(),"",24,-0.1,1.1,24,-0.1,1.1);
  TH2D* jet_ET_dphi = new TH2D(("jet_ET_dphi_"+region).c_str(),"",60,0,60,32,0,M_PI);
  TH1D* calo_hitsgrone[3];
  TH2D* calo_tgrone_eta[3];
  TH2D* jet_eta_frcem = new TH2D(("jet_eta_frcem_"+region).c_str(),"",30,-1.5,1.5,24,-0.1,1.1);
  TH2D* jet_eta_frcem_gr20 = new TH2D(("jet_eta_frcem_gr20_"+region).c_str(),"",30,-1.5,1.5,24,-0.1,1.1);
  TH2D* jet_phi_frcem = new TH2D(("jet_phi_frcem_"+region).c_str(),"",44,-M_PI,M_PI,24,-0.1,1.1);
  TH2D* jet_phi_frcem_gr20 = new TH2D(("jet_phi_frcem_gr20_"+region).c_str(),"",44,-M_PI,M_PI,24,-0.1,1.1);
  TH1D* zhist_nocut = new TH1D(("zhist_nocut_"+region).c_str(),"",400,-200,200);
  TH2D* jet_t_frcem = new TH2D(("jet_t_frcem_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH2D* jet_t_frcem_gr20 = new TH2D(("jet_t_frcem_gr20_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH2D* emat_ohat = new TH2D(("emat_ohat_"+region).c_str(),"",60,-6,6,22,-1.1,1.1);
  TH1D* h_emat = new TH1D(("h_emat_"+region).c_str(),"",60,-6,6);
  TH1D* h_ohat = new TH1D(("h_ohat_"+region).c_str(),"",60,-6,6);
  TH2D* emat_calo_emfrac = new TH2D(("emat_calo_emfrac_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH2D* jet_at_em_at_oh = new TH2D(("jet_at_em_at_oh_"+region).c_str(),"",60,-6,6,60,-6,6);
  TH2D* jet_at_em_at_oh_gr20 = new TH2D(("jet_at_em_at_oh_gr20_"+region).c_str(),"",60,-6,6,60,-6,6);
  TH2D* jet_at_em_frcem = new TH2D(("jet_at_em_frcem_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH2D* jet_at_em_eta = new TH2D(("jet_at_em_eta_"+region).c_str(),"",60,-6,6,22,-1.1,1.1);
  TH2D* jet_at_oh_eta = new TH2D(("jet_at_oh_eta_"+region).c_str(),"",60,-6,6,22,-1.1,1.1);
  TH2D* jet_at_em_frcem_gr20 = new TH2D(("jet_at_em_frcem_gr20_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH2D* jet_at_oh_frcem = new TH2D(("jet_at_oh_frcem_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH2D* jet_at_oh_frcem_gr20 = new TH2D(("jet_at_oh_frcem_gr20_"+region).c_str(),"",60,-6,6,24,-0.1,1.1);
  TH1D* njet_lumi = new TH1D("njet_lumi","",49280-49120,49120-0.5,49280-0.5);
  for(int i=0; i<3; ++i)
    {
      calo_hitsgrone[i] = new TH1D(("calo_hitsgrone_"+to_string(i)+"_"+region).c_str(),"",20,-0.5,19.5);
      calo_tgrone_eta[i] = new TH2D(("calo_tgrone_eta_"+to_string(i)+"_"+region).c_str(),"",60,-6,6,22,-1.1,1.1);
    }

  const int ntrigtypes = 4;
  long long unsigned int trigs[2][ntrigtypes] = {0};
  trigs[0][0]=10;
  trigs[0][1]=18;
  trigs[0][2]=26;
  trigs[0][3]=30;
  //define constants
  const int nmaxjet = 10;
  const int nzvtx = 3;
  const int maxtow = 24576+1536*2;
  const int ntowfield = 6;
  const int mbdchan = 64;
  const int mbdside = 2;
  int nzvtx30 = 0;
  int nzvtxany = 0;
  int nzvtx30ndj = 0;
  int nzvtx30jet = 0;
  int nzvtxanyndj = 0;
  int nzvtxanyjet = 0;
  int nzvtx30jc = 0;
  int nzvtxanyjc = 0;
  //start processing
  long long unsigned int totalentries = 0;
  if(triggerbit==18)
    {
      dodijetcut = 1;
      jetcut = 1;
    }
  else if(triggerbit == 10)
    {
      dodijetcut = 0;
      jetcut = 0;
    }
  for(int g = 0; g<rns.size(); ++g)
    {
      int rn = rns.at(g);
      double lumi = 0;
      for(int i=0; i<42; ++i)
	{
	  if(rn == rnvals[i])
	    {
	      lumi = lumivals[i];
	    }
	}
      int nfile = nfiles.at(g);
  for(int h=0; h<nfile; ++h)
    {
      string filename = "/sphenix/tg/tg01/jets/jocl/multiCol/"+to_string(rn)+"/events_"+tag+"_"+to_string(rn)+"_"+to_string(h)+"_0.root";
      cout << "Processing file " << filename << endl;
      TFile* datfile = TFile::Open(filename.c_str());
      //define branch variables
      if(!datfile) continue;
      int njet, hitsgrone, isdijet;
      double dphilead;
      double jet_e[nmaxjet];
      double jet_eta[nmaxjet];
      double jet_phi[nmaxjet];
      double zvtx[nzvtx];
      double towgrone[maxtow][ntowfield];
      double mbdq[mbdside][mbdchan];
      double frcem[nmaxjet];
      double frcoh[nmaxjet];
      double jet_at[nmaxjet];
      double jet_at_oh[nmaxjet];
      double jet_at_em[nmaxjet];
      long long unsigned int trigvec;
      //define non-branch variables
      double mbdnq, mbdsq, mbdtq;
      double calo_emfrac;//, calo_ohfrac, calo_e;
      double emat[nmaxjet][64];
      double ohat[nmaxjet][64];
      int ncgroe[nmaxjet];
      int ncgroo[nmaxjet];
      //get TTree
      TTree* dattree = (TTree*)datfile->Get("tree");
      //set up branches
      dattree->SetBranchAddress("njet",&njet);
      dattree->SetBranchAddress("hitsgrone",&hitsgrone);
      dattree->SetBranchAddress("isdijet",&isdijet);
      dattree->SetBranchAddress("trigvec",&trigvec);
      dattree->SetBranchAddress("jet_e",jet_e);
      dattree->SetBranchAddress("jet_eta",jet_eta);
      dattree->SetBranchAddress("jet_phi",jet_phi);
      dattree->SetBranchAddress("rzvtx",zvtx);
      dattree->SetBranchAddress("towgrone",towgrone);
      //for above: 0 = E, 1 = eta, 2 = phi, 3 = detector level eta, 4 = calo
      dattree->SetBranchAddress("mbdq",mbdq);
      dattree->SetBranchAddress("frcem",frcem);
      dattree->SetBranchAddress("frcoh",frcoh);
      dattree->SetBranchAddress("dphilead",&dphilead);
      dattree->SetBranchAddress("jet_at",jet_at);
      dattree->SetBranchAddress("emat",emat);
      dattree->SetBranchAddress("ohat",ohat);
      dattree->SetBranchAddress("calo_emfrac",&calo_emfrac);
      dattree->SetBranchAddress("jet_at_em",jet_at_em);
      dattree->SetBranchAddress("jet_at_oh",jet_at_oh);
      dattree->SetBranchAddress("ncgroe",ncgroe);
      dattree->SetBranchAddress("ncgroo",ncgroo);
      
      
      for(int i=0; i<dattree->GetEntries(); ++i)
	{
	  dattree->GetEntry(i);
	  for(int j=0; j<ntrigtypes; ++j)
	    {
	      trigs[1][j] += (trigvec >> trigs[0][j]) & 1;
	    }
	  if(!((trigvec >> triggerbit) & 1)) continue;
	  ++totalentries;
	  if(!std::isnan(zvtx[0]))
	    {
	      ++nzvtxanyndj;
	      if(abs(zvtx[0]) < 30) ++nzvtx30ndj;
	    }
	  
	  if(dodijetcut && (dphilead < 3*M_PI/4 || isdijet == 0)) continue;
	  double ETmax = 0;
	  double ETsub = 0;
	  for(int j=0; j<njet; ++j)
	    {
	      double ET = jet_e[j]/cosh(jet_eta[j]);
	      if(ET > ETmax)
		{
		  ETsub = ETmax;
		  ETmax = ET;
		}
	      else if(ET > ETsub) ETsub = ET;
	    }

	  if(dodijetcut && (isdijet == 0 || ETsub/ETmax < 0.3)) continue;
	  if(!std::isnan(zvtx[0]))
	    {
	      zhist_nocut->Fill(zvtx[0]);
	      ++nzvtxany;
	      if(abs(zvtx[0]) < 30) ++nzvtx30;
	    }
	  if(jetcut && ETmax < 15) continue;
	  mbdnq = 0;
	  mbdsq = 0;
	  mbdtq = 0;
 	  if(!std::isnan(zvtx[0]))
	    {
	      ++nzvtxanyjc;
	      if(abs(zvtx[0]) < 30) ++nzvtx30jc;
	    }
	  //emat_calo_emfrac->Fill(emat,calo_emfrac);
	  for(int j=0; j<njet; ++j)
	    {
	      double ET = jet_e[j]/cosh(jet_eta[j]);
	      if(ET < 15) continue;
	      if(!std::isnan(zvtx[0]))
		{
		  ++nzvtxanyjet;
		  if(abs(zvtx[0]) < 30)
		    {
		      ++nzvtx30jet;
		      if(lumi != 0 && triggerbit == 18) njet_lumi->Fill(rn,1./lumi);
		    }
		}
	      jet_at_em_at_oh->Fill(jet_at_em[j],jet_at_oh[j]);
	      jet_at_em_frcem->Fill(jet_at_em[j],frcem[j]);
	      jet_at_em_eta->Fill(jet_at_em[j],jet_eta[j]);
	      jet_at_oh_eta->Fill(jet_at_oh[j],jet_eta[j]);
	      jet_at_oh_frcem->Fill(jet_at_oh[j],frcem[j]);
	      jet_t_frcem->Fill(jet_at[j],frcem[j]);
	      jet_eta_phi->Fill(jet_eta[j],jet_phi[j]);
	      jet_E_eta->Fill(jet_e[j],jet_eta[j]);
	      jet_E_phi->Fill(jet_e[j],jet_phi[j]);
	      jet_ET_eta->Fill(ET,jet_eta[j]);
	      jet_ET_phi->Fill(ET,jet_phi[j]);
	      jet_E_frcem->Fill(jet_e[j],frcem[j]);
	      frcem_frcoh->Fill(frcem[j],frcoh[j]);
	      jet_eta_frcem->Fill(jet_eta[j],frcem[j]);
	      jet_phi_frcem->Fill(jet_phi[j],frcem[j]);
	      cout << ncgroe[j] << ": ";
	      for(int k = 0; k<ncgroe[j]; ++k)
		{
		  h_emat->Fill(emat[j][k]);
		  cout << emat[j][k] << " ";
		}
	      cout << endl;
	      for(int k=0; k<ncgroo[j]; ++k)
		{
		  h_ohat->Fill(ohat[j][k]);
		}
	      if(ET > 20)
		{
		  jet_at_em_frcem_gr20->Fill(jet_at_em[j],frcem[j]);
		  jet_at_oh_frcem_gr20->Fill(jet_at_oh[j],frcem[j]);
		  jet_at_em_at_oh_gr20->Fill(jet_at_em[j],jet_at_oh[j]);
		  jet_t_frcem_gr20->Fill(jet_at[j],frcem[j]);
		  frcem_frcoh_gr20->Fill(frcem[j],frcoh[j]);
		  jet_eta_frcem_gr20->Fill(jet_eta[j],frcem[j]);
		  jet_phi_frcem_gr20->Fill(jet_phi[j],frcem[j]);
		  jet_eta_phi_gr20->Fill(jet_eta[j],jet_phi[j]);
		}
	    }
	  if(isdijet) jet_ET_dphi->Fill(ETmax,dphilead);
	  if(ETmax > 20)
	    {
	      if(!std::isnan(zvtx[0]))
		{
		  zhist_gr20->Fill(zvtx[0]);
		}
	    }
	  for(int j=0; j<mbdside; ++j)
	    {
	      for(int k=0; k<mbdchan; ++k)
		{
		  mbdtq += mbdq[j][k];
		  if(j==0) mbdn->Fill(mbdq[j][k]);
		  else mbds->Fill(mbdq[j][k]);
		}
	    }
	  mbdt->Fill(mbdtq);
	  int calohits[3] = {0};
	  for(int j=0; j<hitsgrone; ++j)
	    {
	      int hitcalo = -1;
	      if(towgrone[j][4] < 2.5)
		{
		  if(towgrone[j][4] > 1.5) hitcalo = 2;
		  else if(towgrone[j][4] > 0.5) hitcalo = 1;
		  else if(towgrone[j][4] > -0.5) hitcalo = 0;
		}
	      if(hitcalo != -1) ++calohits[hitcalo];
	      else continue;
	      
	      double ET = towgrone[j][0]/cosh(towgrone[j][1]);
	      tow_eta_phi->Fill(towgrone[j][1],towgrone[j][2]);

	      calo_tgrone_eta[hitcalo]->Fill(towgrone[j][5],towgrone[j][3]);
	      
	      tow_deteta_phi->Fill(towgrone[j][3],towgrone[j][2]);
	      tow_E_eta->Fill(towgrone[j][0],towgrone[j][1]);
	      tow_E_deteta->Fill(towgrone[j][0],towgrone[j][3]);
	      if(hitcalo == 0) tow_E_deteta_em->Fill(towgrone[j][0],towgrone[j][3]);
	      else if(hitcalo == 2) tow_E_deteta_oh->Fill(towgrone[j][0],towgrone[j][3]);
	      tow_E_phi->Fill(towgrone[j][0],towgrone[j][2]);
	      tow_ET_eta->Fill(ET,towgrone[j][1]);
	      tow_ET_phi->Fill(ET,towgrone[j][2]);
	      tow_ET_deteta->Fill(ET,towgrone[j][3]);
	      if(towgrone[j][0] > 5)
		{
		  tow_eta_phi_gr5->Fill(towgrone[j][1],towgrone[j][2]);
		  tow_deteta_phi_gr5->Fill(towgrone[j][1],towgrone[j][2]);
		}
	    }
	  for(int j=0; j<3; ++j)
	    {
	      calo_hitsgrone[j]->Fill(calohits[j]);
	    }
	  if(!std::isnan(zvtx[0]))
	    {
	      zhist->Fill(zvtx[0]);
	    }
	}
      datfile->Close();
    }
  }
  TFile* outf = TFile::Open(("/sphenix/user/jocl/projects/multiColStudy/output/hists/hists_"+tag+"_"+(dodijetcut?"dc":"nc")+"_"+region+"_"+to_string(triggerbit)+".root").c_str(),"RECREATE");
  outf->cd();


  for(int i=0; i<3; ++i)
    {
      calo_hitsgrone[i]->Scale(1./totalentries);
      calo_tgrone_eta[i]->Scale(1./totalentries);
    }
  jet_eta_phi->Scale(1./totalentries);
  jet_eta_phi_gr20->Scale(1./totalentries);
  tow_eta_phi->Scale(1./totalentries);
  tow_eta_phi_gr5->Scale(1./totalentries);
  tow_deteta_phi->Scale(1./totalentries);
  tow_deteta_phi_gr5->Scale(1./totalentries);
  jet_E_eta->Scale(1./totalentries);
  tow_E_eta->Scale(1./totalentries);
  tow_E_deteta->Scale(1./totalentries);
  tow_E_deteta_em->Scale(1./totalentries);
  tow_E_deteta_oh->Scale(1./totalentries);
  jet_E_phi->Scale(1./totalentries);
  tow_E_phi->Scale(1./totalentries);
  jet_ET_eta->Scale(1./totalentries);
  tow_ET_eta->Scale(1./totalentries);
  tow_ET_deteta->Scale(1./totalentries);
  jet_ET_phi->Scale(1./totalentries);
  tow_ET_phi->Scale(1./totalentries);
  zhist->Scale(1./totalentries);
  zhist_nocut->Scale(1./totalentries);
  mbdn->Scale(1./totalentries);
  mbds->Scale(1./totalentries);
  mbdt->Scale(1./totalentries);
  jet_E_frcem->Scale(1./totalentries);
  jet_ET_dphi->Scale(1./totalentries);
  jet_eta_frcem->Scale(1./totalentries);
  jet_phi_frcem->Scale(1./totalentries);
  jet_eta_frcem_gr20->Scale(1./totalentries);
  jet_phi_frcem_gr20->Scale(1./totalentries);
  zhist_gr20->Scale(1./totalentries);
  frcem_frcoh->Scale(1./totalentries);
  frcem_frcoh_gr20->Scale(1./totalentries);
  jet_t_frcem->Scale(1./totalentries);
  jet_t_frcem_gr20->Scale(1./totalentries);
  emat_ohat->Scale(1./totalentries);
  h_emat->Scale(1./totalentries);
  h_ohat->Scale(1./totalentries);
  emat_calo_emfrac->Scale(1./totalentries);
  jet_at_em_at_oh->Scale(1./totalentries);
  jet_at_em_eta->Scale(1./totalentries);
  jet_at_oh_eta->Scale(1./totalentries);
  jet_at_em_at_oh_gr20->Scale(1./totalentries);
  jet_at_em_frcem->Scale(1./totalentries);
  jet_at_oh_frcem->Scale(1./totalentries);
  jet_at_em_frcem_gr20->Scale(1./totalentries);
  jet_at_oh_frcem_gr20->Scale(1./totalentries);
  for(int i=0; i<3; ++i)
    {
      calo_hitsgrone[i]->Write();
      calo_tgrone_eta[i]->Write();
    }
  njet_lumi->Write();
  h_emat->Write();
  h_ohat->Write();
  jet_eta_phi->Write();
  jet_eta_phi_gr20->Write();
  tow_eta_phi->Write();
  tow_eta_phi_gr5->Write();
  tow_deteta_phi->Write();
  tow_deteta_phi_gr5->Write();
  jet_E_eta->Write();
  tow_E_eta->Write();
  tow_E_deteta->Write();
  tow_E_deteta_em->Write();
  tow_E_deteta_oh->Write();
  jet_E_phi->Write();
  tow_E_phi->Write();
  jet_ET_eta->Write();
  tow_ET_eta->Write();
  tow_ET_deteta->Write();
  jet_ET_phi->Write();
  tow_ET_phi->Write();
  zhist->Write();
  zhist_nocut->Write();
  mbdn->Write();
  mbds->Write();
  mbdt->Write();
  jet_E_frcem->Write();
  jet_ET_dphi->Write();
  jet_eta_frcem->Write();
  jet_phi_frcem->Write();
  jet_eta_frcem_gr20->Write();
  jet_phi_frcem_gr20->Write();
  zhist_gr20->Write();
  frcem_frcoh->Write();
  frcem_frcoh_gr20->Write();
  jet_t_frcem->Write();
  jet_t_frcem_gr20->Write();
  emat_ohat->Write();
  emat_calo_emfrac->Write();
  jet_at_em_at_oh->Write();
  jet_at_em_frcem->Write();
  jet_at_oh_frcem->Write();
  jet_at_em_at_oh_gr20->Write();
  jet_at_em_frcem_gr20->Write();
  jet_at_oh_frcem_gr20->Write();
  jet_at_em_eta->Write();
  jet_at_oh_eta->Write();
  

  ofstream outtrigs("trigcounts/outtrigs"+region+".txt");
  for(int i=0; i<ntrigtypes; ++i)
    {
      outtrigs << trigs[0][i] << " " << trigs[1][i] << endl;
    }
    
  outf->Write();
  outf->Close();
  cout << "Njet zvtx < 30: " <<  nzvtx30jet << endl;
  cout << "Njet any zvtx:  " << nzvtxanyjet << endl;
  cout << "Nevt dj z < 30: " << nzvtx30 << endl;
  cout << "Nevt dj any z:  " << nzvtxany << endl;
  cout << "Nevt ndj any z: " << nzvtxanyndj << endl;
  cout << "Nevt ndj z<30:  " << nzvtx30ndj << endl;
  cout << "Nevt jc any z:  " << nzvtxanyjc << endl;
  cout << "Nevt jc z < 30: " << nzvtx30jc << endl;
  return 0;
}
