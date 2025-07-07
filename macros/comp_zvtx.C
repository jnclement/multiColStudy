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

bool compfirst(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0] > b[0];
}

vector<vector<double>> make_jet_vector(int njet, double* jet_pt, double* jet_eta, double* jet_phi, int istruth, double zvtx, int sampletype, double* jet_e)
{
  vector<vector<double>> jet_vector = {};
  //cout << endl << endl << "jet pt: " << endl;
  for(int i=0; i<njet; ++i)
    {
      //if(istruth && jet_pt[i] < (sampletype==1?30:15)) continue;
      //if(istruth && sampletype != 1 && jet_pt[i] > 30) continue;
      if(check_bad_jet_eta(jet_eta[i],zvtx,(istruth?0.4:0.4))) continue;
      vector<double> jet = {jet_pt[i],jet_eta[i],jet_phi[i],jet_e[i]};
      jet_vector.push_back(jet);
    }
  //cout << endl << endl;
  if(jet_vector.size() < 1) return {};
  std::sort(jet_vector.begin(),jet_vector.end(),compfirst);
  if(istruth && sampletype == 4)
    {
      if(jet_vector.at(0).at(0) < 52) return {};
    }
  if(istruth && sampletype == 3)
    {
      if(jet_vector.at(0).at(0) < 35 || jet_vector.at(0).at(0) > 52) return {};
    }
  if(istruth && sampletype == 2)
    {
      if(jet_vector.at(0).at(0) < 22 || jet_vector.at(0).at(0) > 35) return {};
    }
  if(istruth && sampletype == 1)
    {
      if(jet_vector.at(0).at(0) > 22 || jet_vector.at(0).at(0) < 14) return {};
    }
  if(istruth && sampletype == 0)
    {
      if(jet_vector.at(0).at(0) > 14) return {};
    }

  return jet_vector;
}

vector<vector<double>> truth_match(vector<vector<double>> truth_jets, vector<vector<double>> reco_jets, TEfficiency* eff)
{
  vector<vector<double>> matches = {};
  for(int i=0; i<truth_jets.size(); ++i)
    {
      //cout << truth_jets.at(i).at(0) << " " << truth_jets.at(i).at(1) << " " << truth_jets.at(i).at(2) << endl;

      int match_index = -1;
      double max_pt = 0;
      vector<double> match = {};
      for(int j=0; j<reco_jets.size(); ++j)
	{
	  double dPhi = abs(truth_jets.at(i).at(2) - reco_jets.at(j).at(2));
	  //cout << "    " << reco_jets.at(j).at(0) << " " << reco_jets.at(j).at(1) << " " << reco_jets.at(j).at(2) << endl;
	  if((abs(dPhi) < 0.3 || 2*M_PI - abs(dPhi) < 0.3) && reco_jets.at(j).at(0) > max_pt)
	    {
	      max_pt = reco_jets.at(j).at(0);
	      match_index = j;
	    }
	}
      if(match_index == -1)
	{
	  eff->Fill(false,truth_jets.at(i).at(0));
	  continue;
	}
      eff->Fill(true,truth_jets.at(i).at(0));
      //cout << "here" << endl;
      match.push_back(truth_jets.at(i).at(0));
      match.push_back(reco_jets.at(match_index).at(0));
      match.push_back(truth_jets.at(i).at(3));
      match.push_back(reco_jets.at(match_index).at(3));
      match.push_back(truth_jets.at(i).at(1));
      match.push_back(reco_jets.at(match_index).at(1));
      matches.push_back(match);
      reco_jets.erase(reco_jets.begin()+match_index);
    }
  return matches;
}

int comp_zvtx(string tag, int rn, int sampletype = 0)
{
  double scalefactor = 4.197e-3;
  if(sampletype == 1) scalefactor = 3.997e-6;
  if(sampletype == 2) scalefactor = 6.218e-8;
  if(sampletype == 3) scalefactor = 2.502e-9;
  if(sampletype == 4) scalefactor = 7.2695e-12;
  const int nmaxjet = 100;
  const int nzvtx = 3;
  TH3D* h3_resp_pT_zvtx = new TH3D("h3_resp_pT_zvtx",";p_{T}^{reco}/p_{T}^{truth};p_{T}^{truth} [GeV];z_{vtx} [cm]",200,0,2,100,0,100,300,-150,150);
  TH3D* h3_resp_pT_zvtx_noz = new TH3D("h3_resp_pT_zvtx_noz",";p_{T}^{reco}/p_{T}^{truth};p_{T}^{truth} [GeV];z_{vtx} [cm]",200,0,2,100,0,100,300,-150,150);
  TH3D* h3_teta_ptd_tz = new TH3D("h3_teta_ptd_tz",";#eta_{jet}^{truth};p_{T,w/z}^{reco} - p_{T,noz}^{reco} [GeV];z_{vtx}^{truth} [cm]",40,-2,2,100,-50,50,300,-150,150);
  TH2D* etapt = new TH2D("etapt",";p_{T}^{truth};Jet #eta",15,10,85,30,-1.5,1.5);
  
  TH3D* h3_resp_E_zvtx = new TH3D("h3_resp_E_zvtx",";E^{reco}/E^{truth};E^{truth} [GeV];z_{vtx} [cm]",200,0,2,100,0,100,300,-150,150);
  TH3D* h3_resp_E_zvtx_noz = new TH3D("h3_resp_E_zvtx_noz",";E^{reco}/E^{truth};E^{truth} [GeV];z_{vtx} [cm]",200,0,2,100,0,100,300,-150,150);
  
  TEfficiency* eff_wz = new TEfficiency("eff_wz",";p_{T} [GeV];Matching Efficiency",25,0,100);
  TEfficiency* eff_nz = new TEfficiency("eff_nz",";p_{T} [GeV];Matching Efficiency",25,0,100);
  TH3D* noz_recoz_corrET = new TH3D("noz_recoz_corrEt",";p_{T,w/z}^{reco} - p_{T,noz}^{reco} [GeV];(p_{T,noz}^{reco} + p_{T,w/z}^{reco})/2 [GeV];z_{vtx} cm",100,-50,50,20,0,100,300,-150,150);
  TH2D* temp2 = new TH2D("temp2",";(p_{T,noz}^{reco} + p_{T,w/z}^{reco})/2 [GeV];p_{T,w/z}^{reco} - p_{T,noz}^{reco} [GeV]",20,0,100,100,-50,50);
  for(int h=rn*100; h<rn*100+100; ++h)
    {
      string filename = "/sphenix/tg/tg01/jets/jocl/multiCol/"+to_string(h)+"/events_"+tag+"_"+to_string(h)+"_0.root";
      //cout << "Processing file " << filename << endl;
      TFile* datfile = TFile::Open(filename.c_str());
      
      if(!datfile) continue;
      TTree* tree = (TTree*)datfile->Get("tree");
      if(!tree) continue;
      int njet, tnjet, njet_noz;
      double jet_pt[nmaxjet], jet_eta[nmaxjet], jet_pt_noz[nmaxjet], jet_eta_noz[nmaxjet], tjet_pt[nmaxjet], tjet_eta[nmaxjet], jet_phi[nmaxjet], tjet_phi[nmaxjet], jet_phi_noz[nmaxjet], tzvtx[nzvtx], rzvtx[nzvtx], jet_e[nmaxjet], jet_e_noz[nmaxjet], tjet_e[nmaxjet];
      
      tree->SetBranchAddress("njet",&njet);
      tree->SetBranchAddress("tnjet",&tnjet);
      tree->SetBranchAddress("njet_noz",&njet_noz);
      tree->SetBranchAddress("jet_et",jet_pt);
      tree->SetBranchAddress("jet_eta",jet_eta);
      tree->SetBranchAddress("tjet_et",tjet_pt);
      tree->SetBranchAddress("tjet_eta",tjet_eta);
      tree->SetBranchAddress("jet_pt_noz",jet_pt_noz);
      tree->SetBranchAddress("jet_eta_noz",jet_eta_noz);
      tree->SetBranchAddress("jet_phi",jet_phi);
      tree->SetBranchAddress("tjet_phi",tjet_phi);
      tree->SetBranchAddress("jet_phi_noz",jet_phi_noz);
      tree->SetBranchAddress("tzvtx",tzvtx);
      tree->SetBranchAddress("rzvtx",rzvtx);
      tree->SetBranchAddress("jet_e_noz",jet_e_noz);
      tree->SetBranchAddress("jet_e",jet_e);
      tree->SetBranchAddress("tjet_e",tjet_e);

      
      for(int i=0; i<tree->GetEntries(); ++i)
	{
	  tree->GetEntry(i);
	  vector<vector<double>> truthjet = make_jet_vector(tnjet, tjet_pt, tjet_eta, tjet_phi,1,tzvtx[0],sampletype,tjet_e);
	  if(truthjet.size() == 0) continue;
	  for(int j=0; j<truthjet.size(); ++j)
	    {
	      etapt->Fill(truthjet.at(j).at(0),truthjet.at(j).at(1),scalefactor);
	    }
	  //cout << "make recojets" << endl;
	  vector<vector<double>> recojets = make_jet_vector(njet, jet_pt, jet_eta, jet_phi,0,rzvtx[0],sampletype,jet_e);
	  //cout << endl<<endl<< "tz/z: " << tzvtx[0] <<" " << rzvtx[0] << endl;
	  //cout <<"recojets:" << endl;
	  for(int j=0; j<recojets.size(); ++j)
	    {
	      //cout << recojets.at(j).at(0) << " " << recojets.at(j).at(1) << endl;
	    }
	  //cout << "make truthjets" << endl;

	  //cout << "make reco noz" << endl;
	  vector<vector<double>> reco_noz = make_jet_vector(njet_noz, jet_pt_noz, jet_eta_noz, jet_phi_noz,0,0,sampletype,jet_e_noz);
	  //cout << "recojets noz:" << endl;
	  for(int j=0; j<reco_noz.size(); ++j)
	    {
	      //cout << reco_noz.at(j).at(0) << " " << reco_noz.at(j).at(1) << endl;
	    }
	  //cout << "make matches" << endl;
	  vector<vector<double>> matches = truth_match(truthjet, recojets, eff_wz);
	  //cout << "make matches noz" << endl;
	  vector<vector<double>> matches_noz = truth_match(truthjet, reco_noz, eff_nz);
	  //cout << "fill" << endl;
	  for(int j=0; j<matches.size(); ++j)
	    {
	      //	      cout << "enter matches"<< endl;
	      for(int k=0; k<matches_noz.size(); ++k)
		{
		  if(abs(matches.at(j).at(0) - matches_noz.at(k).at(0)) < 1e-6)
		    {
		      noz_recoz_corrET->Fill(matches.at(j).at(1)-matches_noz.at(k).at(1),(matches.at(j).at(1)+matches_noz.at(k).at(1))/2,tzvtx[0],scalefactor);
		      temp2->Fill((matches.at(j).at(1)+matches_noz.at(k).at(1))/2,matches.at(j).at(1)-matches_noz.at(k).at(1),scalefactor);
		      h3_teta_ptd_tz->Fill(matches.at(j).at(4),matches.at(j).at(1)-matches_noz.at(k).at(1),tzvtx[0],scalefactor);
		      if((matches.at(j).at(1)+matches_noz.at(k).at(1))/2 > 25 && (matches.at(j).at(1)+matches_noz.at(k).at(1))/2 < 30 && abs(matches.at(j).at(1)-matches_noz.at(k).at(1))>5 && abs(tzvtx[0]) >75) cout << "myinfo: " << matches.at(j).at(1) << " " << matches.at(j).at(5) << " " << matches_noz.at(k).at(1) << " " << matches_noz.at(k).at(5) << " " << rzvtx[0] << " " << matches.at(j).at(0) << endl;
		    }
		}
	      //cout << "matches j 0 " << matches.at(j).at(0) << " matches j 1 " << matches.at(j).at(1) << endl;
	      h3_resp_pT_zvtx->Fill(matches.at(j).at(1)/matches.at(j).at(0),matches.at(j).at(0),tzvtx[0],scalefactor);
	      h3_resp_E_zvtx->Fill(matches.at(j).at(3)/matches.at(j).at(2),matches.at(j).at(2),tzvtx[0],scalefactor);
	    }
	  for(int j=0; j<matches_noz.size(); ++j)
	    {
	      h3_resp_pT_zvtx_noz->Fill(matches_noz.at(j).at(1)/matches_noz.at(j).at(0),matches_noz.at(j).at(0),tzvtx[0],scalefactor);
	      h3_resp_E_zvtx_noz->Fill(matches_noz.at(j).at(3)/matches_noz.at(j).at(2),matches_noz.at(j).at(2),tzvtx[0],scalefactor);
	    }
	  //cout << "done with event" << endl;

	}      
      datfile->Close();
    }
  TFile* outfile = TFile::Open(("./multicolhist/hist_zcomp_"+tag+"_"+to_string(rn)+".root").c_str(),"RECREATE");
  h3_resp_pT_zvtx->Write();
  h3_resp_pT_zvtx_noz->Write();
  eff_nz->Write();
  eff_wz->Write();
  noz_recoz_corrET->Write();
  h3_resp_E_zvtx->Write();
  h3_resp_E_zvtx_noz->Write();
  h3_teta_ptd_tz->Write();
  etapt->Write();
  temp2->Write();
  outfile->Write();
  outfile->Close();
  
  return 0;
}
