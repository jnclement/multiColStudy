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
  return a[0] < b[0];
}

vector<vector<double>> make_jet_vector(int njet, double* jet_pt, double* jet_eta, double* jet_phi)
{
  vector<vector<double>> jet_vector = {};
  for(int i=0; i<njet; ++i)
    {
      vector<double> jet = {jet_pt[i],jet_eta[i],jet_phi[i]};
      jet_vector.push_back(jet);
    }
  std::sort(jet_vector.begin(),jet_vector.end(),compfirst);
  return jet_vector;
}

vector<vector<double>> truth_match(vector<vector<double>> truth_jets, vector<vector<double>> reco_jets)
{
  vector<vector<double>> matches = {};
  for(int i=0; i<truth_jets.size(); ++i)
    {
      int match_index = -1;
      double max_pt = 0;
      vector<double> match = {};
      for(int j=0; j<reco_jets.size(); ++j)
	{
	  double dPhi = abs(truth_jets.at(i).at(2) - reco_jets.at(j).at(2));
	  if(dPhi < 0.3 && reco_jets.at(j).at(0) > max_pt)
	    {
	      max_pt = reco_jets.at(j).at(0);
	      match_index = j;
	    }
	}
      match.push_back(truth_jets.at(i).at(0));
      match.push_back(reco_jets.at(match_index).at(0));
      matches.push_back(match);
      reco_jets.erase(reco_jets.begin()+match_index);
    }
  return matches;
}

int comp_zvtx(string tag, int rn)
{
  
  const int nmaxjet = 100;
  const int nzvtx = 3;
  TH3D* h3_resp_pT_zvtx = new TH3D("h3_resp_pT_zvtx",";p_{T}^{reco}/p_{T}^{truth};p_{T}^{truth} [GeV];z_{vtx} [cm]",100,0,1,100,0,100,300,-150,150);
  TH3D* h3_resp_pT_zvtx_noz = new TH3D("h3_resp_pT_zvtx_noz",";p_{T}^{reco}/p_{T}^{truth};p_{T}^{truth} [GeV];z_{vtx} [cm]",100,0,1,100,0,100,300,-150,150);
  for(int h=rn*100; h<rn*100+100; ++h)
    {
      string filename = "/sphenix/tg/tg01/jets/jocl/multiCol/"+to_string(rn)+"/events_"+tag+"_"+to_string(h)+"_0.root";
      cout << "Processing file " << filename << endl;
      TFile* datfile = TFile::Open(filename.c_str());
      
      if(!datfile) continue;
      TTree* tree = (TTree*)datfile->Get("tree");
      if(!tree) continue;
      int njet, tnjet, njet_noz;
      double jet_pt[nmaxjet], jet_eta[nmaxjet], jet_pt_noz[nmaxjet], jet_eta_noz[nmaxjet], tjet_pt[nmaxjet], tjet_eta[nmaxjet], jet_phi[nmaxjet], tjet_phi[nmaxjet], jet_phi_noz[nmaxjet], tzvtx[3];
      
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
      

      
      for(int i=0; i<tree->GetEntries(); ++i)
	{
	  tree->GetEntry(i);
	  vector<vector<double>> recojets = make_jet_vector(njet, jet_pt, jet_eta, jet_phi);
	  vector<vector<double>> truthjet = make_jet_vector(tnjet, tjet_pt, tjet_eta, tjet_phi);
	  vector<vector<double>> reco_noz = make_jet_vector(njet_noz, jet_pt_noz, jet_eta_noz, jet_phi_noz);
	  
	  vector<vector<double>> matches = truth_match(truthjet, recojets);
	  vector<vector<double>> matches_noz = truth_match(truthjet, reco_noz);
	  
	  for(int j=0; j<matches.size(); ++j)
	    {
	      h3_resp_pT_zvtx->Fill(matches.at(j).at(1)/matches.at(j).at(0),matches.at(j).at(0),tzvtx[0]);
	      h3_resp_pT_zvtx_noz->Fill(matches_noz.at(j).at(1)/matches_noz.at(j).at(0),matches_noz.at(j).at(0),tzvtx[0]);
	    }
	}      
      datfile->Close();
    }
  TFile* outfile = TFile::Open(("./multicolhist/hist_zcomp_"+tag+"_"+to_string(rn)+".root").c_str());
  h3_resp_pT_zvtx->Write();
  h3_resp_pT_zvtx_noz->Write();
  outfile->Write();
  outfile->Close();
  
  return 0;
}
