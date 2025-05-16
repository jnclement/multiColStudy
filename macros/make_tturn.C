int make_tturn(string tag, vector<int> rns, vector<int> nfiles)
{
  //define histograms
  string region = "";
  if(rns[0] < 49216) region = "RegionA";
  else if(rns[0] < 49238) region = "RegionB";
  else if(rns[0] < 49256) region = "RegionC";
  else region = "RegionD";
  TH1D* trigturn = new TH1D(("trigturn_"+region).c_str(),"",10,0,30);
  TH1D* num = new TH1D(("num_"+region).c_str(),"",10,0,30);
  TH1D* den = new TH1D(("den_"+region).c_str(),"",10,0,30);

  //define constants
  const int nmaxjet = 10;
  const int nzvtx = 3;
  const int maxtow = 24576+1536*2;
  const int ntowfield = 6;
  const int mbdchan = 64;
  const int mbdside = 2;
  for(int g = 0; g<rns.size(); ++g)
    {
      int rn = rns.at(g);
      int nfile = nfiles.at(g);
  for(int h=0; h<nfile; ++h)
    {
      string filename = "multicoltree/events_"+tag+"_"+to_string(rn)+"_"+to_string(h)+"_0.root";
      cout << "Processing file " << filename << endl;
      TFile* datfile = TFile::Open(filename.c_str());
      //define branch variables
      {
	cout << "no file " << filename << endl;
	if(!datfile) continue;
      }
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
	  if(!((trigvec >> 10) & 1)) continue;
	  double ETmax = 0;
	  for(int j=0; j<njet; ++j)
	    {
	      double ET = jet_e[j]/cosh(jet_eta[j]);
	      if(ET > ETmax) ETmax = ET;
	    }
	  den->Fill(ETmax);
	  if((trigvec>>18) & 1) num->Fill(ETmax);
	}
      datfile->Close();
    }
  }

  trigturn->Divide(num,den,1,1,"B");
  TFile* outf = TFile::Open(("multicolhist/trigturn_"+tag+"_"+region+".root").c_str(),"RECREATE");
  outf->cd();

  trigturn->Write();
    
  outf->Write();
  outf->Close();
  return 0;
}
