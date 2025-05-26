int make_tturn(string tag, vector<int> rns, vector<int> nfiles)
{
  //define histograms
  string region = "";
  if(rns[0] < 49216) region = "RegionA";
  else if(rns[0] < 49238) region = "RegionB";
  else if(rns[0] < 49256) region = "RegionC";
  else region = "RegionD";
  const int ntrig = 2;
  int triggers[ntrig] = {18,26};
  TH1D* leadhists[ntrig];
  TH1D* spectra[ntrig];
  TH1D* turn[ntrig];
  TH1D* num[ntrig];
  TH1D* den[ntrig]; 


  double lumi = 0;
  double othervals[5];
  int rnval;
  int rn = rns[0];
  int nfile = nfiles[0];
  ifstream inlumi("/sphenix/user/jocl/projects/analysis/LuminosityCounterGoodRuns/run/list_forplot_2.list");
  while(inlumi >> rnval >> lumi >> othervals[0] >> othervals[1] >> othervals[2] >> othervals[3] >> othervals[4])
    {
      if(rnval == rn) break;
    }
  
  
  for(int j=0; j<ntrig; ++j)
    {
      den[j] = new TH1D(("den_"+to_string(triggers[j])+"_"+region).c_str(),"",10,0,j==0?30:10);
      leadhists[j] = new TH1D(("leadhists_"+to_string(triggers[j])+"_"+region).c_str(),"",50,0,50);
      spectra[j] = new TH1D(("spectra_"+to_string(triggers[j])+"_"+region).c_str(),"",50,0,50);
      //trigturn[j] = new TH1D(("trigturn_"+to_string(triggers[j])+"_"+region).c_str(),"",10,0,30);
      num[j] = new TH1D(("num_"+to_string(triggers[j])+"_"+region).c_str(),"",10,0,j==0?30:10);	 
    }
  
  //define constants
  const int nmaxjet = 100;
  const int nzvtx = 3;
  const int maxtow = 24576+1536*2;
  const int ntowfield = 6;
  const int mbdchan = 64;
  const int mbdside = 2;
  for(int g = 0; g<1; ++g)
    {
  for(int h=0; h<nfile; ++h)
    {
      string filename = "multicoltree/events_"+tag+"_"+to_string(rn)+"_"+to_string(h)+"_0.root";
      cout << "Processing file " << filename << endl;
      TFile* datfile = TFile::Open(filename.c_str());
      //define branch variables
      if(!datfile)
      {
	cout << "no file " << filename << endl;
	continue;
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
      long long unsigned int trigvec[3];
      //define non-branch variables
      double mbdnq, mbdsq, mbdtq;
      double calo_emfrac;//, calo_ohfrac, calo_e;
      double emat[nmaxjet][64];
      double ohat[nmaxjet][64];
      int ncgroe[nmaxjet];
      int ncgroo[nmaxjet];
      double cluster_e[nmaxjet];
      double cluster_eta[nmaxjet];
      int ncluster;
      //get TTree
      TTree* dattree = (TTree*)datfile->Get("tree");
      //set up branches
      dattree->SetBranchAddress("cluster_e",cluster_e);
      dattree->SetBranchAddress("cluster_eta",cluster_eta);
      dattree->SetBranchAddress("ncluster",&ncluster);
      dattree->SetBranchAddress("njet",&njet);
      dattree->SetBranchAddress("hitsgrone",&hitsgrone);
      dattree->SetBranchAddress("isdijet",&isdijet);
      dattree->SetBranchAddress("trigvec",trigvec);
      dattree->SetBranchAddress("jet_e",jet_e);
      dattree->SetBranchAddress("jet_eta",jet_eta);
      dattree->SetBranchAddress("jet_phi",jet_phi);
      dattree->SetBranchAddress("rzvtx",zvtx);
      dattree->SetBranchAddress("towgrone",towgrone);
      //for above: 0 = E, 1 = eta, 2 = phi, 3 = detector level eta, ntrig = calo
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
	  double ETmax = 0;
	  double ETmax_clus = 0;
	  for(int j=0; j<njet; ++j)
	    {
	      double ET = jet_e[j];
	      for(int k=0; k<ntrig; ++k)
		{
		  if((trigvec[2]>>triggers[k])&1 && triggers[k] < 24) spectra[k]->Fill(ET/cosh(jet_eta[j]));
		}

	      if(ET > ETmax) ETmax = ET;
	    }
	  for(int j=0; j<ncluster; ++j)
	    {
	      double ET = cluster_e[j];
	      for(int k=0; k<ntrig; ++k)
		{
		  if((trigvec[2]>>triggers[k])&1 && triggers[k] < 32) spectra[k]->Fill(ET/cosh(cluster_eta[j]));
		}
	      if(ET>ETmax_clus) ETmax_clus = ET;
	    }

	  for(int j=0; j<ntrig; ++j)
	    {
	      if((trigvec[2]>>triggers[j])&1 &&triggers[j] < 24) leadhists[j]->Fill(ETmax/cosh(cluster_eta[j]));
	      else if((trigvec[2]>>triggers[j])&1 &&triggers[j]<32) leadhists[j]->Fill(ETmax_clus/cosh(cluster_eta[j]));
	      if(!((trigvec[2] >> 10) & 1)) continue;
	      if(triggers[j]<24) den[j]->Fill(ETmax);
	      else if(triggers[j]<32) den[j]->Fill(ETmax_clus);
	      if((trigvec[1]>>triggers[j]) & 1 && triggers[j] < 24)
		{
		  num[j]->Fill(ETmax);
		}
	      else if((trigvec[1]>>triggers[j]) &1 && triggers[j] < 32)
		{
		  num[j]->Fill(ETmax_clus);
		}
	    }
	}
      datfile->Close();
    }
    }
  /*
  for(int i=0; i<ntrig; ++i)
    {
      trigturn[i]->Divide(num[i],den,1,1,"B");
    }
  */
  cout << "Finished filling" << endl;
  TFile* outf = TFile::Open(("multicolhist/trigturn_"+tag+"_"+region+(rns.size()==1?"_"+to_string(rns[0]):"")+".root").c_str(),"RECREATE");
  outf->cd();

  cout << "Write hists:" << endl;
  //trigturn->Write();
  for(int i=0; i<ntrig; ++i)
    {
      cout << i << endl;
      num[i]->Write();
      cout << "Wrote num" << endl;
      spectra[i]->Scale(1./lumi);
      spectra[i]->Write();
      cout << "Wrote spectrum" << endl;
      leadhists[i]->Scale(1./lumi);
      leadhists[i]->Write();
      cout << "Wrote lead spectrum" << endl;
      den[i]->Write();
      cout << "Write den" << endl;      
    }

  
  cout << "Write file and close" << endl;
  outf->Write();
  outf->Close();
  return 0;
}
