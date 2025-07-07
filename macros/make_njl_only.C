double find_etsub_over_etlead(double* jets, int njets)
{
  double lead = 0;
  double sublead = 0;

  for(int i=0; i<njets; ++i)
    {
      if(jets[i] > lead)
	{
	  sublead = lead;
	  lead = jets[i];
	}
      else if(jets[i] > sublead)
	{
	  sublead = jets[i];
	}
    }
  return sublead/lead;
}

int make_njl_only(string tag, vector<int> rns, vector<int> nfiles, int triggerbit = 22, int dodijetcut = 1)
{
  TH2D* njl_rn_jetet = new TH2D("njl_rn_jetet",";Run Number;E_{T}^{jet} [GeV]",53900-47200,47200,53900,100,0,100);
  const int nmaxjet = 100;
  const int nzvtx = 3;
  const int maxtow = 24576+1536*2;
  const int ntowfield = 6;
  const int mbdchan = 64;
  const int mbdside = 2;
  //start processing

  int whichlumi = 2;
  if(triggerbit==18) whichlumi = 0;
  else if(triggerbit==22) whichlumi = 2;
  double lumi[6] = {0};
  int rnval;
  int rn = rns[0];
  int nfile = nfiles[0];
  ifstream inlumi("/sphenix/user/jocl/projects/LumiList/lumi_fornjl.list");
  while(inlumi >> rnval >> lumi[0] >> lumi[1] >> lumi[2] >> lumi[3] >> lumi[4] >> lumi[5])
    {
      if(rnval == rn) break;
    }
  cout << "RN/LUMI: " << rnval << " " << lumi[whichlumi] << endl;
  for(int h=0; h<nfile; ++h)
    {
      string filename = "multicoltree/events_"+tag+"_"+to_string(rn)+"_"+to_string(h)+"_0.root";
      cout << "Processing file " << filename << endl;
      TFile* datfile = TFile::Open(filename.c_str());
      //define branch variables
      if(!datfile) continue;
      int njet, isdijet;
      double dphilead;
      double jet_e[nmaxjet];
      long long unsigned int trigvec[3];
      double zvtx[3];
      TTree* dattree = (TTree*)datfile->Get("tree");
      if(!dattree) continue;
      //set up branches
      dattree->SetBranchAddress("njet_noz",&njet);
      dattree->SetBranchAddress("isdijet",&isdijet);
      dattree->SetBranchAddress("trigvec",trigvec);
      dattree->SetBranchAddress("jet_pt_noz",jet_e);
      dattree->SetBranchAddress("dphilead",&dphilead);
      dattree->SetBranchAddress("rzvtx",zvtx);
      
      for(int i=0; i<dattree->GetEntries(); ++i)
	{
	  cout << zvtx[0] << " " << zvtx[1] << endl;
	  //if(std::isnan(zvtx[0])) cout << "nan zvtx! hooray! " << zvtx[0] << endl;
	  dattree->GetEntry(i);
	  if(!((trigvec[2] >> triggerbit) & 1)) continue;	  
	  if(dphilead < 3*M_PI/4 || isdijet == 0 ||  find_etsub_over_etlead(jet_e,njet) < 0.3) continue;
	  for(int j=0; j<njet; ++j)
	    {
	      if(lumi[whichlumi] != 0 && triggerbit == 22) njl_rn_jetet->Fill(rn,jet_e[j],1./lumi[whichlumi]);
	    }
	}
      datfile->Close();
    }
  TFile* outf = TFile::Open(("multicolhist/hists_njlonly_"+tag+"_"+(rns.size()==1?"_"+to_string(rns[0])+"_":"_")+to_string(triggerbit)+".root").c_str(),"RECREATE");
  outf->cd();

  outf->cd();

  njl_rn_jetet->Write();

  outf->Close();

  return 0;
}
