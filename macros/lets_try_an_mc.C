float get_truth_pt(TRandom* rand)
{
  float pt = 0;
  while(pt < 10 || pt > 80)
    {
      pt = rand->Exp(10);
    }
  return pt;
}

float get_zvtx(TRandom* rand)
{
  float z = 1e3;
  while(abs(z)>150)
    {
      z = rand->Gaus(0,50);
    }
  return z;
}

float get_real_eta(TRandom* rand, float pt, float w)
{
  return rand->Gaus(0,w/pt);
}

float get_z0_eta(float z, float eta, float r)
{
  float thetajet = 2*atan(exp(-eta));
  float zjet = r/tan(thetajet);
  float newz = z+zjet;
  float newtheta = atan2(r,newz);
  float neweta = -log(tan(newtheta/2));
  return neweta;
}

float get_truth_E(float pt, float eta)
{
  return pt*cosh(eta);
}

float get_reco_E(TRandom* rand, float truthE, float a, float b, float scale)
{
  float c = b/sqrt(truthE);
  float sig = sqrt(a*a+c*c);
  float m = rand->Gaus(1,sig);
  return truthE*scale*m;
}

float get_reco_pt(float recoE, float eta)
{
  return recoE/cosh(eta);
}

float get_reco_pt_z0(float recoE, float z0eta)
{
  return recoE/cosh(z0eta);
}

vector<TGraphErrors*> get_th2d_mean_tgraph(TH2D* h)
{
  TGraphErrors* means = new TGraphErrors();
  TGraphErrors* jes = new TGraphErrors();
  TGraphErrors* jer = new TGraphErrors();
  TGraphErrors* getrms = new TGraphErrors();
  getrms->SetMarkerStyle(20);
  getrms->SetMarkerColor(kRed);
  getrms->SetLineColor(kRed);
  getrms->SetMarkerSize(1);

  means->SetMarkerStyle(20);
  means->SetMarkerColor(kBlack);
  means->SetLineColor(kBlack);
  means->SetMarkerSize(1.5);
  means->SetLineWidth(2);

  int nx = h->GetNbinsX();
  int npoints = 0;
  //  TCanvas* newcan = new TCanvas("newcan2","",1000,1000);
  for(int i=1; i<=nx; ++i)
    {
      TH1D* projy = h->ProjectionY("_px",i,i);

      TF1* gaus = new TF1("gaus","gaus",projy->GetXaxis()->GetXmin(),projy->GetXaxis()->GetXmax());//projy->GetMean()-projy->GetStdDev(),projy->GetMean()+2*projy->GetStdDev( \
));//                                                                                                                                                
      projy->Fit(gaus,"QRIN");
      //projy->GetYaxis()->SetRangeUser(1e-12,1);
      //projy->Draw("PE");
      //gPad->SetLogy();
      //gPad->SaveAs(("../output/plots/"+std::string(h->GetName())+"_"+to_string(i)+".png").c_str());
      //gPad->SetLogy(0);
      double x = h->GetXaxis()->GetBinCenter(i);
      double mean = gaus->GetParameter(1);
      double err = gaus->GetParameter(2);
      if(mean == 0 || err == 0) continue;
      cout << x << " " <<  mean << " " << err << " " << gaus->GetParError(1) << " " << gaus->GetParError(2) << endl;
      jes->AddPoint(x,mean);
      jes->SetPointError(npoints,0,gaus->GetParError(1));
      jer->AddPoint(x,err/mean);
      jer->SetPointError(npoints,0,gaus->GetParError(2)/mean);
      means->AddPoint(x,mean);
      means->SetPointError(npoints,0,err);
      getrms->AddPoint(x,projy->GetMean());
      getrms->SetPointError(npoints,0,projy->GetStdDev());
      ++npoints;
      //means->AddPointError(x,mean,0,err);                                                                                                          

      //delete gaus;                                                                                                                                 
      //delete projy;                                                                                                                                
    }
  means->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  jes->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  jer->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  getrms->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  return {means,jes,jer,getrms};
}

int lets_try_an_mc()
{
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("","",1000,1000);
  const int ntype = 3;
  float r[ntype] = {125,150,175};
  float w[ntype] = {17.5,20,22.5};
  float a[ntype] = {0.125,0.135,0.145};
  float b[ntype] = {0.6,0.65,0.7};
  float scale[ntype] = {0.675,0.7,0.725};
  TRandom3* rand = new TRandom3();
  for(int g=0; g<ntype; ++g)
    {
      for(int h=0; h<ntype; ++h)
	{
	  for(int i=0; i<ntype; ++i)
	    {
	      for(int j=0; j<ntype; ++j)
		{
		  for(int k=0; k<ntype; ++k)
		    {
		      cout<<g<<h<<i<<j<<k<<endl;
		      TH2D* hist = new TH2D(("hist"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)).c_str(),";(p_{T,noz}^{reco} + p_{T,w/z}^{reco})/2 [GeV];p_{T,noz}^{reco} - p_{T,w/z}^{reco} [GeV]",15,10,85,100,-50,50);
		      for(int l=0; l<1000000; ++l)
			{
			  float pt = get_truth_pt(rand);
			  float z = get_zvtx(rand);
			  float eta;
			  float z0eta = 1e3;
			  while(abs(z0eta) > 0.7)
			    {
			      eta = get_real_eta(rand, pt, w[g]);
			      z0eta = get_z0_eta(z, eta, r[h]);
			    }
			  float recoE = get_reco_E(rand, get_truth_E(pt, eta), a[i], b[j], scale[k]);
			  float recoPt = get_reco_pt(recoE, eta);
			  float recoPtZ0 = get_reco_pt_z0(recoE, z0eta);
			  hist->Fill((recoPt + recoPtZ0)/2,recoPt - recoPtZ0);
			}
		      vector<TGraphErrors*> graphs = get_th2d_mean_tgraph(hist);
		      graphs[0]->SetMarkerStyle(20);
		      graphs[0]->SetLineWidth(2);
		      gPad->SetLogz();
		      hist->Draw("COLZ");
		      graphs[0]->Draw("SAME PE");
		      gPad->SaveAs(("../output/plots/mctest"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)+".png").c_str());
		      //goto end;
		    }
		}
	    }
	}
    }
 end:
  return 0;
}
