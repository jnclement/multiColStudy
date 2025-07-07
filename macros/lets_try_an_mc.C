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

float get_truth_pt(TRandom* rand)
{
  float pt = 0;
  while(pt < 10 || pt > 80)
    {
      pt = rand->Exp(30);
    }
  return pt;
}

float get_zvtx(TRandom* rand, float max, float getgt)
{
  float z = getgt?0:1e3;
  while((getgt?abs(z)<max:abs(z)>max) || abs(z) > 150)
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
  TCanvas* newcan = new TCanvas("newcan2","",1000,1000);
  for(int i=1; i<=nx; ++i)
    {
      TH1D* projy = h->ProjectionY("_px",i,i);

      TF1* gaus = new TF1("gaus","gaus",projy->GetXaxis()->GetXmin(),projy->GetXaxis()->GetXmax());//projy->GetMean()-projy->GetStdDev(),projy->GetMean()+2*projy->GetStdDev( \
));//                                                                                                                                                
      projy->Fit(gaus,"QRIN");
      projy->GetYaxis()->SetRangeUser(1e-6,1);
      projy->Draw("PE");
      gPad->SetLogy();
      gPad->SaveAs(("../output/plots/"+std::string(h->GetName())+"_"+to_string(i)+".png").c_str());
      gPad->SetLogy(0);
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

      delete gaus;                                                                                                                                 
      delete projy;                                                                                                                                
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
  can->cd();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  const int ntype = 3;
  int nend = 2;
  int nstart = 1;
  float r[ntype] = {125,150,175};
  float w[ntype] = {17.5,25,22.5};
  float a[ntype] = {0.125,0.135,0.145};
  float b[ntype] = {0.6,0.65,0.7};
  float scale[ntype] = {0.675,0.7,0.725};
  TRandom3* rand = new TRandom3();
  for(int g=nstart; g<nend; ++g)
    {
      for(int h=nstart; h<nend; ++h)
	{
	  for(int i=nstart; i<nend; ++i)
	    {
	      for(int j=nstart; j<nend; ++j)
		{
		  for(int k=nstart; k<nend; ++k)
		    {
		      cout<<g<<h<<i<<j<<k<<endl;
		      TH2D* hist = new TH2D(("hist"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)).c_str(),";(p_{T,noz}^{reco} + p_{T,w/z}^{reco})/2 [GeV];p_{T,w/z}^{reco} - p_{T,noz}^{reco} [GeV]",15,10,85,100,-50,50);
		      TH2D* etapt = new TH2D(("etapt"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)).c_str(),";p_{T}^{truth};Jet #eta",15,10,85,30,-1.5,1.5);
		      //TH2D* etaz0pt = new TH2D(("etaz0pt"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)).c_str(),";p_{T}^{truth};Jet #eta_{z=0}",15,10,85,30,-1.5,1.5);
		      for(int l=0; l<10000000; ++l)
			{
			  if(l%10000 == 0) cout << l << endl;
			  float pt = get_truth_pt(rand);
			  float z = get_zvtx(rand,150,0);
			  float eta = 1e3;
			  float z0eta = 1e3;
			  while(check_bad_jet_eta(eta,z,0.4))
			    {
			      eta = get_real_eta(rand, pt, w[g]);
			      z0eta = get_z0_eta(z, eta, 93);
			      
			    }

			  etapt->Fill(pt,eta);
			  //etaz0pt->Fill(pt,z0eta);
			  float recoE = get_reco_E(rand, get_truth_E(pt, eta), a[i], b[j], scale[k]);
			  float recoPt = get_reco_pt(recoE, eta);
			  float recoPtZ0 = get_reco_pt_z0(recoE, z0eta);
			  //if((recoPt + recoPtZ0)/2 < 30 && (recoPt + recoPtZ0)/2 > 25 && abs(recoPt - recoPtZ0) > 5 && abs(z) > 75) cout << recoPt << " " << eta << " " << recoPtZ0 << " " << z0eta << " "  << z << endl;
			  hist->Fill((recoPt + recoPtZ0)/2,recoPt - recoPtZ0);
			}
		      hist->Scale(1./hist->Integral());
		      etapt->Scale(1./hist->Integral());
		      vector<TGraphErrors*> graphs = get_th2d_mean_tgraph(hist);
		      vector<TGraphErrors*> graphseta = get_th2d_mean_tgraph(etapt);
		      //vector<TGraphErrors*> graphsz0eta = get_th2d_mean_tgraph(etaz0pt);
		      
		      graphs[0]->SetMarkerStyle(20);
		      graphs[0]->SetLineWidth(2);
		      gPad->SetLogz();
		      hist->Draw("COLZ");
		      graphs[0]->Draw("SAME PE");
		      gPad->SaveAs(("../output/plots/mctest"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)+".png").c_str());

		      etapt->Draw("COLZ");
		      graphseta[0]->SetMarkerStyle(20);
		      graphseta[0]->SetLineWidth(2);
		      graphseta[0]->Draw("SAME PE");
		      gPad->SaveAs(("../output/plots/etapt"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)+".png").c_str());
		      /*
		      etaz0pt->Draw("COLZ");
		      graphsz0eta[0]->Draw("SAME PE");
		      gPad->SaveAs(("../output/plots/etaz0pt"+to_string(g)+to_string(h)+to_string(i)+to_string(j)+to_string(k)+".png").c_str());
		      */
		      //goto end;
		    }
		}
	    }
	}
    }
 end:
  return 0;
}
