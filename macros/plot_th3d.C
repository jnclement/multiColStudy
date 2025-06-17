#include <dlUtility.h>
#include <TGraphErrors.h>

int set_axis_range_based_on_axis(TH3D* h, int axis, int from, int to)
{
  if(axis==0) h->GetXaxis()->SetRange(from,to);
  else if(axis==1) h->GetYaxis()->SetRange(from,to);
  else if(axis==2) h->GetZaxis()->SetRange(from,to);
  else return 1;
  return 0;
}

int get_bin_number_for_axis(TH3D* h, int axis, double value, int getlo)
{
  int bin = -1;
  if(axis==0) bin = h->GetXaxis()->FindBin(value);
  else if(axis==1) bin = h->GetYaxis()->FindBin(value);
  else if(axis==1) bin = h->GetZaxis()->FindBin(value);
  else return bin;
  bin += (getlo?0:1);
  return bin;
}
vector<TObject*> make_projections(TH3D* h, int axis, double from, double to)
{

  string projdir;
  vector<TObject*> outhists = {};

  
  if(axis==0) projdir = "yzo";
  else if(axis==1) projdir = "xzo";
  else if(axis==2) projdir = "xyo";
  else return {};

  TH2D* proj1;
  TH2D* proj2;
  
  if(to > from)
    {
      set_axis_range_based_on_axis(h, axis, get_bin_number_for_axis(h, axis, from, 1), get_bin_number_for_axis(h, axis, to, 0));
      proj1 = (TH2D*)(h->Project3D(projdir.c_str()));
    }
  else
    {
      set_axis_range_based_on_axis(h, axis, 0, get_bin_number_for_axis(h, axis, to, 0));
      proj1 = (TH2D*)(h->Project3D(projdir.c_str()));

      set_axis_range_based_on_axis(h, axis, get_bin_number_for_axis(h, axis, from, 1), 9999999);
      proj2 = (TH2D*)(h->Project3D(projdir.c_str()));
      proj1->Add(proj2);
    }

  proj1->GetZaxis()->SetTitle("Counts");

  
  outhists.push_back(proj1);
  outhists.push_back(proj1->ProjectionX());
  outhists.push_back(proj1->ProjectionY());
  
  return outhists;
}

TGraphErrors* get_th2d_mean_tgraph(TH2D* h)
{
  TGraphErrors* means = new TGraphErrors();
  means->SetMarkerStyle(20);
  means->SetMarkerColor(kBlack);
  means->SetLineColor(kBlack);
  means->SetMarkerSize(1.5);
  means->SetLineWidth(2);
  
  int nx = h->GetNbinsX();

  for(int i=1; i<=nx; ++i)
    {
      TH1D* projy = h->ProjectionY("_py",i,i);
      if(projy->GetEntries() < 5) continue;
      TF1* gaus = new TF1("gaus","gaus",projy->GetXaxis()->GetXmin(),projy->GetXaxis()->GetXmax());
      projy->Fit(gaus,"QNR");
      double x = h->GetXaxis()->GetBinCenter(i);
      double mean = gaus->GetParameter(1);
      double err = gaus->GetParameter(2);
      means->AddPoint(x,mean);
      means->SetPointError(0,err);
      //means->AddPointError(x,mean,0,err);

      delete gaus;
      delete projy;
    }
  return means;
}

int format_th1d(TH1D* hist, int color)
{
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);

  hist->SetMarkerSize(1.5);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(20);

  return 0;
}

int format_th2d(TH2D* hist)
{
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetZaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetTitleSize(0.03);
  return 0;
}

int alltext(int hasz, string etcuttext = "", string zcuttext = "")
{
  sphenixtext(0.85,0.96);
  sqrt_s_text(0.85,0.92);
  antikt_text(0.4,0.6,0.92);
  //dijet_cut_text(0.6,0.96);
  drawText(etcuttext.c_str(),0.4,0.96);
  drawText(zcuttext.c_str(),0.4,0.92);
  string hasztext = "";
  if(hasz == 0) hasztext = "z_{vtx} assumed 0";
  else if(hasz == 1) hasztext = "Reco z_{vtx} must exist";
  drawText(hasztext.c_str(),0.6,0.96);
  return 0;
}

TGraphErrors* getRatioGraphBinomial(TGraphErrors* gNum, TGraphErrors* gDen)
{
  int n = gNum->GetN();
  TGraphErrors* gRatio = new TGraphErrors();
  
  for (int i = 0; i < n; ++i) {
    double x1, y1, x2, y2;
    gNum->GetPoint(i, x1, y1);
    gDen->GetPoint(i, x2, y2);
    
    if (x1 != x2) {
      std::cerr << "Warning: X points do not match at i = " << i << std::endl;
      continue;
    }
    
    if (y2 == 0) continue;  // avoid divide by zero
    
    double ratio = y1 / y2;
    double binomError = sqrt(ratio * (1 - ratio) / y2);
    
    gRatio->SetPoint(i, x1, ratio);
    gRatio->SetPointError(i, 0, binomError);  // x error = 0
  }
  
  gRatio->SetMarkerStyle(20);
  gRatio->SetMarkerColor(kBlack);
  gRatio->SetLineColor(kBlack);
  gRatio->SetTitle("Ratio with Binomial Errors");
  
  return gRatio;
}

int overlay_w_ratio_tgraph(TGraphErrors* wz, TGraphErrors* nz, TCanvas* can, TLegend* leg, string etcuttext, string zcuttext, string histtype, string xy)
{
  can->cd(1);
  TGraphErrors* ratio = getRatioGraphBinomial(wz,nz);
  ratio->GetYaxis()->SetTitle("Ratio Reco z / Zero z");
  //ratio->Divide(wz,nz,1,1,"B");
  
  double ymax = wz->GetMaximum();
  double nzmax = nz->GetMaximum();
  if(nzmax > ymax) ymax = nzmax;
  ymax *= 1.5; 
  wz->GetYaxis()->SetRangeUser(0,ymax);
  nz->GetYaxis()->SetRangeUser(0,ymax);

  wz->Draw("PE");
  nz->Draw("SAME PE");
  can->cd(2);
  ratio->Draw("PE");
  can->cd(0);
  alltext(-1, etcuttext, zcuttext);
  can->cd(1);
  leg->Draw();
  TLine* oneline = new TLine(wz->GetXaxis()->GetXmin(),1,nz->GetXaxis()->GetXmax(),1);
  can->cd(2);
  oneline->Draw();

  can->SaveAs(("../output/plots/"+histtype + "_"+xy+".png").c_str());

  delete oneline;
  delete ratio;
  
  return 0;
}

int overlay_w_ratio_th1d(TH1D* wz, TH1D* nz, TCanvas* can, TLegend* leg, string etcuttext, string zcuttext, string histtype, string xy)
{
  can->cd(1);
  TH1D* ratio = (TH1D*)wz->Clone();
  ratio->Divide(wz,nz,1,1,"B");
  ratio->GetYaxis()->SetTitle("Ratio Reco z / Zero z");
  
  double ymax = wz->GetMaximum();
  double nzmax = nz->GetMaximum();
  if(nzmax > ymax) ymax = nzmax;
  ymax *= 1.5; 
  wz->GetYaxis()->SetRangeUser(0,ymax);
  nz->GetYaxis()->SetRangeUser(0,ymax);

  wz->Draw("PE");
  nz->Draw("SAME PE");
  can->cd(2);
  ratio->Draw("PE");
  can->cd(0);
  alltext(-1, etcuttext, zcuttext);
  can->cd(1);
  leg->Draw();
  TLine* oneline = new TLine(wz->GetXaxis()->GetXmin(),1,nz->GetXaxis()->GetXmax(),1);
  can->cd(2);
  oneline->Draw();

  can->SaveAs(("../output/plots/"+histtype + "_"+xy+".png").c_str());

  delete oneline;

  return 0;
}

int draw_all(string histtype, vector<TObject*> wz, vector<TObject*> nz, string etcuttext = "", string zcuttext = "")
{

  TH2D* wz2 = (TH2D*)wz.at(0);
  TH1D* wzx = (TH1D*)wz.at(1);
  TH1D* wzy = (TH1D*)wz.at(2);

  TH2D* nz2 = (TH2D*)nz.at(0);
  TH1D* nzx = (TH1D*)nz.at(1);
  TH1D* nzy = (TH1D*)nz.at(2);

  format_th2d(wz2);
  format_th2d(nz2);
  format_th1d(wzx,kBlack);
  format_th1d(wzy,kBlack);
  format_th1d(nzx,kRed);
  format_th1d(nzy,kRed);

  TGraphErrors* wzg = get_th2d_mean_tgraph(wz2);
  TGraphErrors* nzg = get_th2d_mean_tgraph(nz2);

  TCanvas* can = new TCanvas("","",1000,1000);
  can->cd();
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.2);

  wz2->Draw("COLZ");
  wzg->Draw("SAME PE");
  alltext(1,etcuttext,zcuttext);
  can->SaveAs(("../output/plots/"+histtype + "_wz.png").c_str());

  nz2->Draw("COLZ");
  nzg->Draw("SAME PE");
  alltext(0,etcuttext,zcuttext);
  can->SaveAs(("../output/plots/"+histtype + "_nz.png").c_str());

  delete can;
  can = new TCanvas("","",1000,1500);
  ratioPanelCanvas(can,0.3);

  TLegend* leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->AddEntry(wzx,"With Reco z_{vtx}","p");
  leg->AddEntry(nzx,"z_{vtx} Assumed 0","p");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  overlay_w_ratio_th1d(wzx, nzx, can, leg, etcuttext, zcuttext, histtype, "x");
  overlay_w_ratio_th1d(wzy, nzy, can, leg, etcuttext, zcuttext, histtype, "y");
  overlay_w_ratio_tgraph(wzg, nzg, can, leg, etcuttext, zcuttext, histtype, "mean");

  delete leg;
  delete can;
  delete wzg;
  delete nzg;
  
  return 0;
}

int get_and_draw(TH3D* hw, TH3D* hn, int axis, double from, double to, string etcuttext = "", string zcuttext = "")
{

  string histtype = "proj";
  if(axis == 0) histtype += "x";
  else if(axis==1) histtype += "y";
  else if(axis==2) histtype += "z";
  else return 1;

  histtype += "_";
  histtype += to_string(from);
  histtype += "_";
  histtype += to_string(to);
  
  vector<TObject*> wz = make_projections(hw, axis, from, to);
  vector<TObject*> nz = make_projections(hn, axis, from, to);
  
  draw_all(histtype,wz,nz,etcuttext,zcuttext);

  return 0;
}

int plot_th3d(string filename)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile* file = TFile::Open(filename.c_str());
  if(!file) return 1;
  TH3D* hw = (TH3D*)file->Get("h3_resp_pT_zvtx");
  TH3D* hn = (TH3D*)file->Get("h3_resp_pT_zvtx_noz");
  if(!hw || !hn) return 2;
  
  get_and_draw(hw, hn, 1, 20, 100, "E_{T}^{jet} > 20 GeV","All z_{vtx}");
  get_and_draw(hw, hn, 2, -30, 30, "E_{T}^{jet} > 15 GeV","|z_{vtx}|<30 cm");
  get_and_draw(hw, hn, 2, 60, -60, "E_{T}^{jet} > 15 GeV","|z_{vtx}|>60 cm");
  
  
  return 0;
}
