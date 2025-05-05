#include <dlUtility.h>
const int nhistplot = 4;
string region[nhistplot] = {"RegionA", "RegionB","RegionC","RegionD"};
int format_th1d(TH1D* hist, string xtitle, string ytitle, int n = 1)
{
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);


  hist->SetMarkerSize(1.5);
  if(n==0)
    {
      hist->SetMarkerColor(kBlack);
      hist->SetLineColor(kBlack);
      hist->SetLineWidth(2);
      hist->SetMarkerStyle(20);
    }
  else if(n==1)
    {
      hist->SetMarkerColor(kCyan+2);
      hist->SetLineColor(kCyan+2);
      hist->SetLineWidth(2);
      hist->SetMarkerStyle(21);
    }
  else if(n==2)
    {
      hist->SetMarkerColor(kMagenta+2);
      hist->SetLineColor(kMagenta+2);
      hist->SetLineWidth(2);
      hist->SetMarkerStyle(53);
    }
  else
    {
      hist->SetMarkerColor(kYellow+2);
      hist->SetLineColor(kYellow+2);
      hist->SetLineWidth(2);
      hist->SetMarkerStyle(54);
    }

  return 0;
}

int draw_overlay_with_ratio_th1d(TH1D** hists, string histbasename, TCanvas* can, int dijetcut)
{
  can->cd();
  can->SetRightMargin(0.2);
  can->SetBottomMargin(0.15);
  can->SetLeftMargin(0.15);
  can->SetTopMargin(0.15);
  TH1D* ratio[nhistplot];
  static int has_divided = 0;
  if(!has_divided)
    {
      ratioPanelCanvas(can,0.3);
      has_divided = 1;
    }
  can->cd(1);
  float ymax = 0;
  TLegend* leg = new TLegend(0.4,0.65,0.6,0.85);
  for(int i=0; i<nhistplot; ++i)
    {
      if(hists[i]->GetMaximum() > ymax) ymax = hists[i]->GetMaximum();
      leg->AddEntry(hists[i],region[i].c_str(),"p");
      ratio[i] = new TH1D((histbasename+"_"+region[0]+"_ratio_"+region[i]).c_str(),"",hists[0]->GetNbinsX(),hists[0]->GetXaxis()->GetXmin(),hists[0]->GetXaxis()->GetXmax());
      if(i>0) ratio[i]->Divide(hists[0],hists[i]);
      ratio[i]->GetYaxis()->SetRangeUser(0,2);
      format_th1d(ratio[i],hists[i]->GetXaxis()->GetTitle(),"Ratio of Region A to Region B,C,D",i);
      ratio[i]->GetXaxis()->SetTitleSize(ratio[i]->GetXaxis()->GetTitleSize()/0.5);
      ratio[i]->GetYaxis()->SetTitleSize(ratio[i]->GetYaxis()->GetTitleSize()/0.5);
      ratio[i]->GetXaxis()->SetLabelSize(ratio[i]->GetXaxis()->GetLabelSize()/0.5);
      ratio[i]->GetYaxis()->SetLabelSize(ratio[i]->GetYaxis()->GetLabelSize()/0.5);
    }
  
  hists[0]->GetYaxis()->SetRangeUser(0,1.5*ymax);
  for(int i=0; i<nhistplot; ++i)
    {
      can->cd(1);
      hists[i]->GetXaxis()->SetLabelSize(0);
      hists[i]->GetXaxis()->SetTitleSize(0);
      hists[i]->GetYaxis()->SetRangeUser(0,1.5*ymax);
      hists[i]->Draw(i==0?"PE":"SAME PE");
      can->cd(2);
      if(i>0) ratio[i]->Draw(i==1?"PE":"SAME PE");
    }

  can->cd(1);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  float minet = 15;
  if(std::string(hists[0]->GetName()).find("gr15") != std::string::npos) minet = 20;
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  TLine* oneline = new TLine(hists[0]->GetXaxis()->GetXmin(),1,hists[0]->GetXaxis()->GetXmax(),1);
  can->cd(2);
  oneline->Draw();
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hists[0]->GetName())+"_overlay_"+(dijetcut?"dc":"nc")+".pdf").c_str());

  ymax = 0;
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i]->Scale(1./hists[i]->Integral());
      hists[i]->GetYaxis()->SetTitle("Counts / N_{bit18}^{scaled} (Normalized to #Sigma(bins)=1)");
      if(hists[i]->GetMaximum() > ymax) ymax = hists[i]->GetMaximum();
    }
  
  hists[0]->GetYaxis()->SetRangeUser(0,1.5*ymax);
  for(int i=0; i<nhistplot; ++i)
    {
      can->cd(1);
      hists[i]->GetXaxis()->SetLabelSize(0);
      hists[i]->GetXaxis()->SetTitleSize(0);
      hists[i]->GetYaxis()->SetRangeUser(0,1.5*ymax);
      if(i>0) ratio[i]->Divide(hists[0],hists[i]);
      hists[i]->Draw(i==0?"PE":"SAME PE");
      can->cd(2);
      ratio[i]->GetYaxis()->SetRangeUser(0,2);
      if(i>0) ratio[i]->Draw(i==1?"PE":"SAME PE");
    }

  can->cd(1);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->cd(2);
  oneline->Draw();

  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hists[0]->GetName())+"_overlay_"+(dijetcut?"dc":"nc")+"_normed.pdf").c_str());


  
  delete leg;
  for(int i=0; i<nhistplot; ++i)
    {
      delete ratio[i];
    }
  return 0;
}

int draw_th1d(TH1D* hist, TCanvas* can, int dijetcut)
{
  can->cd();
  can->SetRightMargin(0.2);
  can->SetBottomMargin(0.15);
  can->SetLeftMargin(0.15);
  can->SetTopMargin(0.15);

  float minet = 15;
  if(std::string(hist->GetName()).find("gr15") != std::string::npos) minet = 20;
  hist->Draw("PE");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+(dijetcut?"dc":"nc")+".pdf").c_str());

  gPad->SetLogy();
  hist->Draw("PE");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+(dijetcut?"dc":"nc")+"_log.pdf").c_str());
  gPad->SetLogy(0);
  return 0;
}

int draw_ratio_th1d(TH1D** hists, string histbasename, string* region, TCanvas* can, int dijetcut)
{
  int nxbin = hists[0]->GetNbinsX();
  double xmax = hists[0]->GetXaxis()->GetXmax();
  double xmin = hists[0]->GetXaxis()->GetXmin();
  TH1D* ratio = new TH1D((histbasename+"_"+region[0]+"_to_"+region[1]).c_str(),"",nxbin,xmin,xmax);

  ratio->Divide(hists[0],hists[1]);

  format_th1d(ratio, hists[0]->GetXaxis()->GetTitle(), ("Ratio of "+region[0]+" to "+region[1]).c_str());
  ratio->GetYaxis()->SetRangeUser(0,2);
  draw_th1d(ratio, can, dijetcut);
  delete ratio;
  return 0;
}

int get_and_draw_th1d(string histbasename, string* region, TFile* histfile, string xtitle, string ytitle, TCanvas* can, TCanvas* ratcan, int dijetcut)
{
  TH1D* hists[nhistplot];
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i] = (TH1D*)histfile->Get((histbasename+"_"+region[i]).c_str());
      if(std::string(hists[i]->GetName()).find("zhist") != std::string::npos) hists[i]->Rebin(10);
      format_th1d(hists[i], xtitle, ytitle,i);
      //draw_th1d(hists[i], can, dijetcut);
    }
  draw_overlay_with_ratio_th1d(hists,histbasename,ratcan,dijetcut);
  //draw_ratio_th1d(hists, histbasename, region, can, dijetcut);
  
  return 0;
}
int format_th2d(TH2D* hist, string xtitle, string ytitle, string ztitle)
{
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetZaxis()->SetTitle(ztitle.c_str());
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetZaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetTitleSize(0.03);

  return 0;
}

int draw_th2d(TH2D* hist, TCanvas* can, int dijetcut)
{
  can->cd();
  can->SetRightMargin(0.2);
  can->SetBottomMargin(0.15);
  can->SetLeftMargin(0.15);
  can->SetTopMargin(0.15);
  
  hist->Draw("COLZ");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  float minet = 15;
  if(std::string(hist->GetName()).find("gr15") != std::string::npos) minet = 20;
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+(dijetcut?"dc":"nc")+".pdf").c_str());

  gPad->SetLogz();
  hist->Draw("COLZ");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+(dijetcut?"dc":"nc")+"_log.pdf").c_str());
  gPad->SetLogz(0);
  return 0;
}

int draw_ratio_th2d(TH2D** hists, string histbasename, string* region, TCanvas* can, int dijetcut)
{
  int nxbin = hists[0]->GetNbinsX();
  int nybin = hists[0]->GetNbinsY();
  double xmax = hists[0]->GetXaxis()->GetXmax();
  double ymax = hists[0]->GetYaxis()->GetXmax();
  double xmin = hists[0]->GetXaxis()->GetXmin();
  double ymin = hists[0]->GetYaxis()->GetXmin();
  TH2D* ratio = new TH2D((histbasename+"_"+region[0]+"_to_"+region[1]).c_str(),"",nxbin,xmin,xmax,nybin,ymin,ymax);

  ratio->Divide(hists[0],hists[1]);
  format_th2d(ratio, hists[0]->GetXaxis()->GetTitle(), hists[0]->GetYaxis()->GetTitle(), ("Ratio of "+region[0]+" to "+region[1]).c_str());
  ratio->GetZaxis()->SetRangeUser(0,2);
  draw_th2d(ratio, can, dijetcut);
  delete ratio;
  return 0;
}

int get_and_draw_th2d(string histbasename, string* region, TFile* histfile, string xtitle, string ytitle, string ztitle, TCanvas* can, TCanvas* ratcan, int dijetcut)
{
  TH2D* hists[nhistplot];
  TH1D* projx[nhistplot];
  TH1D* projy[nhistplot];
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i] = (TH2D*)histfile->Get((histbasename+"_"+region[i]).c_str());
      if(std::string(hists[i]->GetName()).find("frcem") != std::string::npos) hists[i]->Rebin2D(5,5);
      format_th2d(hists[i], xtitle, ytitle, ztitle);
      //draw_th2d(hists[i], can,dijetcut);
      projx[i] = hists[i]->ProjectionX();
      projy[i] = hists[i]->ProjectionY();

      format_th1d(projx[i],hists[i]->GetXaxis()->GetTitle(),hists[i]->GetZaxis()->GetTitle(),i);
      format_th1d(projy[i],hists[i]->GetYaxis()->GetTitle(),hists[i]->GetZaxis()->GetTitle(),i);

      //draw_th1d(projx[i],can,dijetcut);
      //draw_th1d(projy[i],can,dijetcut);
    }

  //draw_ratio_th1d(projx, histbasename+"_projx", region, can, dijetcut);
  //draw_ratio_th1d(projy, histbasename+"_projy", region, can, dijetcut);
  //draw_ratio_th2d(hists, histbasename, region, can, dijetcut);
  draw_overlay_with_ratio_th1d(projx, histbasename+"_projx",ratcan,dijetcut);
  draw_overlay_with_ratio_th1d(projy, histbasename+"_projy",ratcan,dijetcut);
  return 0;
}

int plot()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile* histfile;
  histfile = TFile::Open("../output/hists/hadded.root");
  int dijetcut = 1;
  const int nhist = 24;

  TCanvas* can = new TCanvas("can","",1000,1000);
  TCanvas* ratcan = new TCanvas("ratcan","",1000,1500);
  string histnames[nhist] = {"jet_eta_phi","jet_eta_phi_gr15","tow_eta_phi","tow_eta_phi_gr5","tow_deteta_phi","tow_deteta_phi_gr5","jet_E_eta","tow_E_eta","tow_E_deteta","jet_E_phi","tow_E_phi","jet_ET_eta","tow_ET_eta","tow_ET_deteta","jet_ET_phi","tow_ET_phi","jet_E_frcem","jet_ET_dphi","jet_eta_frcem","jet_eta_frcem_gr15","jet_phi_frcem","jet_phi_frcem_gr15","frcem_frcoh","frcem_frcoh_gr15"};

  string xtitles[nhist] = {"Jet #eta","Jet #eta","Tower #eta","Tower #eta","Tower Detector #eta","Tower Detector #eta","jet E [GeV]","Tower E [GeV]","Tower E [GeV]","Jet E [GeV]","Tower E [GeV]","Jet E_{T} [GeV]","Tower E_{T} [GeV]","Tower E_{T} [GeV]","Jet E_{T} [GeV]","Tower E_{T} [GeV]","Jet E [GeV]","Jet E_{T} [GeV]","Jet #eta","Jet #eta","Jet #phi","Jet #phi","E_{jet}^{EM}/E_{jet}","E_{jet}^{EM}/E_{jet}"};

  string ytitles[nhist] = {"Jet #phi","Jet #phi","Tower #phi","Tower #phi","Tower #phi","Tower #phi","Jet #eta","Tower #eta","Tower Detector #eta","Jet #phi","Tower #phi","Jet #eta","Tower #eta","Tower Detector #eta","Jet #phi","Tower #phi","E_{jet}^{EM}/E_{jet}","#Delta#phi","E_{jet}^{EM}/E_{jet}","E_{jet}^{EM}/E_{jet}","E_{jet}^{EM}/E_{jet}","E_{jet}^{EM}/E_{jet}","E_{jet}^{OH}/E_{jet}","E_{jet}^{OH}/E_{jet}"};

  for(int i=0; i<nhist; ++i)
    {
      get_and_draw_th2d(histnames[i],region,histfile,xtitles[i],ytitles[i],"Counts / N_{bit18}^{scaled}",can,ratcan,dijetcut);
    }

  const int nth1d = 9;
  string th1dnames[nth1d] = {"zhist","mbdn","mbds","mbdt","calo_hitsgrone_0","calo_hitsgrone_1","calo_hitsgrone_2","zhist_gr15","zhist_nocut"};
  string th1dxtitl[nth1d] = {"z_{vtx} [cm]","MBD Charge [Arb.]","MBD Charge [Arb.]","MBD Charge [Arb.]","Number of Towers with E > 1 GeV","Number of Towers with E > 1 GeV","Number of Towers with E > 1 GeV","z_{vtx} [cm]","z_{vtx} [cm]"};

  for(int i=0; i<nth1d; ++i)
    {
      get_and_draw_th1d(th1dnames[i],region,histfile,th1dxtitl[i],"Counts / N_{bit18}^{scaled}",can,ratcan,dijetcut);
    }
  delete can;
  delete ratcan;
  return 0;
}
