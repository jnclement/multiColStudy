#include <dlUtility.h>

int format_th1d(TH1D* hist, string xtitle, string ytitle)
{
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.5);
  hist->SetMarkerColor(kCyan+3);

  return 0;
}

int draw_th1d(TH1D* hist, TCanvas* can)
{
  can->cd();
  can->SetRightMargin(0.15);
  can->SetBottomMargin(0.15);
  can->SetLeftMargin(0.15);
  can->SetTopMargin(0.15);

  hist->Draw("PE");
  sqrt_s_text(0.9,0.93);
  sphenixtext(0.9,0.97);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+".pdf").c_str());

  gPad->SetLogz();
  hist->Draw("PE");
  sqrt_s_text(0.9,0.93);
  sphenixtext(0.9,0.97);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_log.pdf").c_str());
  
  return 0;
}

int draw_ratio_th1d(TH1D** hists, string histbasename, int* rn, TCanvas* can)
{
  int nxbin = hists[0]->GetNbinsX();
  float xmax = hists[0]->GetXaxis()->GetXmax();
  float xmin = hists[0]->GetXaxis()->GetXmin();
  TH1D* ratio = new TH1D((histbasename+"_"+to_string(rn[0])+"_to_"+to_string(rn[1])).c_str(),"",nxbin,xmin,xmax);

  ratio->Divide(hists[0],hists[1]);

  format_th1d(ratio, hists[0]->GetXaxis()->GetTitle(), ("Ratio of "+to_string(rn[0])+" to "+to_string(rn[1])).c_str());
  
  draw_th1d(ratio, can);
    
  return 0;
}

int get_and_draw_th1d(string histbasename, int* rn, TFile** histfile, string xtitle, string ytitle, TCanvas* can)
{
  TH1D* hists[2];
  for(int i=0; i<2; ++i)
    {
      hists[i] = (TH1D*)histfile[i]->Get((histbasename+"_"+to_string(rn[i])).c_str());
      format_th1d(hists[i], xtitle, ytitle);
      draw_th1d(hists[i], can);
    }

  draw_ratio_th1d(hists, histbasename, rn, can);

  return 0;
}
int format_th2d(TH2D* hist, string xtitle, string ytitle, string ztitle)
{
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetZaxis()->SetTitle(ztitle.c_str());
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetLabelSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);

  return 0;
}

int draw_th2d(TH2D* hist, TCanvas* can)
{
  can->cd();
  can->SetRightMargin(0.15);
  can->SetBottomMargin(0.15);
  can->SetLeftMargin(0.15);
  can->SetTopMargin(0.15);

  hist->Draw("COLZ");
  sqrt_s_text(0.9,0.93);
  sphenixtext(0.9,0.97);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+".pdf").c_str());

  gPad->SetLogz();
  hist->Draw("COLZ");
  sqrt_s_text(0.9,0.93);
  sphenixtext(0.9,0.97);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_log.pdf").c_str());
  
  return 0;
}

int draw_ratio_th2d(TH2D** hists, string histbasename, int* rn, TCanvas* can)
{
  int nxbin = hists[0]->GetNbinsX();
  int nybin = hists[0]->GetNbinsY();
  float xmax = hists[0]->GetXaxis()->GetXmax();
  float ymax = hists[0]->GetYaxis()->GetXmax();
  float xmin = hists[0]->GetXaxis()->GetXmin();
  float ymin = hists[0]->GetYaxis()->GetXmin();
  TH2D* ratio = new TH2D((histbasename+"_"+to_string(rn[0])+"_to_"+to_string(rn[1])).c_str(),"",nxbin,xmin,xmax,nybin,ymin,ymax);

  ratio->Divide(hists[0],hists[1]);

  format_th2d(ratio, hists[0]->GetXaxis()->GetTitle(), hists[0]->GetYaxis()->GetTitle(), ("Ratio of "+to_string(rn[0])+" to "+to_string(rn[1])).c_str());
  
  draw_th2d(ratio, can);
    
  return 0;
}

int get_and_draw_th2d(string histbasename, int* rn, TFile** histfile, string xtitle, string ytitle, string ztitle, TCanvas* can)
{
  TH2D* hists[2];
  TH1D* projx[2];
  TH1D* projy[2];
  for(int i=0; i<2; ++i)
    {
      hists[i] = (TH2D*)histfile[i]->Get((histbasename+"_"+to_string(rn[i])).c_str());
      format_th2d(hists[i], xtitle, ytitle, ztitle);
      draw_th2d(hists[i], can);
      projx[i] = hists[i]->ProjectionX();
      projy[i] = hists[i]->ProjectionY();

      format_th1d(projx[i],hists[i]->GetXaxis()->GetTitle(),hists[i]->GetZaxis()->GetTitle());
      format_th1d(projy[i],hists[i]->GetYaxis()->GetTitle(),hists[i]->GetZaxis()->GetTitle());

      draw_th1d(projx[i],can);
      draw_th1d(projy[i],can);
    }

  draw_ratio_th1d(projx, histbasename+"_projx", rn, can);
  draw_ratio_th1d(projy, histbasename+"_projy", rn, can);
  draw_ratio_th2d(hists, histbasename, rn, can);

  return 0;
}

int plot(string histfilename0, string histfilename1, int rn0, int rn1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile* histfile[2];
  histfile[0] = TFile::Open(histfilename0.c_str());
  histfile[1] = TFile::Open(histfilename1.c_str());
  int rn[2] = {rn0, rn1};

  const int nhist = 18;

  TCanvas* can = new TCanvas("","",1000,1000);
  
  string histnames[nhist] = {"jet_eta_phi","jet_eta_phi_gr15","tow_eta_phi","tow_eta_phi_gr5","tow_deteta_phi","tow_deteta_phi_gr5","jet_E_eta","tow_E_eta","tow_E_deteta","jet_E_phi","tow_E_phi","jet_ET_eta","tow_ET_eta","tow_ET_deteta","jet_ET_phi","tow_ET_phi","jet_E_frcem","jet_ET_dphi"};

  string xtitles[nhist] = {"Jet #eta","Jet #eta","Tower #eta","Tower #eta","Tower Detector #eta","Tower Detector #eta","jet E [GeV]","Tower E [GeV]","Tower E [GeV]","Jet E [GeV]","Tower E [GeV]","Jet E_{T} [GeV]","Tower E_{T} [GeV]","Tower E_{T} [GeV]","Jet E_{T} [GeV]","Tower E_{T} [GeV]","Jet E [GeV]","Jet E_{T} [GeV]"};

  string ytitles[nhist] = {"Jet #phi","Jet #phi","Tower #phi","Tower #phi","Tower #phi","Tower #phi","Jet #eta","Tower #eta","Tower Detector #eta","Jet #phi","Tower #phi","Jet #eta","Tower #eta","Tower Detector #eta","Jet #phi","Tower #phi","E_{jet}^{EM}/E_{jet}","#Delta#phi"};

  for(int i=0; i<nhist; ++i)
    {
      get_and_draw_th2d(histnames[i],rn,histfile,xtitles[i],ytitles[i],"Event Normalized Counts",can);
    }

  const int nth1d = 7;
  string th1dnames[nth1d] = {"zhist","mbdn","mbds","mbdt","calo_hitsgrone_0","calo_hitsgrone_1","calo_hitsgrone_2"};
  string th1dxtitl[nth1d] = {"z_{vtx} [cm]","MBD Charge [Arb.]","MBD Charge [Arb.]","MBD Charge [Arb.]","Number of Towers with E > 1 GeV","Number of Towers with E > 1 GeV","Number of Towers with E > 1 GeV"};

  for(int i=0; i<nth1d; ++i)
    {
      get_and_draw_th1d(th1dnames[i],rn,histfile,th1dxtitl[i],"Event Normalized Counts",can);
    }
  
  return 0;
}
