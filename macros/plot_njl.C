#include <dlUtility.h>
const int nhistplot = 4;//4;
string region[nhistplot] = {"RegionA","RegionB","RegionC","RegionD"};//,"RegionA","RegionD"};//,"RegionA","RegionD","RegionA","RegionD"};
int triggers[nhistplot] = {18,18,18,18};//{18, 18,26,26};//, 17, 17, 19, 19};
int bit = 18;
int colors[nhistplot] = {kRed+1,kBlue+1,kGreen+1,kMagenta+1};
int markers[nhistplot] = {20,53,21,54};//,21,54,33,56};
double singlegaus(double* x, double* par)
{
  double exparg = (x[0]-par[1])/par[2];
  double result = par[0]*exp(-exparg*exparg);
  return result;
}
long long unsigned int globntot[nhistplot] = {0};
double globlumi[nhistplot] = {0};

vector<double> binom_error_hack(TH1D* num, TH1D* den)
{
  TH1D* numcp = (TH1D*)num->Clone();
  TH1D* dencp = (TH1D*)den->Clone();
  TH1D* dvded = new TH1D((std::string("dvd_")+num->GetName()+"_"+den->GetName()).c_str(),"",numcp->GetNbinsX(),numcp->GetBinLowEdge(1),numcp->GetBinLowEdge(numcp->GetNbinsX()+1));
  for(int i=1; i<numcp->GetNbinsX()+1; ++i)
    {
      
      if(abs(numcp->GetBinContent(i) - dencp->GetBinContent(i)) < 0.000001 && numcp->GetBinContent(i) > 1)
	{
	  numcp->SetBinContent(i,numcp->GetBinContent(i)-1);
	  numcp->SetBinError(i,sqrt(numcp->GetBinContent(i)-1));
	}
      else if(abs(numcp->GetBinContent(i) - dencp->GetBinContent(i)) < 0.000001)
	{
	  numcp->SetBinContent(i,1);
	  numcp->SetBinError(i,1);
	  dencp->SetBinContent(i,dencp->GetBinContent(i)+1);
	  dencp->SetBinError(i,sqrt(dencp->GetBinContent(i)+1));
	}
    }
  vector<double> retvec = {};
  dvded->Divide(numcp,dencp,1,1,"B");
  for(int i=1; i<numcp->GetNbinsX()+1; ++i)
    {
      cout << "bincontent num/den/binerror dvded " << i << " " << numcp->GetBinContent(i) << " " << dencp->GetBinContent(i) << " " << dvded->GetBinError(i) << endl;
      retvec.push_back(dvded->GetBinError(i));
    }
  return retvec;
}


double doublegaus(double* x, double* par)
{
  return singlegaus(x, par)+singlegaus(x,par+3);
}

TF1* fit_doublegaus(TH1D* hist, int n)
{
  TF1* dgaus = new TF1((hist->GetName()+to_string(n)).c_str(),doublegaus,-2,2,6);
  dgaus->SetParameters(0.25,0,1,0.075,-1,0.5);
  dgaus->SetParNames("A1","M1","W1","A2","M2","W2");
  hist->Fit(dgaus,"RI0");
  return dgaus;
}

double getMaxXatY(TH2D* hist, int ybin)
{
  int maxx = -1;
  double maxval = -9999999;
  int nxbin = hist->GetNbinsX();

  for(int i = 1; i < nxbin; ++i)
    {
      double content = hist->GetBinContent(i, ybin);
      if(content > maxval)
	{
	  maxval = content;
	  maxx = i;
	}
    }
  return hist->GetXaxis()->GetBinCenter(maxx);
}

int correct_eta_timing(TH2D* hist)
{

  int nbinx = hist->GetNbinsX();
  int nbiny = hist->GetNbinsY();
  float wxbin = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin())/nbinx;
  for(int i=1; i<nbiny+1; ++i)
    {
      double maxxy = getMaxXatY(hist, i);
      int nbinoff = maxxy>0?floor(abs(maxxy/wxbin)):-floor(abs(maxxy/wxbin));
      cout << "nbinoff: " << nbinoff << endl;
      int from = maxxy>0?1:nbinx;
      for(int j = from; (maxxy>0?-j:j)>(maxxy>0?-nbinx-1:0); (maxxy>0?j++:j--))
	{
	  cout << i << " " << j << " " << hist->GetBinContent(j,i) <<  " " << hist->GetBinContent(j+nbinoff,i) << endl;
	  hist->SetBinContent(j,i,hist->GetBinContent(j+nbinoff,i));
	  hist->SetBinError(j,i,hist->GetBinError(j+nbinoff,i));
	  cout << i << " " << j << " " << hist->GetBinContent(j,i) <<  " " << hist->GetBinContent(j+nbinoff,i) << endl;
	}
    }
  return 0;
}

int format_th1d(TH1D* hist, string xtitle, string ytitle, int n = 1)
{
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);


  hist->SetMarkerSize(1.5);
  hist->SetMarkerColor(colors[n]);
  hist->SetLineColor(colors[n]);
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(markers[n]);
  
  return 0;
}

int draw_overlay_with_ratio_th1d(TH1D** hists, string histbasename, TCanvas* can, int dijetcut, int divfunc = 0, TF1** func = NULL)
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
  TLegend* leg = new TLegend(0.6,0.65,0.8,0.85);
  //cout <<" overlay and ratio time" << endl;
  for(int i=0; i<nhistplot; ++i)
    {
      if(std::string(hists[i]->GetName()).find("spec") != std::string::npos || std::string(hists[i]->GetName()).find("lead") != std::string::npos) hists[i]->GetYaxis()->SetTitle("Counts / \\mathscr{L}_{int} [pb]");
      if(std::string(hists[i]->GetName()).find("lead") != std::string::npos || std::string(hists[i]->GetName()).find("lead") != std::string::npos) hists[i]->GetYaxis()->SetTitle("Counts / \\mathscr{L}_{int} [pb]");
      if(hists[i]->GetMaximum() > ymax) ymax = hists[i]->GetMaximum();
      leg->AddEntry(hists[i],(region[i]+" Trig"+to_string(triggers[i])).c_str(),"p");
      ratio[i] = new TH1D((histbasename+"_"+region[0]+"_ratio_"+region[i]+"_trigger_"+to_string(triggers[i])).c_str(),"",hists[i]->GetNbinsX(),hists[i]->GetXaxis()->GetXmin(),hists[i]->GetXaxis()->GetXmax());
      if(std::string(hists[i]->GetName()).find("mbdt") != std::string::npos)
	{
	  ymax = 0.2;
	  hists[i]->Rebin(5);
	}
      if(i>0 && !divfunc) ratio[i]->Divide(hists[i],hists[0]);
      else if(divfunc)
	{
	  ratio[i] = (TH1D*)hists[i]->Clone();
	  ratio[i]->Divide(func[i]);
	}
      ratio[i]->GetYaxis()->SetRangeUser(0,2);
      format_th1d(ratio[i],hists[i]->GetXaxis()->GetTitle(),"Ratio to Fit",i);
      ratio[i]->GetXaxis()->SetTitleSize(ratio[i]->GetXaxis()->GetTitleSize()/0.5);
      ratio[i]->GetYaxis()->SetTitleSize(ratio[i]->GetYaxis()->GetTitleSize()/0.5);
      ratio[i]->GetXaxis()->SetLabelSize(ratio[i]->GetXaxis()->GetLabelSize()/0.5);
      ratio[i]->GetYaxis()->SetLabelSize(ratio[i]->GetYaxis()->GetLabelSize()/0.5);
    }
  //cout << "did firstloop" << endl;
  TF1* dgaus[nhistplot];
  hists[0]->GetYaxis()->SetRangeUser(0,1.5*ymax);
  for(int i=0; i<nhistplot; ++i)
    {
      can->cd(1);
      if(std::string(hists[i]->GetName()).find("mbd") != std::string::npos) ymax = hists[i]->GetBinContent(11)*1.5;

      hists[i]->GetXaxis()->SetLabelSize(0);
      hists[i]->GetXaxis()->SetTitleSize(0);
      hists[i]->GetYaxis()->SetRangeUser(0,1.5*ymax);

      hists[i]->Draw(i==0?"PE":"SAME PE");
      can->cd(2);
      if(divfunc) ratio[i]->GetListOfFunctions()->Clear();
      if(i>0 && !divfunc) ratio[i]->Draw(i==1?"PE":"SAME PE");
      else if(divfunc) ratio[i]->Draw(i==0?"PE":"SAME PE");
      cout << "DREW HIST!!!" << endl;
    }

  can->cd(1);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  float minet = 15;
  if(triggers[0] == 18)
    {
  antikt_text(0.4,0.3,0.92);

  if(std::string(hists[0]->GetName()).find("gr20") != std::string::npos) minet = 20;
  //et_cut_text(minet,0.015,0.96);
  //dijet_cut_text(0.3,0.96);
    }
  TLine* oneline = new TLine(hists[0]->GetXaxis()->GetXmin(),1,hists[0]->GetXaxis()->GetXmax(),1);
  can->cd(2);
  oneline->Draw();
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hists[0]->GetName())+"_"+to_string(bit)+"_overlay_"+(divfunc?"divfunc_":"")+(dijetcut?"dc":"nc")+".png").c_str());


  can->cd(1);
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i]->GetYaxis()->SetRangeUser(hists[0]->GetBinContent(hists[0]->FindLastBinAbove(0))*0.5,1.5*ymax);
    }
  gPad->SetLogy();
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i]->Draw(i==0?"PE":"SAME PE");
    }
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  if(bit != 10)
    {
  antikt_text(0.4,0.3,0.92);

  if(std::string(hists[0]->GetName()).find("gr20") != std::string::npos) minet = 20;
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
    }
  leg->Draw();
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hists[0]->GetName())+"_"+to_string(bit)+"_overlay_"+(dijetcut?"dc":"nc")+"_log.png").c_str());
  gPad->SetLogy(0);
  ymax = 0;
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i]->Scale(1./hists[i]->Integral());
      hists[i]->GetYaxis()->SetTitle("Counts / N_{bit10}^{scaled} (Normalized to #Sigma(bins)=1)");

      if(hists[i]->GetMaximum() > ymax) ymax = hists[i]->GetMaximum();
      if(std::string(hists[i]->GetName()).find("mbd") != std::string::npos) ymax = hists[i]->GetBinContent(11)*1.5;
      if(std::string(hists[i]->GetName()).find("mbdt") != std::string::npos)
	{
	  ymax = 0.2;
	}
    }

  hists[0]->GetYaxis()->SetRangeUser(0,1.5*ymax);
  /*
  for(int i=0; i<nhistplot; ++i)
    {
      dgaus[i] = fit_doublegaus(hists[i],i);
    }
  */
  for(int i=0; i<nhistplot; ++i)
    {
      can->cd(1);
      hists[i]->GetXaxis()->SetLabelSize(0);
      hists[i]->GetXaxis()->SetTitleSize(0);
      hists[i]->GetYaxis()->SetRangeUser(0,1.5*ymax);
      if(std::string(hists[i]->GetName()).find("mbdt") != std::string::npos)
	{
	  ymax = 0.2;
	}
      if(i>0) ratio[i]->Divide(hists[i],hists[0]);
      //dgaus[i]->SetLineColor(colors[i]);
      hists[i]->Draw(i==0?"PE":"SAME PE");
      //dgaus[i]->Draw("SAME");
      can->cd(2);
      if(i>0) ratio[i]->GetYaxis()->SetRangeUser(0,2);
      if(i>0) ratio[i]->Draw(i==1?"PE":"SAME PE");
    }


  
  can->cd(1);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  if(bit != 10)
    {
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
    }
  can->cd(2);
  oneline->Draw();

  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hists[0]->GetName())+"_"+to_string(bit)+"_overlay_"+(dijetcut?"dc":"nc")+"_normed.png").c_str());



  
  
  delete leg;
  for(int i=0; i<nhistplot; ++i)
    {
      delete ratio[i];
    }
  return 0;
}

int draw_th1d(TH1D* hist, TCanvas* can, int dijetcut, int trigger = bit)
{
  can->cd();
  can->SetRightMargin(0.2);
  can->SetBottomMargin(0.15);
  can->SetLeftMargin(0.15);
  can->SetTopMargin(0.15);

  float minet = 15;
  if(std::string(hist->GetName()).find("gr20") != std::string::npos) minet = 20;
  hist->Draw("PE");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+to_string(trigger)+"_"+(dijetcut?"dc":"nc")+".png").c_str());

  gPad->SetLogy();
  hist->Draw("PE");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+to_string(trigger)+"_"+(dijetcut?"dc":"nc")+"_log.png").c_str());
  gPad->SetLogy(0);
  return 0;
}

int draw_ratio_th1d(TH1D** hists, string histbasename, string* region, TCanvas* can, int dijetcut)
{
  int nxbin = hists[0]->GetNbinsX();
  double xmax = hists[0]->GetXaxis()->GetXmax();
  double xmin = hists[0]->GetXaxis()->GetXmin();
  TH1D* ratio = new TH1D((histbasename+"_"+region[0]+"_to_"+region[1]).c_str(),"",nxbin,xmin,xmax);

  ratio->Divide(hists[1],hists[0]);

  format_th1d(ratio, hists[0]->GetXaxis()->GetTitle(), ("Ratio of "+region[1]+" to "+region[0]).c_str());
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
      if(std::string(hists[i]->GetName()).find("spec") == std::string::npos && std::string(hists[i]->GetName()).find("lead") == std::string::npos) hists[i]->Scale(1./globntot[i]);
      else hists[i]->Scale(1./globlumi[i]);
      if(std::string(hists[i]->GetName()).find("zhist") != std::string::npos) hists[i]->Rebin(10);

      format_th1d(hists[i], xtitle, ytitle,i);
      if(std::string(hists[i]->GetName()).find("trigturn") != std::string::npos) hists[i]->GetYaxis()->SetTitle("Efficiency");;
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
  if(std::string(hist->GetName()).find("gr20") != std::string::npos) minet = 20;
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+to_string(bit)+"_"+(dijetcut?"dc":"nc")+".png").c_str());
  gPad->SetLogz();
  hist->Draw("COLZ");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  antikt_text(0.4,0.3,0.92);
  et_cut_text(minet,0.015,0.96);
  dijet_cut_text(0.3,0.96);
  //can->SaveAs(("/sphenix/user/jocl/projects/multiColStudy/output/plots/"+string(hist->GetName())+"_"+to_string(bit)+"_"+(dijetcut?"dc":"nc")+"_log.png").c_str());
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

  ratio->Divide(hists[1],hists[0]);
  format_th2d(ratio, hists[0]->GetXaxis()->GetTitle(), hists[0]->GetYaxis()->GetTitle(), ("Ratio of "+region[1]+" to "+region[0]).c_str());
  ratio->GetZaxis()->SetRangeUser(0,2);
  draw_th2d(ratio, can, dijetcut);
  delete ratio;
  return 0;
}

double smart_hist_max(TH1D* h)
{
  int nbins = h->GetNbinsX();

  std::vector<double> contents;
  for (int i = 1; i <= nbins; ++i) {
    double val = h->GetBinContent(i);
    if(val>0) contents.push_back(val);
  }
  
  // Sort in descending order
  std::sort(contents.begin(), contents.end(), std::greater<double>());
  
  // Use, for example, the 95th percentile as the effective maximum
  int index = contents.size() * 0.05;  // ignore top 5%
  double smart_max = contents[index];
  
  return smart_max;
}

double smart_hist_min(TH1D* h)
{
  int nbins = h->GetNbinsX();

  std::vector<double> contents;
  for (int i = 1; i <= nbins; ++i) {
    double val = h->GetBinContent(i);
    if(val>0) contents.push_back(val);
  }
  
  // Sort in descending order
  std::sort(contents.begin(), contents.end(), std::less<double>());
  
  // Use, for example, the 95th percentile as the effective maximum
  int index = contents.size() * 0.05;  // ignore top 5%
  double smart_min = contents[index];
  
  return smart_min;
}


int get_and_draw_th2d(string histbasename, string* region, TFile* histfile, string xtitle, string ytitle, string ztitle, TCanvas* can, TCanvas* ratcan, int dijetcut)
{
  TH2D* hists[nhistplot];
  TH1D* projx[nhistplot];
  TH1D* projy[nhistplot];
  for(int i=0; i<nhistplot; ++i)
    {
      hists[i] = (TH2D*)histfile->Get((histbasename+"_"+region[i]).c_str());
      cout << globntot[i] << endl;
      if(std::string(hists[i]->GetName()).find("spec") == std::string::npos && std::string(hists[i]->GetName()).find("lead") == std::string::npos) hists[i]->Scale(1./globntot[i]);
      else hists[i]->Scale(1./globlumi[i]);
      //if(std::string(hists[i]->GetName()).find("frcem") != std::string::npos) hists[i]->Rebin2D(5,5);
      //if(std::string(hists[i]->GetName()).find("calo") != std::string::npos) hists[i]->Rebin2D(5,5);
      if(std::string(hists[i]->GetName()).find("tgrone_eta_2") != std::string::npos)
	{
	  //correct_eta_timing(hists[i]);
	}
      /*
      else
	{
	  cout << histbasename << " " << hists[i]->GetName() << endl;
	  return 0;
	}
      */
      format_th2d(hists[i], xtitle, ytitle, ztitle);
      //if(std::string(hists[i]->GetName()).find("tgrone") != std::string::npos)
      draw_th2d(hists[i], can,dijetcut);

      projx[i] = hists[i]->ProjectionX((histbasename+"_projx"+to_string(i)).c_str(),1,hists[i]->GetNbinsX());
      projy[i] = hists[i]->ProjectionY((histbasename+"_projy"+to_string(i)).c_str(),1,hists[i]->GetNbinsY());
      if(std::string(hists[i]->GetName()).find("frcem") != std::string::npos) cout << "frcem integral " << hists[i]->GetName() << ": " << projx[i]->Integral() << endl;
      format_th1d(projx[i],hists[i]->GetXaxis()->GetTitle(),hists[i]->GetZaxis()->GetTitle(),i);
      format_th1d(projy[i],hists[i]->GetYaxis()->GetTitle(),hists[i]->GetZaxis()->GetTitle(),i);

      //draw_th1d(projx[i],can,dijetcut);
      //draw_th1d(projy[i],can,dijetcut);
    }

  //draw_ratio_th1d(projx, histbasename+"_projx", region, can, dijetcut);
  //draw_ratio_th1d(projy, histbasename+"_projy", region, can, dijetcut);
  draw_ratio_th2d(hists, histbasename, region, can, dijetcut);
  draw_overlay_with_ratio_th1d(projx, histbasename+"_projx",ratcan,dijetcut);
  draw_overlay_with_ratio_th1d(projy, histbasename+"_projy",ratcan,dijetcut);
  if(std::string(hists[0]->GetName()).find("jet_at_em_at_oh") != std::string::npos)
    {
      for(int i=0; i<nhistplot; ++i)
	{
	  cout << "jet at EM OH means: " << projx[i]->GetMean() << " " << projy[i]->GetMean() << endl;
	}
    }
  return 0;
}

int plot_njl(int tb)
{
  cout << "test" << endl;
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  bit = tb;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile* histfile;
  histfile = TFile::Open(("../output/hists/hadded"+to_string(bit)+".root").c_str());

  TCanvas* lumican = new TCanvas("","",2000,1000);

  string histname1 = "njet_lumi";
  string histname2 = "njet_lumi_Ecut";

  if(bit==26)
    {
      histname1 = "nclus_lumi";
      histname2 = "nclus_lumi_Ecut";
    }
  
  TH1D* njet_lumi = (TH1D*)histfile->Get(histname1.c_str());

  double njlsmin = smart_hist_min(njet_lumi);
  double njlsmax = smart_hist_max(njet_lumi);
  njet_lumi->GetYaxis()->SetRangeUser(0.75*njlsmin,1.25*njlsmax);

  
  //njet_lumi->GetYaxis()->SetRangeUser(0,7e4);
  njet_lumi->SetMarkerSize(1);
  njet_lumi->SetMarkerStyle(20);
  lumican->cd();
  njet_lumi->GetYaxis()->SetTitle(("N_{jet}/\\mathscr{L}_{int,corr}^{trig"+to_string(tb)+"}").c_str());
  njet_lumi->GetXaxis()->SetTitle("Run Number");
  njet_lumi->Draw("PE");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  if(bit!=26) antikt_text(0.4,0.3,0.92);
  if(bit!=26) et_cut_text(10,0.015,0.96);
  if(bit!=26) dijet_cut_text(0.3,0.96);
  lumican->SaveAs(("../output/plots/"+to_string(bit)+"_njet_lumi.png").c_str());

  TH1D* njet_lumi_Ecut = (TH1D*)histfile->Get(histname2.c_str());

  njlsmin = smart_hist_min(njet_lumi_Ecut);
  njlsmax = smart_hist_max(njet_lumi_Ecut);
  njet_lumi_Ecut->GetYaxis()->SetRangeUser(0.75*njlsmin,1.25*njlsmax);
  
  //njet_lumi_Ecut->GetYaxis()->SetRangeUser(0,1.5e4);
  njet_lumi_Ecut->SetMarkerSize(1);
  njet_lumi_Ecut->SetMarkerStyle(20);
  lumican->cd();
  njet_lumi_Ecut->GetYaxis()->SetTitle(("N_{jet}/\\mathscr{L}_{int,corr}^{trig"+to_string(tb)+"}").c_str());
  njet_lumi_Ecut->GetXaxis()->SetTitle("Run Number");
  njet_lumi_Ecut->Draw("PE");
  sphenixtext(0.65,0.96);
  sqrt_s_text(0.65,0.92);
  if(bit!=26) antikt_text(0.4,0.3,0.92);
  if(bit!=26) et_cut_text(15,0.015,0.96);
  if(bit!=26) dijet_cut_text(0.3,0.96);
  lumican->SaveAs(("../output/plots/"+to_string(bit)+"_njet_lumi_Ecut.png").c_str());


  return 0;
}
