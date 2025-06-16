#include <dlUtility.h>

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
}
vector<TH1D*> make_projections(TH3D* h, int axis, double to, double from, double lo1, double hi1, double lo2, double hi2)
{
  TCanvas* can = new TCanvas("","",2000,1000);
  string projdir;
  vector<TH1D*> outhists = {};

  
  if(axis==0) projdir = "yzo";
  else if(axis==1) projdir = "xzo";
  else if(axis==2) projdir = "xyo";
  else return 1;

  TH2D* proj1, proj2;
  
  if(to > from)
    {
      set_axis_range_based_on_axis(h, axis, get_bin_number_for_axis(h, axis, from, 1), get_bin_number_for_axis(h, axis, to, 0));
      proj1 = h->Project3D(projdir.c_str());
    }
  else
    {
      set_axis_range_based_on_axis(h, axis, 0, get_bin_number_for_axis(h, axis, to, 0));
      proj1 = h->Project3D(projdir.c_str());

      set_axis_range_based_on_axis(h, axis, get_bin_number_for_axis(h, axis, from, 1), 9999999);
      proj2 = h->Project3D(projdir.c_str());
      proj1->Add(proj2);
    }

  proj1->GetZaxis()->SetTitle("Counts");

  

  outhists.push_back(proj1->ProjectionX());
  outhists.push_back(proj1->ProjectionY());
  
  return outhists;
}

int plot_th3d()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  return 0;
}
