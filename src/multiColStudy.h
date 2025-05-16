#ifndef MULTICOLSTUDY_H
#define MULTICOLSTUDY_H
#include <fun4all/SubsysReco.h>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include <globalvertex/GlobalVertex.h>
class PHCompositeNode;
class CentralityInfo;
class multiColStudy : public SubsysReco
{
 public:

  multiColStudy(const std::string &filename, const std::string &name = "multiColStudy", const int debug = 0, int issim = 1);

  virtual ~multiColStudy();
  
  int Init(PHCompositeNode *topNode) override;
  
  int InitRun(PHCompositeNode *topNode) override;
  
  int process_event(PHCompositeNode *topNode) override;
  
  int ResetEvent(PHCompositeNode *topNode) override;
  
  int EndRun(const int runnumber) override;
  
  int End(PHCompositeNode *topNode) override;
  
  int Reset(PHCompositeNode * /*topNode*/) override;
  
  void Print(const std::string &what = "ALL") const override;
  
  
  private:

  std::string _name;
  int _debug;
  int _issim;
  std::string _filename;
  
  static const int _maxzvtx = 3;
  static const int _maxjet = 100;
  static const int _mbdside = 2;
  static const int _etow = 24576;
  static const int _htow = 1536;
  static const int _gronefield = 6;
  static const int _ncalotype = 3;
  static const int _maxjetcomp = 2000;
  TFile* _f;
  TTree* _tree;
  int _njet;
  double _jet_e[_maxjet];
  double _jet_et[_maxjet];
  double _jet_eta[_maxjet];
  double _jet_phi[_maxjet];
  double _jet_at[_maxjet];
  double _jet_at_em[_maxjet];
  double _jet_at_oh[_maxjet];
  double _frcem[_maxjet];
  double _frcoh[_maxjet];
  int _ncgroe[_maxjet];
  int _ncgroo[_maxjet];
  double _dphilead;
  int _isdijet;

  double _towgrone[_etow+2*_htow][_gronefield];

  double _ohat[_maxjet][64];
  double _emat[_maxjet][64];
  double _calo_emfrac;
  double _calo_ohfrac;
  double _calo_e;
  
  int _tnjet;
  double _tjet_e[_maxjet];
  double _tjet_et[_maxjet];
  double _tjet_eta[_maxjet];
  double _tjet_phi[_maxjet];
  double _tdphilead;
  int _tisdijet;

  int _nzvtx;
  int _tnzvtx;
  double _rzvtx[_maxzvtx];
  double _tzvtx[_maxzvtx];

  double _mbdq[_mbdside][64];

  int _hitsgrone;

  long long unsigned int _trigvec;
  GlobalVertex::VTXTYPE _vtxtype = GlobalVertex::MBD;
};

#endif // MULTICOLSTUDY
