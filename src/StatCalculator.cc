#include <cmath>
#include <queue>
#include <boost/thread.hpp>

#include "StatCalculator.hh"
#include "HypothesisTest.hh"

using namespace std;
using namespace boost;

class CalcMatrix
{
  uint nthr_;
  uint thrno_;
  ZscoreMatrix& zm_;
  const map<uint,bool>& used_;
  const map<uint,vector<uint> >& bin2elm_;
  const Heatmap& hm_;
  const Bed& elm_;
  const Bed& gap_;
  uint ncnk_;
  uint seed_;
  uint nsmp_;
  bool verbose_;

public:
  CalcMatrix(uint nthr, uint thrno, ZscoreMatrix& zm, 
	     const map<uint,bool>& used, 
	     const map<uint,vector<uint> >& bin2elm, 
	     const Heatmap& hm, const Bed& elm, const Bed& gap, 
	     uint ncnk, uint seed, uint nsmp, bool verbose)
    : nthr_(nthr), thrno_(thrno), zm_(zm), 
      used_(used), 
      bin2elm_(bin2elm), 
      hm_(hm), elm_(elm), gap_(gap), 
      ncnk_(ncnk), seed_(seed), nsmp_(nsmp), verbose_(verbose) {}

  void operator()()
  {
    uint cnt = 0;
    // make sure incrementing them coordinately
    map<uint,vector<uint> >::const_iterator itr, jtr;
    uint i, j;
    for (itr=bin2elm_.begin(), i=0; itr!=bin2elm_.end(); ++itr, ++i) {
      //if (used_[itr->first]) continue; // wrong since [] is not const-compatible
      if ( ( used_.find(itr->first) )->second ) continue; 
      for (jtr=bin2elm_.begin(), j=0; jtr!=itr; ++jtr, ++j) {
	//if (used_[jtr->first]) continue; 
	if ( ( used_.find(jtr->first) )->second ) continue;
	if (cnt++ % nthr_ == thrno_) {
	  if (verbose_) cout << "i" << flush;
	  Bed tmp;
	  const vector<uint>& ielm = itr->second;
	  const vector<uint>& jelm = jtr->second;
	  for (uint idx=0; idx!=ielm.size(); ++idx) 
	    tmp.push_back(elm_, ielm[idx]);
	  for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
	    tmp.push_back(elm_, jelm[jdx]);
	  HypothesisTest htest(ncnk_, seed_, false, hm_, tmp, gap_);
	  double stat, p, z;
	  htest.test_stat();
	  htest.null_dist(nsmp_, false);
	  htest.get_stat(stat, p, z);
	  zm_(i,j) = zm_(j,i) = z;
	}
      }
    }
  }

};

void ZscoreMatrix::
calculate(const map<uint,bool>& used, 
	  const map<uint,vector<uint> >& bin2elm, 
	  const Heatmap& hm, const Bed& elm, const Bed& gap, 
	  uint ncnk, uint seed, uint nsmp, bool verbose, uint nthr)
{
  assert(bin2elm.size()==used.size());
  resize(bin2elm.size(), bin2elm.size());

  if (nthr > 1) {
    vector<thread*> thr(nthr );
    vector<CalcMatrix*> calc(nthr);
    for (uint t=0; t!=nthr; ++t) {
      calc[t] = new CalcMatrix(nthr, t, *this, used, bin2elm,
			       hm, elm, gap, ncnk, seed, nsmp, verbose);
      thr[t] = new boost::thread(*calc[t]);
    }
    for (uint t=0; t!=nthr; ++t) {
      thr[t]->join();
      delete calc[t];
      delete thr[t];
    }
  } 
  else {
    CalcMatrix calc(1, 0, *this, used, bin2elm, 
		    hm, elm, gap, ncnk, seed, nsmp, verbose);
    calc();
  }

  if (verbose) cout << endl << flush;
}

void ZscoreMatrix::
print()
{
  cout << "Zscore Matrix:" << endl << flush;
  for (uint i=0; i!=Nrow_; ++i) {
    for (uint j=0; j!=Ncol_-1; ++j) 
      cout << Zscore_[i][j] << "\t" << flush;      
    cout << Zscore_[i][Ncol_-1] << endl << flush;
  }
}

class CalcArray
{
  uint nthr_;
  uint thrno_;
  ZscoreArray& za_;
  const map<uint,bool>& used_;
  const Bed& preset_; 
  const map<uint,vector<uint> >& bin2elm_; 		 
  const Heatmap& hm_;
  const Bed& elm_;
  const Bed& gap_;
  uint ncnk_;
  uint seed_;
  uint nsmp_;
  bool verbose_;

public:
  CalcArray(uint nthr, uint thrno, ZscoreArray& za,
	    const map<uint,bool>& used, const Bed& preset, 
	    const map<uint,vector<uint> >& bin2elm,
	    const Heatmap& hm, const Bed& elm, const Bed& gap, 
	    uint ncnk, uint seed, uint nsmp, bool verbose)
    : nthr_(nthr), thrno_(thrno), za_(za), 
      used_(used), preset_(preset),
      bin2elm_(bin2elm),
      hm_(hm), elm_(elm), gap_(gap), 
      ncnk_(ncnk), seed_(seed), nsmp_(nsmp), verbose_(verbose) {}

  void operator()()
  {
    uint cnt = 0;
    // make sure incrementing them coordinately
    map<uint,vector<uint> >::const_iterator itr;
    uint i;
    for (itr=bin2elm_.begin(), i=0; itr!=bin2elm_.end(); ++itr, ++i) {
      //if (used_[itr->first]) continue; // wrong since [] is not const-compatible
      if ( ( used_.find(itr->first) )->second ) continue; 
      if (cnt++ % nthr_ == thrno_) {
	if (verbose_) cout << "." << flush;
	Bed tmp = preset_;
	const vector<uint>& ielm = itr->second;
	for (uint idx=0; idx!=ielm.size(); ++idx) 
	  tmp.push_back(elm_, ielm[idx]);
	HypothesisTest htest(ncnk_, seed_, false, hm_, tmp, gap_);
	double stat, p, z;
	htest.test_stat();
	htest.null_dist(nsmp_, false);
	htest.get_stat(stat, p, z);
	za_(i) = z;
      }
    }
  }

};

void ZscoreArray::
calculate(const map<uint,bool>& used, const Bed& preset, 
	  const map<uint,vector<uint> >& bin2elm, 
	  const Heatmap& hm, const Bed& elm, const Bed& gap, 
	  uint ncnk, uint seed, uint nsmp, bool verbose, uint nthr)
{
  assert(bin2elm.size()==used.size());
  resize(bin2elm.size());

  if (nthr > 1) {
    vector<thread*> thr(nthr );
    vector<CalcArray*> calc(nthr);
    for (uint t=0; t!=nthr; ++t) {
      calc[t] = new CalcArray(nthr, t, *this, used, preset, bin2elm, 
			      hm, elm, gap, ncnk, seed, nsmp, verbose);
      thr[t] = new boost::thread(*calc[t]);
    }
    for (uint t=0; t!=nthr; ++t) {
      thr[t]->join();
      delete calc[t];
      delete thr[t];
    }
  } 
  else {
    CalcArray calc(1, 0, *this, used, preset, bin2elm, 
		   hm, elm, gap, ncnk, seed, nsmp, verbose);
    calc();
  }

  if (verbose) cout << endl << flush;
}

void ZscoreArray::
print()
{
  cout << "Zscore Array:" << endl << flush;
  for (uint i=0; i!=Size_-1; ++i) 
    cout << Zscore_[i] << "\t" << flush;      
  cout << Zscore_[Size_-1] << endl << flush;
}
