#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/random.hpp>

#include "HypothesisTest.hh"
#include "StatCalculator.hh"

using namespace std;
using namespace boost;

ChromosomeChunk::
ChromosomeChunk(const vector<uint>& elm, uint start, uint stop) 
  : Elm_(elm), Nelm_(elm.size()), Start_(start), Stop_(stop)
{
  if (Nelm_ == 0) { 
    Dist_.resize(0);
    Dsum_ = 0;
    ShfElm_.resize(0);
    ShfDist_.resize(0);
    return;
  }

  Dist_.resize(Nelm_-1);
  for (uint i=0; i!=Nelm_-1; ++i)
    Dist_[i] = Elm_[i+1] - Elm_[i];
  Dsum_ = 0;
  for (uint i=0; i!=Nelm_-1; ++i)
    Dsum_ += Dist_[i];
  ShfElm_.assign(Elm_.begin(), Elm_.end());
  ShfDist_.assign(Dist_.begin(), Dist_.end());
}

void ChromosomeChunk::
shuffle(RandomGenerator& rng)
{
  if (Nelm_ == 0) return;

  uniform_int<> unif_dist(Start_, Stop_-1-Dsum_);
  variate_generator<RandomGenerator&, uniform_int<> > unif_samp(rng, unif_dist);
  RandomShuffler rns(rng);
  ShfElm_[0] = unif_samp();
  random_shuffle(ShfDist_.begin(), ShfDist_.end(), rns);

  for (uint i=1; i!=Nelm_; ++i) 
    ShfElm_[i] = ShfElm_[i-1] + ShfDist_[i-1];
}

HypothesisTest::
HypothesisTest(uint ncnk, uint seed, bool verbose, 
	       const Heatmap& hm, const Bed& elm, const Bed& gap)
  : Nchr_(hm.nchr()), Ncnk_(ncnk), 
    TestStat_(0.0), NullStat_(0.0), Pval_(0.0), 
    Hm_(hm), Elm_(elm), Gap_(gap), seed_(seed), rng_(seed)
{
  assert(Elm_.size() > 1);

  ShortArm_.resize(Nchr_);
  LongArm_.resize(Nchr_);

  // convert chromosome-wise nucleotide coordinates into 
  // whole-genome bin coordinates, separately for each chromosome
  vector<vector<uint> > GapBinStart(Nchr_);
  vector<vector<uint> > GapBinStop(Nchr_);
  vector<vector<uint> > ElmBin(Nchr_);

  for (uint i=0; i!=Gap_.size(); ++i) {
    uint chr_st, bin_st, chr_sp, bin_sp;
    if (Hm_.pos2bin(Gap_.chr(i), Gap_.start(i), chr_st, bin_st) && 
	Hm_.pos2bin(Gap_.chr(i), Gap_.stop(i), chr_sp, bin_sp)) {
      // In principle, gaps (i.e. centromeres) can be identified from 
      // a heatmap automatically by detecting Hm_.is_null()==true regions 
      // around chromosome midpoints. However, such regions do not always 
      // exist in actual heatmaps produced by frag_merge.py. Therefore, 
      // we re-define centromeres by ourselves using gap.txt. 
      // In the mirnylib and hiclib libraries, bowtie2 indexes include 
      // centromeres; thus, reads mapped in centromeres need to be 
      // removed afterward in frag_merge.py (more specifically, in 
      // mirnylib/genome.py) using gap.txt so that output heatmaps have 
      // Hm_.is_null()==true regions. I think this process might be 
      // broken in the mirnylib library. 
      /*
      for (uint j=bin_st-5; j!=bin_sp+5; ++j) {
	if (j==bin_st) cout << "start " << flush;
	cout << j << ":" << Hm_.is_null(j) << " " << flush;
	if (j==bin_sp) cout << "stop " << flush;
      }
      */
      if (verbose) cout << "registered gap" << i << " in bins" << bin_st << "," << bin_sp+1 << ":" << endl 
			<< Gap_.line(i);
      assert(chr_st == chr_sp);
      GapBinStart[chr_st].push_back(bin_st);
      GapBinStop[chr_sp].push_back(bin_sp+1);
    }
    else {
      if (verbose) cout << "skipped gap" << i << " in undefined regions:" << endl 
			<< Gap_.line(i);
    }
  }
  for (uint i=0; i!=Nchr_; ++i) {
    assert(GapBinStart[i].size() == 1);
    assert(GapBinStop[i].size() == 1);
  }

  for (uint i=0; i!=Elm_.size(); ++i) {
    uint start = (Elm_.strand(i) == "-") ? Elm_.stop(i) : Elm_.start(i); // "." is regarded as plus.
    uint chr, bin;
    if (Hm_.pos2bin(Elm_.chr(i), start, chr, bin)) {
      if (Hm_.is_null(bin)) {
	if (verbose) cout << "skipped element" << i << " in null regions:" << endl 
			  << Elm_.line(i);
      }
      else {
	if (bin >= GapBinStart[chr][0] && bin < GapBinStop[chr][0]) {
	  if (verbose) cout << "skipped element" << i << " in gap regions:" << endl 
			    << Elm_.line(i);
	}
	else {
	  if (verbose) cout << "registered element" << i << " in bin" << bin << ":" << endl 
			    << Elm_.line(i);
	  ElmBin[chr].push_back(bin);
	  Bin2Elm_[bin].push_back(i);
	}
      }
    }
    else {
      if (verbose) cout << "skipped element" << i << " in undefined regions:" << endl 
			<< Elm_.line(i);
    }
  }
  for (uint i=0; i!=Nchr_; ++i) {
    sort(ElmBin[i].begin(), ElmBin[i].end());
  }


  // separate elements into chromosome chunks
  for (uint i=0; i!=Nchr_; ++i) {
    uint short_st = Hm_.chr_start(i);
    uint short_sp = GapBinStart[i][0];
    uint long_st = GapBinStop[i][0]; 
    uint long_sp = Hm_.chr_stop(i);
    uint short_d = short_sp - short_st;
    uint long_d = long_sp - long_st;
    // at least one bin for each chunk 
    assert(short_st + Ncnk_ <= short_sp);
    assert(long_st + Ncnk_ <= long_sp);

    // short arm
    for (uint j=0; j!=Ncnk_; ++j) {
      vector<uint> elm;
      uint cnk_st = short_st + j * short_d / Ncnk_;
      uint cnk_sp = short_st + (j+1) * short_d / Ncnk_;
      for (uint k=0; k!=ElmBin[i].size(); ++k) 
	if (ElmBin[i][k] >= cnk_st && ElmBin[i][k] < cnk_sp)
	  elm.push_back(ElmBin[i][k]);
      ShortArm_[i].push_back(ChromosomeChunk(elm, cnk_st, cnk_sp));
    }
    // long arm
    for (uint j=0; j!=Ncnk_; ++j) {
      vector<uint> elm;
      uint cnk_st = long_st + j * long_d / Ncnk_;
      uint cnk_sp = long_st + (j+1) * long_d / Ncnk_;
      for (uint k=0; k!=ElmBin[i].size(); ++k) 
	if (ElmBin[i][k] >= cnk_st && ElmBin[i][k] < cnk_sp)
	  elm.push_back(ElmBin[i][k]);
      LongArm_[i].push_back(ChromosomeChunk(elm, cnk_st, cnk_sp));
    }
  }

  //print();
}

bool HypothesisTest::
test_stat(bool shuffled)
{
  bool verbose = false;

  vector<uint> elm;
  double sum = 0.0;
  uint n = 0;

  // collect elements
  for (uint i=0; i!=Nchr_; ++i) {
    for (uint j=0; j!=Ncnk_; ++j) {
      if (shuffled) {
	elm.insert(elm.end(), ShortArm_[i][j].ShfElm_.begin(), ShortArm_[i][j].ShfElm_.end());
	elm.insert(elm.end(), LongArm_[i][j].ShfElm_.begin(), LongArm_[i][j].ShfElm_.end());
      }
      else {
	elm.insert(elm.end(), ShortArm_[i][j].Elm_.begin(), ShortArm_[i][j].Elm_.end());
	elm.insert(elm.end(), LongArm_[i][j].Elm_.begin(), LongArm_[i][j].Elm_.end());
      }
    }
  }

  for (uint i=0; i!=elm.size(); ++i) {
    if (verbose) cout << "." << flush;
    vector<Heatmap::ValueType> row; 
    Hm_.get_row(row, elm[i]);
    for (uint j=i; j!=elm.size(); ++j) {
      if (uabs(elm[j], elm[i]) < 2) continue; 
      sum += row[elm[j]];
      n++;
    }
  }
  if (verbose) cout << endl << flush;

  double stat = (n==0) ? 0.0 : sum/n;
  if (shuffled) NullStat_ = stat;
  else TestStat_ = stat;

  if (verbose) cout << "#element: " << elm.size() << ", "
		    << "#comparison: " << n << ", "
		    << "statistics: " << stat << endl;

  return true;
}

bool HypothesisTest::
null_stat()
{
  return test_stat(true);
}

bool HypothesisTest::
null_dist(uint nsmp, bool verbose)
{
  NullDist_.clear();
  NullDist_.resize(nsmp);

  for (uint idx=0; idx!=nsmp; ++idx) {
    if (verbose) cout << "Monte Carlo sample " << idx << endl << flush;
    for (uint i=0; i!=Nchr_; ++i) {
      // short arm
      for (int j=0; j!=Ncnk_; ++j) { // use int for j--
	bool redo = false;
	ShortArm_[i][j].shuffle(rng_);
	for (uint k=0; k!=ShortArm_[i][j].Nelm_; ++k) {
	  if (Hm_.is_null(ShortArm_[i][j].ShfElm_[k])) {
	    redo = true;
	    break;
	  }
	}
	if (redo) {
	  if (verbose) cout << "r" << flush;
	  j--;
	  continue;
	}
      }
      // long arm
      for (int j=0; j!=Ncnk_; ++j) { // use int for j--
	bool redo = false;
	LongArm_[i][j].shuffle(rng_);
	for (uint k=0; k!=LongArm_[i][j].Nelm_; ++k) {
	  if (Hm_.is_null(LongArm_[i][j].ShfElm_[k])) {
	    redo = true;
	    break;
	  }
	}
	if (redo) {
	  if (verbose) cout << "r" << flush;
	  j--;
	  continue;
	}
      }
    }

    if (verbose) cout << endl << flush;
    null_stat();
    NullDist_[idx] = NullStat_;
  }

  sort(NullDist_.begin(), NullDist_.end());

  return true;
}

bool HypothesisTest::
output_stat(const string& filename) 
{
  calc_pval();
  calc_zval();

  ofstream ofs(filename.c_str());
  if (!ofs) return false;

  ofs << "Pval:" << endl << Pval_ << endl;
  ofs << "Zval:" << endl << Zval_ << endl;
  ofs << "TestStat:" << endl << TestStat_ << endl;

  return true;
}

bool HypothesisTest::
output_dist(const string& filename) 
{
  ofstream ofs(filename.c_str());
  if (!ofs) return false;

  ofs << "NullDist:" << endl;
  for (uint i=0; i!=NullDist_.size(); ++i) 
    ofs << NullDist_[i] << endl;
  ofs << endl;

  return true;
}

bool HypothesisTest::
get_stat(double& stat, double& p, double& z) 
{
  calc_pval();
  calc_zval();

  stat = TestStat_;
  p = Pval_;
  z = Zval_;

  return true;
}

void HypothesisTest::
print()
{
  cout << "print HypothesisTest class:" << endl << endl << flush;

  cout << "Nchr:" << endl << Nchr_ << endl << flush;  

  cout << "Ncnk:" << endl << Ncnk_ << endl << flush;  

  for (uint i=0; i!=Nchr_; ++i) {
    cout << "Chromosome " << i << endl << endl << flush;

    cout << "Short arm" << endl << flush;
    for (uint j=0; j!=Ncnk_; ++j) {
      cout << "Chunk " << j << " " 
	   << ShortArm_[i][j].Start_ << "-" << ShortArm_[i][j].Stop_ << endl << flush;
      for (uint k=0; k!=ShortArm_[i][j].Nelm_; ++k) 
	cout << ShortArm_[i][j].Elm_[k] << " " << flush;
      cout << endl << flush;
    }
    cout << "Long arm" << endl << flush;
    for (uint j=0; j!=Ncnk_; ++j) {
      cout << "Chunk " << j << " " 
	   << LongArm_[i][j].Start_ << "-" << LongArm_[i][j].Stop_ << endl << flush;
      for (uint k=0; k!=LongArm_[i][j].Nelm_; ++k) 
	cout << LongArm_[i][j].Elm_[k] << " " << flush;
      cout << endl << flush;
    }
  }
}

void HypothesisTest::
search_subset(Bed& subset, double& zmax, uint nsmp, bool verbose, uint nthr)
{
  static bool init = true;
  static map<uint,bool> used;
  static ZscoreMatrix zmat;

  // make sure incrementing them coordinately
  map<uint,vector<uint> >::iterator itr, jtr;
  uint i, j;

  if (init) {
    for (itr=Bin2Elm_.begin(); itr!=Bin2Elm_.end(); ++itr) {
      assert(used.find(itr->first)==used.end());
      used[itr->first] = false;
    }
    zmat.calculate(used, Bin2Elm_, Hm_, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);
    //zmat.print();
    init = false;
  }

  bool last;
  map<uint,vector<uint> >::iterator itrdx, jtrdx;  

  progress("the initial pair of bins");
  last = true;
  itrdx = jtrdx = Bin2Elm_.begin();
  for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
    if (used[itr->first]) continue;
    for (jtr=Bin2Elm_.begin(), j=0; jtr!=itr; ++jtr, ++j) {
      if (used[jtr->first]) continue;
      if (zmat(i,j) > zmax) {
	zmax = zmat(i,j);
	itrdx = itr;
	jtrdx = jtr;
	last = false;
      }
    }
  }
  if (last) return; 
  const vector<uint>& ielm = itrdx->second;
  const vector<uint>& jelm = jtrdx->second;
  for (uint idx=0; idx!=ielm.size(); ++idx) 
    subset.push_back(Elm_, ielm[idx]);
  for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
    subset.push_back(Elm_, jelm[jdx]);
  used[itrdx->first] = true;
  used[jtrdx->first] = true;
  if (verbose) {
    cout << "z = " << zmax << endl;
    cout << "bin" << itrdx->first << endl;
    for (uint idx=0; idx!=ielm.size(); ++idx) 
      cout << Elm_.line(ielm[idx]) << flush;
    cout << "bin" << jtrdx->first << endl;
    for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
      cout << Elm_.line(jelm[jdx]) << flush;
  }

  for (uint k=2; k<Bin2Elm_.size(); ++k) {
    progress("the next bin");
    last = true;
    itrdx = jtrdx = Bin2Elm_.begin();

    ZscoreArray zarr;
    zarr.calculate(used, subset, Bin2Elm_, Hm_, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);
    for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
      if (used[itr->first]) continue;
      if (zarr(i) >= zmax) {
	zmax = zarr(i);
	itrdx = itr;
	last = false;
      }
    }
    if (last) return;
    const vector<uint>& ielm = itrdx->second;
    for (uint idx=0; idx!=ielm.size(); ++idx) 
      subset.push_back(Elm_, ielm[idx]);
    used[itrdx->first] = true;
    if (verbose) {
      cout << "z = " << zmax << endl;
      cout << "bin" << itrdx->first << endl;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	cout << Elm_.line(ielm[idx]) << flush;
    }
  }

  return;
}

void HypothesisTest::
diff_search_subset(Bed& subset, double& zmax, uint nsmp, bool verbose, uint nthr, 
		   const Heatmap& ctl)
{
  static bool init = true;
  static map<uint,bool> used;
  static ZscoreMatrix zmat, zmat_ctl;

  // make sure incrementing them coordinately
  map<uint,vector<uint> >::iterator itr, jtr;
  uint i, j;

  if (init) {
    for (itr=Bin2Elm_.begin(); itr!=Bin2Elm_.end(); ++itr) {
      assert(used.find(itr->first)==used.end());
      used[itr->first] = false;
    }

    assert(Hm_.compatible(ctl));
    zmat.calculate(used, Bin2Elm_, Hm_, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);
    zmat_ctl.calculate(used, Bin2Elm_, ctl, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);
    //zmat.print();
    //zmat_ctl.print();
    for (i=0; i!=Bin2Elm_.size(); ++i) {
      for (j=0; j!=i; ++j) {
	double z = zmat(i,j) - zmat_ctl(i,j);
	zmat(i,j) = zmat(j,i) = z;
      }
    }
    //zmat.print();

    init = false;
  }

  bool last;
  map<uint,vector<uint> >::iterator itrdx, jtrdx;  

  progress("the initial pair of bins");
  last = true;
  itrdx = jtrdx = Bin2Elm_.begin();
  for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
    if (used[itr->first]) continue;
    for (jtr=Bin2Elm_.begin(), j=0; jtr!=itr; ++jtr, ++j) {
      if (used[jtr->first]) continue;
      if (zmat(i,j) > zmax) {
	zmax = zmat(i,j);
	itrdx = itr;
	jtrdx = jtr;
	last = false;
      }
    }
  }
  if (last) return; 
  const vector<uint>& ielm = itrdx->second;
  const vector<uint>& jelm = jtrdx->second;
  for (uint idx=0; idx!=ielm.size(); ++idx) 
    subset.push_back(Elm_, ielm[idx]);
  for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
    subset.push_back(Elm_, jelm[jdx]);
  used[itrdx->first] = true;
  used[jtrdx->first] = true;
  if (verbose) {
    cout << "z = " << zmax << endl;
    cout << "bin" << itrdx->first << endl;
    for (uint idx=0; idx!=ielm.size(); ++idx) 
      cout << Elm_.line(ielm[idx]) << flush;
    cout << "bin" << jtrdx->first << endl;
    for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
      cout << Elm_.line(jelm[jdx]) << flush;
  }

  for (uint k=2; k<Bin2Elm_.size(); ++k) {
    progress("the next bin");
    last = true;
    itrdx = jtrdx = Bin2Elm_.begin();

    ZscoreArray zarr;
    ZscoreArray zarr_ctl;
    zarr.calculate(used, subset, Bin2Elm_, Hm_, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);
    zarr_ctl.calculate(used, subset, Bin2Elm_, ctl, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);
    for (i=0; i!=Bin2Elm_.size(); ++i) {
      double z = zarr(i) - zarr_ctl(i);
      zarr(i) = z;
    }

    for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
      if (used[itr->first]) continue;
      if (zarr(i) >= zmax) {
	zmax = zarr(i);
	itrdx = itr;
	last = false;
      }
    }
    if (last) return;
    const vector<uint>& ielm = itrdx->second;
    for (uint idx=0; idx!=ielm.size(); ++idx) 
      subset.push_back(Elm_, ielm[idx]);
    used[itrdx->first] = true;
    if (verbose) {
      cout << "z = " << zmax << endl;
      cout << "bin" << itrdx->first << endl;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	cout << Elm_.line(ielm[idx]) << flush;
    }
  }

  return;
}

bool HypothesisTest::
look_interactive(uint nsmp, bool verbose, uint nthr) 
{
  int mode = -1;
  cout << "(1) look at a heatmap value for a position pair"          << endl
       << "(2) look at a z-score for a position pair"                << endl
       << "(3) look at heatmap values for one-to-all position pairs" << endl
       << "(4) look at z-scores for one-to-all position pairs"       << endl
       << "(5) output a submatrix of the heatmap"                    << endl
       << "(6) exit"                                                 << endl
       << "select a mode:"                                           << flush;
  cin  >> mode;
  cout << endl << flush;

  cout << "input format:" << endl 
       << "chr1:30000000 (chromosome-wise nucleotide coordinate)" << endl
       << "bin:30 (whole-genome bin coordinate)" << endl << flush;

  string label_f, label_s, buf_f, buf_s;
  uint pos_f, bin_f, pos_s, bin_s;
  if (mode==1) {
    cout << "input the 1st position:" << flush;
    cin  >> buf_f;
    if (! parse_position(buf_f, label_f, pos_f, bin_f)) return false; 
    format_position(buf_f, label_f, pos_f, bin_f);
    cout << "input the 2nd position:" << flush;
    cin  >> buf_s;
    if (! parse_position(buf_s, label_s, pos_s, bin_s)) return false; 
    format_position(buf_s, label_s, pos_s, bin_s);
    cout << endl << flush;

    cout << "heatmap value for a pair " << buf_f << " " << buf_s << endl;
    cout << Hm_(bin_f, bin_s) << endl << flush;
    cout << endl << flush;
    return true;
  }
  else if (mode==2) {
    cout << "input the 1st position:" << flush;
    cin  >> buf_f;
    if (! parse_position(buf_f, label_f, pos_f, bin_f)) return false; 
    format_position(buf_f, label_f, pos_f, bin_f);
    cout << "input the 2nd position:" << flush;
    cin  >> buf_s;
    if (! parse_position(buf_s, label_s, pos_s, bin_s)) return false; 
    format_position(buf_s, label_s, pos_s, bin_s);
    cout << endl << flush;

    Bed tmp;
    tmp.push_back(label_f, pos_f, pos_f+1, "1st", 0, "+");
    tmp.push_back(label_s, pos_s, pos_s+1, "2nd", 0, "+");
    HypothesisTest htest(Ncnk_, seed_, false, Hm_, tmp, Gap_);
    double stat, p, z;
    htest.test_stat();
    htest.null_dist(nsmp, false);
    htest.get_stat(stat, p, z);

    cout << "z-score for a pair " << buf_f << " " << buf_s << endl;
    cout << z << endl << flush;
    cout << endl << flush;
    return true;
  }
  else if (mode==3) {
    cout << "input the position:" << flush;
    cin  >> buf_f;
    if (! parse_position(buf_f, label_f, pos_f, bin_f)) return false; 
    format_position(buf_f, label_f, pos_f, bin_f);
    cout << endl << flush;

    cout << "heatmap values for pairs " << buf_f << " to all" << endl;
    map<uint,vector<uint> >::iterator itr;
    for (itr=Bin2Elm_.begin(); itr!=Bin2Elm_.end(); ++itr) {
      const vector<uint>& ielm = itr->second;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	cout << Hm_(bin_f, itr->first) << "\t" << Elm_.line(ielm[idx]) << flush;
    }
    cout << endl << flush;
    return true;
  }
  else if (mode==4) {
    cout << "input the position:" << flush;
    cin  >> buf_f;
    if (! parse_position(buf_f, label_f, pos_f, bin_f)) return false; 
    format_position(buf_f, label_f, pos_f, bin_f);
    cout << endl << flush;

    Bed tmp;
    ZscoreArray zarr;
    map<uint,bool> used;
    // make sure incrementing them coordinately
    map<uint,vector<uint> >::iterator itr;
    uint i;

    tmp.push_back(label_f, pos_f, pos_f+1, "1st", 0, "+");
    for (itr=Bin2Elm_.begin(); itr!=Bin2Elm_.end(); ++itr) {
      assert(used.find(itr->first)==used.end());
      used[itr->first] = false;
    }
    zarr.calculate(used, tmp, Bin2Elm_, Hm_, Elm_, Gap_, Ncnk_, seed_, nsmp, verbose, nthr);

    cout << "z-scores for pairs " << buf_f << " to all" << endl;
    for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
      const vector<uint>& ielm = itr->second;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	cout << zarr(i) << "\t" << Elm_.line(ielm[idx]) << flush;
    }
    cout << endl << flush;
    return true;
  }
  else if (mode==5) {
    string filename;
    uint bin_st_row, bin_sp_row, bin_st_col, bin_sp_col;
    cout << "input the start position of rows:" << flush;
    cin  >> buf_f;
    if (! parse_position(buf_f, label_f, pos_f, bin_st_row)) return false; 
    cout << "input the stop position of rows:" << flush;
    cin  >> buf_s;
    if (! parse_position(buf_s, label_s, pos_s, bin_sp_row)) return false; 
    cout << "input the start position of cols:" << flush;
    cin  >> buf_f;
    if (! parse_position(buf_f, label_f, pos_f, bin_st_col)) return false; 
    cout << "input the stop position of cols:" << flush;
    cin  >> buf_s;
    if (! parse_position(buf_s, label_s, pos_s, bin_sp_col)) return false; 

    if (! (bin_st_row < bin_sp_row && bin_st_col < bin_sp_col && 
	   bin_sp_row <= Hm_.nbin() && bin_sp_col <= Hm_.nbin()) ) return false;

    cout << "input the output filename:" << flush;
    cin  >> filename;
    cout << endl << flush;

    return Hm_.output_val(filename, bin_st_row, bin_sp_row, bin_st_col, bin_sp_col);
  }
  else {
    return false;
  }

  return false; // never come here
}

bool HypothesisTest::
parse_position(const string buf, 
	       string& label, uint& pos, uint& bin)
{
  string::size_type idx = buf.find(":");
  if (idx==string::npos || idx==buf.size()-1) return false;

  string head = buf.substr(0, idx);
  stringstream ss(buf.substr(idx+1));
  if (head=="bin") {
    ss >> bin;
    return Hm_.bin2pos(bin, label, pos);
  }
  else {
    label = head;
    ss >> pos;
    return Hm_.pos2bin(label, pos, bin);
  }

  return false; // never come here
}

bool HypothesisTest::
format_position(string& buf, 
		const string label, const uint pos, const uint bin)
{
  stringstream ss;
  ss << label << ":" << pos << "(" << bin << ")";
  buf = ss.str();

  return true;
}

void HypothesisTest::
calc_pval()
{
  uint nsmp = NullDist_.size();

  Pval_ = 1.0 / (nsmp + 1);
  for (uint i=0; i!=nsmp; ++i) {
    if (NullDist_[i] > TestStat_) {
      Pval_ = (double) (nsmp - i + 1) / (nsmp + 1);
      break;
    }
  }
}

void HypothesisTest::
calc_zval()
{
  uint nsmp = NullDist_.size();
  double mean = 0.0;
  double sd = 0.0;

  for (uint i=0; i!=nsmp; ++i) 
    mean += NullDist_[i];
  mean /= nsmp;

  for (uint i=0; i!=nsmp; ++i) {
    double d = NullDist_[i] - mean;
    sd += d * d;
  }
  sd = sqrt(sd / nsmp);

  if (sd == 0.0) {
    if (TestStat_ > mean) Zval_ = ZvalMax_;
    else if (TestStat_ < mean) Zval_ = - ZvalMax_;
    else Zval_ = 0.0;
  }
  else {
    Zval_ = (TestStat_ - mean) / sd;
  }
}


void HypothesisTest::
search_subset_single(Bed& subset, double& zmax, uint nsmp, bool verbose)
{
  static bool init = true;
  static map<uint,bool> used;
  static vector<vector<double> > zmat;

  // make sure incrementing them coordinately
  map<uint,vector<uint> >::iterator itr, jtr;
  uint i, j;

  if (init) {
    /*
    cout << "Bin2Elm_:" << endl << flush;
    for (itr=Bin2Elm_.begin(); itr!=Bin2Elm_.end(); ++itr) {
      cout << "bin" << itr->first << endl << flush;
      const vector<uint>& ielm = itr->second;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	cout << Elm_.line(ielm[idx]) << flush;
    }
    */
    for (itr=Bin2Elm_.begin(); itr!=Bin2Elm_.end(); ++itr) {
      assert(used.find(itr->first)==used.end());
      used[itr->first] = false;
    }
    zmat.resize(Bin2Elm_.size(), vector<double>(Bin2Elm_.size(), - INF));
    for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
      for (jtr=Bin2Elm_.begin(), j=0; jtr!=itr; ++jtr, ++j) {
	if (verbose) cout << itr->first << "," << jtr->first << " " << flush;
	Bed tmp;
	const vector<uint>& ielm = itr->second;
	const vector<uint>& jelm = jtr->second;
	for (uint idx=0; idx!=ielm.size(); ++idx) 
	  tmp.push_back(Elm_, ielm[idx]);
	for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
	  tmp.push_back(Elm_, jelm[jdx]);
	HypothesisTest htest(Ncnk_, seed_, false, Hm_, tmp, Gap_);
	double stat, p, z;
	htest.test_stat();
	htest.null_dist(nsmp, false);
	htest.get_stat(stat, p, z);
	zmat[i][j] = zmat[j][i] = z;
      }
    }
    if (verbose) cout << endl << flush;
    /*
    cout << "Zscore Matrix:" << endl << flush;
    for (i=0; i!=Bin2Elm_.size(); ++i) {
      for (j=0; j!=Bin2Elm_.size()-1; ++j) 
	cout << zmat[i][j] << "\t" << flush;      
      cout << zmat[i][Bin2Elm_.size()-1] << endl << flush;
    }
    */
    init = false;
  }

  bool last;
  map<uint,vector<uint> >::iterator itrdx, jtrdx;  

  progress("the initial pair of bins");
  last = true;
  itrdx = jtrdx = Bin2Elm_.begin();
  for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
    if (used[itr->first]) continue;
    for (jtr=Bin2Elm_.begin(), j=0; jtr!=itr; ++jtr, ++j) {
      if (used[jtr->first]) continue;
      if (zmat[i][j] > zmax) {
	zmax = zmat[i][j];
	itrdx = itr;
	jtrdx = jtr;
	last = false;
      }
    }
  }
  if (last) return; 
  const vector<uint>& ielm = itrdx->second;
  const vector<uint>& jelm = jtrdx->second;
  for (uint idx=0; idx!=ielm.size(); ++idx) 
    subset.push_back(Elm_, ielm[idx]);
  for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
    subset.push_back(Elm_, jelm[jdx]);
  used[itrdx->first] = true;
  used[jtrdx->first] = true;
  if (verbose) {
    cout << "z = " << zmax << endl;
    cout << "bin" << itrdx->first << endl;
    for (uint idx=0; idx!=ielm.size(); ++idx) 
      cout << Elm_.line(ielm[idx]) << flush;
    cout << "bin" << jtrdx->first << endl;
    for (uint jdx=0; jdx!=jelm.size(); ++jdx) 
      cout << Elm_.line(jelm[jdx]) << flush;
  }

  for (uint k=2; k<Bin2Elm_.size(); ++k) {
    progress("the next bin");
    last = true;
    itrdx = jtrdx = Bin2Elm_.begin();
    for (itr=Bin2Elm_.begin(), i=0; itr!=Bin2Elm_.end(); ++itr, ++i) {
      if (used[itr->first]) continue;
      if (verbose) cout << "." << flush;
      Bed tmp = subset;
      const vector<uint>& ielm = itr->second;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	tmp.push_back(Elm_, ielm[idx]);
      HypothesisTest htest(Ncnk_, seed_, false, Hm_, tmp, Gap_);
      double stat, p, z;
      htest.test_stat();
      htest.null_dist(nsmp, false);
      htest.get_stat(stat, p, z);
      if (z >= zmax) {
	zmax = z;
	itrdx = itr;
	last = false;
      }
    }
    if (verbose) cout << endl << flush;
    if (last) return;
    const vector<uint>& ielm = itrdx->second;
    for (uint idx=0; idx!=ielm.size(); ++idx) 
      subset.push_back(Elm_, ielm[idx]);
    used[itrdx->first] = true;
    if (verbose) {
      cout << "z = " << zmax << endl;
      cout << "bin" << itrdx->first << endl;
      for (uint idx=0; idx!=ielm.size(); ++idx) 
	cout << Elm_.line(ielm[idx]) << flush;
    }
  }
}
