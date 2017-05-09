#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Heatmap.hh"

using namespace std;

// fixed file names in a heatmap directory
const string MatrixFile = "heatmap";
const string ChrIdxFile = "chromosomeIndex";
const string PosIdxFile = "positionIndex";
const string ChrLabelFile = "genomeIdxToLabel";
const string ChrStartFile = "chromosomeStarts";
const string NbinFile = "binNumber";
const string ResolFile = "resolution";
const string StatFile = "statistics";

string 
make_filename(const string& dir, const string& file, const string& suffix) 
{
  if (suffix == "") return dir + "/" + file;
  else return dir + "/" + file + "_" + suffix;
}

bool Heatmap::
load(const string& in_dir, const string& suffix) 
{
  string filename, line;
  ifstream ifs;
  stringstream ss;
  uint n;

  // binNumber
  filename = in_dir + "/" + NbinFile;
  ifs.open(filename.c_str());
  if (!ifs) return false;
  getline(ifs, line); 
  ss << line;
  ss >> Nbin_;
  clear_ss(ss);
  ifs.close();
  assert(Nbin_ > 1);

  // resolution
  filename = in_dir + "/" + ResolFile;
  ifs.open(filename.c_str());
  if (!ifs) return false;
  getline(ifs, line); 
  ss << line;
  ss >> Resol_;
  clear_ss(ss);
  ifs.close();

  // chromosomeIndex
  ChrIdx_.clear();
  ChrIdx_.resize(Nbin_);
  n = 0;
  filename = in_dir + "/" + ChrIdxFile;
  ifs.open(filename.c_str());
  if (!ifs) return false;
  while (getline(ifs, line) != NULL) {
    ss << line;
    ChrIdx_[n++] = ss2a2i(ss);
    clear_ss(ss);
  }
  ifs.close();
  assert(n == Nbin_);

  // positionIndex
  PosIdx_.clear();
  PosIdx_.resize(Nbin_);
  n = 0;
  filename = in_dir + "/" + PosIdxFile;
  ifs.open(filename.c_str());
  if (!ifs) return false;
  while (getline(ifs, line) != NULL) {
    ss << line;
    PosIdx_[n++] = ss2a2i(ss);
    clear_ss(ss);
  }
  ifs.close();
  assert(n == Nbin_);

  // genomeIdxToLabel
  ChrLabel_.clear();
  filename = in_dir + "/" + ChrLabelFile;
  ifs.open(filename.c_str());
  if (!ifs) return false;
  getline(ifs, line); 
  vector<uint> qts;
  for (uint i=0; i<line.size(); ++i) {
    if (line.find("'", i) == string::npos) break;
    i = line.find("'", i);
    qts.push_back(i);
  }
  assert(qts.size() % 2 == 0);
  for (uint i=0; i<qts.size(); i=i+2) 
    ChrLabel_.push_back("chr" + line.substr(qts[i]+1, qts[i+1]-qts[i]-1));
  Nchr_ = ChrLabel_.size();
  ifs.close();

  // chromosomeStarts
  ChrStart_.clear();
  ChrStart_.resize(Nchr_);
  n = 0;
  filename = in_dir + "/" + ChrStartFile;
  ifs.open(filename.c_str());
  if (!ifs) return false;
  while (getline(ifs, line) != NULL) {
    ss << line;
    ChrStart_[n++] = ss2a2i(ss);
    clear_ss(ss);
  }
  ifs.close();
  assert(n == Nchr_);

  // heatmap
  reset_matrix(Nbin_, Nbin_);
  n = 0;
  filename = make_filename(in_dir, MatrixFile, suffix);
  ifs.open(filename.c_str());
  if (!ifs) return false;
  while (getline(ifs, line) != NULL) {
    vector<ValueType> row(Nbin_);
    ValueType v = 0.0;
    uint m = 0;
    ss << line;
    while (true) {
      v = ss2a2f(ss);
      if (!ss) break;
      row[m++] = v;
    }
    assert(m == Nbin_);
    set_row(row, n++);
    clear_ss(ss);
  }
  ifs.close();
  assert(n == Nbin_);

  //print();
  return true;
}

bool Heatmap::
output_stat(const string& out_dir, const string& suffix) const
{
  string filename = make_filename(out_dir, StatFile, suffix);
  ofstream ofs(filename.c_str());
  if (!ofs) return false;

  ofs << "chromosome-wise statistics:" << endl;
  ofs << "IntraMean:" << endl;
  for (uint chri=0; chri!=Nchr_; ++chri) {
    ofs << "chromosome " << chri << endl;
    for (uint i=0; i!=Nbin_-1; ++i) 
      ofs << IntraMean_[chri][i] << endl;
  }
  ofs << "IntraStdev:" << endl;
  for (uint chri=0; chri!=Nchr_; ++chri) {
    ofs << "chromosome " << chri << endl;
    for (uint i=0; i!=Nbin_-1; ++i) 
      ofs << IntraStdev_[chri][i] << endl;
    ofs << endl;
  }
  ofs << "InterMean:" << endl;
  for (uint chri=0; chri!=Nchr_; ++chri) {
    for (uint chrj=0; chrj!=Nchr_-1; ++chrj) 
      ofs << InterMean_[chri][chrj] << "\t";
    ofs << InterMean_[chri][Nchr_-1] << endl;
  }
  ofs << "InterStdev:" << endl;
  for (uint chri=0; chri!=Nchr_; ++chri) {
    for (uint chrj=0; chrj!=Nchr_-1; ++chrj) 
      ofs << InterStdev_[chri][chrj] << "\t";
    ofs << InterStdev_[chri][Nchr_-1] << endl;
  }

  ofs << "all-chromosome statistics:" << endl;
  ofs << "IntraMean:" << endl;
  for (uint i=0; i!=Nbin_-1; ++i) 
    ofs << AllIntraMean_[i] << endl;
  ofs << endl;
  ofs << "IntraStdev:" << endl;
    for (uint i=0; i!=Nbin_-1; ++i) 
      ofs << AllIntraStdev_[i] << endl;
    ofs << endl;
  ofs << "InterMean:" << endl << AllInterMean_ << endl;
  ofs << "InterStdev:" << endl << AllInterStdev_ << endl;

  return true;
}

bool Heatmap::
output_val(const string& out_dir, const string& suffix) const 
{
  string filename = make_filename(out_dir, MatrixFile, suffix);

  return output_val(filename, 0, Nbin_, 0, Nbin_);
}

bool Heatmap::
output_val(const string& filename, 
	   uint st_row, uint sp_row, uint st_col, uint sp_col) const  
{
  const uint prec = 5;
  
  ofstream ofs(filename.c_str());
  if (!ofs) return false;

  for (uint i=st_row; i<sp_row; ++i) {
    vector<ValueType> row;
    get_row(row, i);
    for (uint j=st_col; j<sp_col; ++j) {
      ofs << fixed << setprecision(prec) << scientific << row[j];
      if (j == sp_col-1) ofs << endl;
      else ofs << "\t";
    }
  }

  return true;
}

bool Heatmap::
output_val_chrsep(const string& out_dir, const string& suffix) 
{
  const string na = "na";
  const uint prec = 5;

  for (uint i=0; i<Nbin_; ++i) {
    if (is_null(i)) continue; 
    ValueType max = DefVal_;
    vector<ValueType> row;
    get_row(row, i);
    for (uint j=0; j<Nbin_; ++j) {
      if (is_null(j)) continue; 
      if (row[j] > max) max = row[j];
    }
    row[i] = max;
    set_row(row, i);
  }

  for (uint chri=0; chri<Nchr_; ++chri) {
    for (uint chrj=chri; chrj<Nchr_; ++chrj) {
      string sfx = suffix + "_" + ChrLabel_[chri] + "_" + ChrLabel_[chrj];
      string filename = make_filename(out_dir, MatrixFile, sfx);
      ofstream ofs(filename.c_str());
      if (!ofs) return false;

      for (uint i=chr_start(chri); i<chr_stop(chri); ++i) {
	vector<ValueType> row;
	get_row(row, i);
	for (uint j=chr_start(chrj); j<chr_stop(chrj); ++j) {
	  if (is_null(i) || is_null(j)) ofs << na;
	  else ofs << fixed << setprecision(prec) << scientific << row[j];
	  if (j == chr_stop(chrj)-1) ofs << endl;
	  else ofs << "\t";
	}
      }
    }
  }

  return true;
}

void Heatmap::
normalize(float aggr, bool chrwise, bool verbose) 
{
  // IntraMean_[chri][i], IntraStdev_[chri][i], 
  // AllIntraMean_[i] and AllIntraStdev_[i] are calculated by 
  // aggregating neighbor bins when the number of observed interactions
  // IntraCnt[chri][i] < AggrThsh_ * IntraTotalCnt[chri],  
  // AllIntraCnt[i] < AggrThsh_ * AllIntraTotalCnt.
  const float AggrThsh = aggr;

  IntraMean_.clear();
  IntraMean_.resize(Nchr_, vector<ValueType>(Nbin_-1, 0.0));
  IntraStdev_.clear();
  IntraStdev_.resize(Nchr_, vector<ValueType>(Nbin_-1, 0.0));
  InterMean_.clear();
  InterMean_.resize(Nchr_, vector<ValueType>(Nchr_, 0.0));
  InterStdev_.clear();
  InterStdev_.resize(Nchr_, vector<ValueType>(Nchr_, 0.0));
  AllIntraMean_.clear();
  AllIntraMean_.resize(Nbin_-1, 0.0);
  AllIntraStdev_.clear();
  AllIntraStdev_.resize(Nbin_-1, 0.0);
  AllInterMean_ = 0.0;
  AllInterStdev_ = 0.0;

  // chromosome-wise statistics
  vector<vector<double> > IntraSum(Nchr_, vector<double>(Nbin_-1, 0.0)); 
  vector<vector<double> > IntraSqd(Nchr_, vector<double>(Nbin_-1, 0.0)); 
  vector<vector<uint> > IntraCnt(Nchr_, vector<uint>(Nbin_-1, 0)); 
  vector<uint> IntraTotalCnt(Nchr_, 0); 
  vector<vector<double> > InterSum(Nchr_, vector<double>(Nchr_, 0.0)); 
  vector<vector<double> > InterSqd(Nchr_, vector<double>(Nchr_, 0.0)); 
  vector<vector<uint> > InterTotalCnt(Nchr_, vector<uint>(Nchr_, 0)); 
  // all-chromosome statistics
  vector<double> AllIntraSum(Nbin_-1, 0.0);
  vector<double> AllIntraSqd(Nbin_-1, 0.0);
  vector<uint> AllIntraCnt(Nbin_-1, 0);
  uint AllIntraTotalCnt = 0;
  double AllInterSum = 0.0;
  double AllInterSqd = 0.0;
  uint AllInterTotalCnt = 0;

  progress("calculating mean");
  for (uint i=0; i<Nbin_; ++i) {
    if (verbose) cout << "." << flush;
    if (is_null(i)) continue;
    vector<ValueType> row;
    get_row(row, i);
    for (uint j=i; j<Nbin_; ++j) {
      if (is_null(j)) continue;
      if (ChrIdx_[i] == ChrIdx_[j]) {
	IntraSum[ChrIdx_[i]][j-i] += row[j];
	IntraCnt[ChrIdx_[i]][j-i]++;
	IntraTotalCnt[ChrIdx_[i]]++;
	AllIntraSum[j-i] += row[j];
	AllIntraCnt[j-i]++;
	AllIntraTotalCnt++;
      }
      else {
	InterSum[ChrIdx_[i]][ChrIdx_[j]] += row[j];
	InterSum[ChrIdx_[j]][ChrIdx_[i]] = InterSum[ChrIdx_[i]][ChrIdx_[j]];
	InterTotalCnt[ChrIdx_[i]][ChrIdx_[j]]++;
	InterTotalCnt[ChrIdx_[j]][ChrIdx_[i]] = InterTotalCnt[ChrIdx_[i]][ChrIdx_[j]];
	AllInterSum += row[j];
	AllInterTotalCnt++;
      }
    }
  }
  if (verbose) cout << endl << flush;
  // aggregating neighbor bins for intra-chromosomal interactions
  for (uint chri=0; chri<Nchr_; ++chri) {
    for (uint i=0; i<Nbin_-1; ++i) {
      double mean = 0.0;
      double sum = 0.0;
      uint cnt = 0;
      for (uint j=i; j<Nbin_-1; ++j) {
	uint jdx = Nbin_-1-j-1;
	sum += IntraSum[chri][jdx];
	cnt += IntraCnt[chri][jdx];
	if (cnt > AggrThsh * IntraTotalCnt[chri] || jdx == 0) {
	  mean = sum / cnt;
	  for (uint k=i; k<=j; ++k) {
	    uint kdx = Nbin_-1-k-1;
	    IntraMean_[chri][kdx] = mean;
	  }
	  if (verbose) cout << "aggregated bins for the chromosome " << chri 
			    << " to " << jdx << endl << flush;
	  i = j;
	  break;
	}
      }
    }
  }
  for (uint i=0; i<Nbin_-1; ++i) {
    ValueType mean = 0.0;
    ValueType sum = 0.0;
    uint cnt = 0;
    for (uint j=i; j<Nbin_-1; ++j) {
      uint jdx = Nbin_-1-j-1;
      sum += AllIntraSum[jdx];
      cnt += AllIntraCnt[jdx];
      if (cnt > AggrThsh * AllIntraTotalCnt || jdx == 0) {
	mean = sum / cnt;
	for (uint k=i; k<=j; ++k) {
	  uint kdx = Nbin_-1-k-1;
	  AllIntraMean_[kdx] = mean;
	}
	if (verbose) cout << "aggregated bins for all chromosomes"  
			  << " to " << jdx << endl << flush;
	i = j;
	break;
      }
    }
  }
  // for inter-chromosomal interactions
  for (uint chri=0; chri<Nchr_; ++chri) {
    for (uint chrj=chri+1; chrj<Nchr_; ++chrj) {
      InterMean_[chri][chrj] = InterSum[chri][chrj] / InterTotalCnt[chri][chrj];
      InterMean_[chrj][chri] = InterMean_[chri][chrj];
    }
  }
  AllInterMean_ = AllInterSum / AllInterTotalCnt;

  progress("calculating stdev");
  for (uint i=0; i<Nbin_; ++i) {    
    if (verbose) cout << "." << flush;
    if (is_null(i)) continue;
    vector<ValueType> row;
    get_row(row, i);
    for (uint j=i; j<Nbin_; ++j) {
      if (is_null(j)) continue;
      if (ChrIdx_[i] == ChrIdx_[j]) {
	double d = row[j] - IntraMean_[ChrIdx_[i]][j-i];
	IntraSqd[ChrIdx_[i]][j-i] += d * d;
	d = row[j] - AllIntraMean_[j-i];
	AllIntraSqd[j-i] += d * d;
      }
      else {
	double d = row[j] - InterMean_[ChrIdx_[i]][ChrIdx_[j]];
	InterSqd[ChrIdx_[i]][ChrIdx_[j]] += d * d;
	InterSqd[ChrIdx_[j]][ChrIdx_[i]] = InterSqd[ChrIdx_[i]][ChrIdx_[j]];
	d = row[j] - AllInterMean_;
	AllInterSqd += d * d;
      }
    }
  }
  if (verbose) cout << endl << flush;
  // aggregating neighbor bins for intra-chromosomal interactions
  for (uint chri=0; chri<Nchr_; ++chri) {
    for (uint i=0; i<Nbin_-1; ++i) {
      double stdev = 0.0;
      double sqd = 0.0;
      uint cnt = 0;
      for (uint j=i; j<Nbin_-1; ++j) {
	uint jdx = Nbin_-1-j-1;
	sqd += IntraSqd[chri][jdx];
	cnt += IntraCnt[chri][jdx];
	if (cnt > AggrThsh * IntraTotalCnt[chri] || jdx == 0) {
	  stdev = sqrt(sqd / cnt);
	  for (uint k=i; k<=j; ++k) {
	    uint kdx = Nbin_-1-k-1;
	    IntraStdev_[chri][kdx] = stdev;
	  }
	  if (verbose) cout << "aggregated bins for the chromosome " << chri 
			    << " to " << jdx << endl << flush;
	  i = j;
	  break;
	}
      }
    }
  }
  for (uint i=0; i<Nbin_-1; ++i) {
    ValueType stdev = 0.0;
    ValueType sqd = 0.0;
    uint cnt = 0;
    for (uint j=i; j<Nbin_-1; ++j) {
      uint jdx = Nbin_-1-j-1;
      sqd += AllIntraSqd[jdx];
      cnt += AllIntraCnt[jdx];
      if (cnt > AggrThsh * AllIntraTotalCnt || jdx == 0) {
	stdev = sqrt(sqd / cnt);
	for (uint k=i; k<=j; ++k) {
	  uint kdx = Nbin_-1-k-1;
	  AllIntraStdev_[kdx] = stdev;
	}
	if (verbose) cout << "aggregated bins for all chromosomes"
			  << " to " << jdx << endl << flush;
	i = j;
	break;
      }
    }
  }
  // for inter-chromosomal interactions
  for (uint chri=0; chri<Nchr_; ++chri) {
    for (uint chrj=chri+1; chrj<Nchr_; ++chrj) {
      InterStdev_[chri][chrj] = sqrt(InterSqd[chri][chrj] / InterTotalCnt[chri][chrj]);
      InterStdev_[chrj][chri] = InterStdev_[chri][chrj];
    }
  }
  AllInterStdev_ = sqrt(AllInterSqd / AllInterTotalCnt);

  /*
  for (uint chri=0; chri<Nchr_; ++chri) {
    for (uint i=0; i<Nbin_-1; ++i) {
      ValueType stdev = 0.0;
      ValueType sqd = 0.0;
      uint cnt = 0;
      for (uint j=i; j<Nbin_-1; ++j) {
	uint jdx = Nbin_-1-j-1;
	sqd += Sqd[jdx];
	cnt += Cnt[jdx];
	if (cnt > AggrThsh * IntraCnt || jdx == 0) {
	  stdev = sqrt(sqd / cnt);
	  for (uint k=i; k<=j; ++k) {
	    uint kdx = Nbin_-1-k-1;
	    Stdev_[kdx] = stdev;
	  }
	  if (verbose) cout << "aggregated bins to " << jdx << endl << flush;
	  i = j;
	  break;
	}
      }
    }
  }
  InterStdev_ = sqrt(InterSqd / InterCnt);
  */

  if (chrwise) { 
    progress("normalizing with chromosome-wise statistics");
    for (uint i=0; i<Nbin_; ++i) {
      if (verbose) cout << "." << flush;
      if (is_null(i)) continue;
      vector<ValueType> new_row(Nbin_, ValueType(DefVal_));
      vector<ValueType> row;
      get_row(row, i);
      for (uint j=0; j<Nbin_; ++j) {
	if (is_null(j)) continue;
	uint d = uabs(i, j);
	if (ChrIdx_[i] == ChrIdx_[j]) {
	  if (d < 2) 
	    new_row[j] = 0.0;
	  else 
	    new_row[j] = (row[j] - IntraMean_[ChrIdx_[i]][d]) 
	      / IntraStdev_[ChrIdx_[i]][d];
	}
	else {
	  new_row[j] = (row[j] - InterMean_[ChrIdx_[i]][ChrIdx_[j]]) 
	    / InterStdev_[ChrIdx_[i]][ChrIdx_[j]];
	}
      }
      set_row(new_row, i);
    }
    if (verbose) cout << endl << flush;
  }
  else { // chrwise == false
    progress("normalizing with all-chromosome statistics");
    for (uint i=0; i<Nbin_; ++i) {
      if (verbose) cout << "." << flush;
      if (is_null(i)) continue;
      vector<ValueType> new_row(Nbin_, ValueType(DefVal_));
      vector<ValueType> row;
      get_row(row, i);
      for (uint j=0; j<Nbin_; ++j) {
	if (is_null(j)) continue;
	uint d = uabs(i, j);
	if (ChrIdx_[i] == ChrIdx_[j]) {
	  if (d < 2) 
	    new_row[j] = 0.0;
	  else 
	    new_row[j] = (row[j] - AllIntraMean_[d]) / AllIntraStdev_[d];
	}
	else {
	  new_row[j] = (row[j] - AllInterMean_) / AllInterStdev_;
	}
      }
      set_row(new_row, i);
    }
    if (verbose) cout << endl << flush;
  }

}

void Heatmap::
logarithm() 
{
  for (uint i=0; i!=Nbin_; ++i) {
    vector<ValueType> new_row(Nbin_, DefVal_);
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j!=Nbin_; ++j) 
      new_row[j] = (my_row[j] > 0.0) ? log(my_row[j]) : - INF;
    set_row(new_row, i);
  }
}

void Heatmap::
scale() 
{
  double sum = 0.0;
  for (uint i=0; i!=Nbin_; ++i) {
    if (is_null(i)) continue;
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j!=Nbin_; ++j) 
      sum += my_row[j];
    break;
  }
  assert(sum != 0.0);

  for (uint i=0; i!=Nbin_; ++i) {
    vector<ValueType> new_row(Nbin_, DefVal_);
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j!=Nbin_; ++j) 
      new_row[j] = my_row[j] / sum;
    set_row(new_row, i);
  }
}

void Heatmap::
scaletj() 
{
  double numer = 0.0;
  vector<double> denom(Nbin_, 0.0);

  for (uint i=0; i<Nbin_; ++i) {
    if (is_null(i)) continue;
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j<Nbin_; ++j) 
      if (i != j) numer += my_row[j];
  }
  numer /= 2;
  assert(numer != 0.0);

  for (uint i=0; i<Nbin_; ++i) {
    if (is_null(i)) continue;
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j<Nbin_; ++j) 
      //if (i != j) denom[i] += my_row[j];
      denom[i] += my_row[j];
  }

  /*
  for (uint i=0; i<Nbin_; ++i) {
    if (is_null(i)) continue;
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j<Nbin_; ++j) 
      if (i == j) numer += 2 * my_row[j];
      else numer += my_row[j];
  }
  numer /= 2;
  assert(numer != 0.0);

  for (uint i=0; i<Nbin_; ++i) {
    if (is_null(i)) continue;
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j<Nbin_; ++j) 
      denom[i] += my_row[j];
  }
  */

  for (uint i=0; i!=Nbin_; ++i) {
    vector<ValueType> new_row(Nbin_, DefVal_);
    vector<ValueType> my_row;
    get_row(my_row, i);
    for (uint j=0; j!=Nbin_; ++j) 
      if (denom[i] != 0.0 && denom[j] != 0.0)
	new_row[j] = (my_row[j] * numer) / (denom[i] * denom[j]);
    set_row(new_row, i);
  }
}

void Heatmap::
subtract(const Heatmap& hm) 
{
  assert(compatible(hm));

  for (uint i=0; i!=Nbin_; ++i) {
    vector<ValueType> new_row(Nbin_, DefVal_);
    vector<ValueType> my_row, row;
    get_row(my_row, i);
    hm.get_row(row, i);
    for (uint j=0; j!=Nbin_; ++j) 
      new_row[j] = my_row[j] - row[j];
    set_row(new_row, i);
  }
}

void Heatmap::
divide(const Heatmap& hm) 
{
  assert(compatible(hm));

  for (uint i=0; i!=Nbin_; ++i) {
    vector<ValueType> new_row(Nbin_, DefVal_);
    vector<ValueType> my_row, row;
    get_row(my_row, i);
    hm.get_row(row, i);
    for (uint j=0; j!=Nbin_; ++j) {
      if (my_row[j]==0.0 && row[j]==0.0) {
	new_row[j] = 1.0;
      }
      else if (my_row[j]==0.0) {
	new_row[j] = 0.0;
      }
      else if (row[j]==0.0) {
	new_row[j] = (my_row[j] > 0.0) ? INF : - INF;
      }
      else {
	new_row[j] = my_row[j] / row[j];
      }
    }
    set_row(new_row, i);
  }
}

bool Heatmap::
compatible(const Heatmap& hm) const
{
  if (Nbin_ != hm.nbin() || Nchr_ != hm.nchr() || Resol_ != hm.resol())
    return false;

  for (uint i=0; i!=Nchr_; ++i) 
    if (ChrLabel_[i] != hm.chr_label(i))
      return false;

  for (uint i=0; i!=Nchr_; ++i) 
    if (ChrStart_[i] != hm.chr_start(i))
      return false;

  return true;
}

void Heatmap::
print() 
{
  cout << "print Heatmap class:" << endl << endl << flush;

  cout << "Nbin:" << endl << Nbin_ << endl << flush;  

  cout << "Resol:" << endl << Resol_ << endl << flush;  

  cout << "ChrIdx:" << endl;  
  for (uint i=0; i!=Nbin_; ++i) 
    cout << ChrIdx_[i] << endl << flush;
  cout << endl << flush;

  cout << "PosIdx:" << endl;  
  for (uint i=0; i!=Nbin_; ++i) 
    cout << PosIdx_[i] << endl << flush;
  cout << endl << flush;

  cout << "ChrLabel:" << endl;  
  for (uint i=0; i!=Nchr_; ++i) 
    cout << ChrLabel_[i] << endl << flush;
  cout << endl << flush;

  cout << "ChrStart:" << endl;  
  for (uint i=0; i!=Nchr_; ++i) 
    cout << ChrStart_[i] << endl << flush;
  cout << endl << flush;

  cout << "Matrix:" << endl;
  for (uint i=0; i!=Nbin_; ++i) {
    vector<ValueType> row;
    get_row(row, i);
    for (uint j=0; j!=Nbin_-1; ++j) 
      cout << row[j] << "\t" << flush;
    cout << row[Nbin_-1] << endl << flush;
  }
  cout << endl << flush;
}

bool Heatmap::
pos2bin(const string label, const uint pos, uint& chr, uint& bin) const 
{
  bin = Nbin_;
  chr = Nchr_;

  for (uint i=0; i<Nchr_; ++i) {
    if (label == ChrLabel_[i]) {
      chr = i;
      break;
    }
  }
  if (chr == Nchr_) return false; // unknown chromosome label

  uint start = ChrStart_[chr];
  for (uint i=start; i<Nbin_ && ChrIdx_[i]==ChrIdx_[start]; ++i) {
    if (pos < PosIdx_[i] + Resol_) {
      bin = i;
      break;
    }
  }
  if (bin == Nbin_) return false; // out-of-range coordinate

  assert(bin >= chr_start(chr) && bin < chr_stop(chr));
  return true;
}

bool Heatmap::
pos2bin(const string label, const uint pos, uint& bin) const 
{
  uint chr;
  return pos2bin(label, pos, chr, bin);
}

bool Heatmap::
bin2pos(const uint bin, string& label, uint& pos) const
{
  if (! (bin < Nbin_)) return false;
  label = ChrLabel_[ChrIdx_[bin]];
  pos = PosIdx_[bin];

  return true;
}
