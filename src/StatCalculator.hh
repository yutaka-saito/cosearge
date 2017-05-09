#ifndef __INC_STAT_CALCULATOR_HH__
#define __INC_STAT_CALCULATOR_HH__

#include <vector>
#include <map>
#include <cassert>

#include "Utility.hh"
#include "Heatmap.hh"
#include "Bed.hh"

class ZscoreMatrix 
{
public:
  typedef double ValueType;

public:
  ZscoreMatrix() : Nrow_(0), Ncol_(0) {}
  ~ZscoreMatrix() {}

  void resize(uint n, uint m) {
    Zscore_.clear();
    Zscore_.resize(n, std::vector<ValueType>(m, DefVal_));
    Nrow_ = n;
    Ncol_ = m;
  }

  double& operator()(uint i, uint j) {
    assert(i < Nrow_ && j < Ncol_);
    return Zscore_[i][j];
  }
  const double& operator()(uint i, uint j) const {
    assert(i < Nrow_ && j < Ncol_);
    return Zscore_[i][j];
  }

  uint nrow() const {
    return Nrow_;
  }
  uint ncol() const {
    return Ncol_;
  }

  void calculate(const std::map<uint,bool>& used, 
		 const std::map<uint,std::vector<uint> >& bin2elm, 
		 const Heatmap& hm, const Bed& elm, const Bed& gap, 
		 uint ncnk, uint seed, uint nsmp, bool verbose, uint nthr);

  void print();

private:
  std::vector<std::vector<ValueType> > Zscore_;
  uint Nrow_;
  uint Ncol_;

private:
  static const ValueType DefVal_ = - INF;
};

class ZscoreArray
{
public:
  typedef double ValueType;

public:
  ZscoreArray() : Size_(0) {}
  ~ZscoreArray() {}

  void resize(uint n) {
    Zscore_.clear();
    Zscore_.resize(n, DefVal_);
    Size_ = n;
  }

  double& operator()(uint i) {
    assert(i < Size_);
    return Zscore_[i];
  }
  const double& operator()(uint i) const {
    assert(i < Size_);
    return Zscore_[i];
  }

  uint size() const {
    return Size_;
  }

  void calculate(const std::map<uint,bool>& used, const Bed& preset, 
		 const std::map<uint,std::vector<uint> >& bin2elm, 
		 const Heatmap& hm, const Bed& elm, const Bed& gap, 
		 uint ncnk, uint seed, uint nsmp, bool verbose, uint nthr);

  void print();

private:
  std::vector<ValueType> Zscore_;
  uint Size_;

private:
  static const ValueType DefVal_ = - INF;
};

#endif 
