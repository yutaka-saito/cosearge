#ifndef __INC_SPARSE_MATRIX_HH__
#define __INC_SPARSE_MATRIX_HH__

#include <cassert>
#include <vector>
#include <list>

struct Cell
{
  typedef float ValueType;

  Cell(const Cell& c) : col(c.col), val(c.val) {}
  Cell(uint i, ValueType v) : col(i), val(v) {}

  uint col;
  ValueType val;
};

struct SparseMatrixIterator
{
  typedef std::list<Cell>::iterator IteratorType;

  SparseMatrixIterator() {}

  SparseMatrixIterator(const SparseMatrixIterator& smi) : itr(smi.itr) {}
  SparseMatrixIterator(const IteratorType& smi_itr) : itr(smi_itr) {}

  SparseMatrixIterator operator++() { 
    ++itr; 
    return *this;
  }
  SparseMatrixIterator operator++(int n) { 
    SparseMatrixIterator tmp = *this;
    itr++; 
    return tmp;
  }

  IteratorType operator->() {
    return itr;
  }

  SparseMatrixIterator operator=(const SparseMatrixIterator& smi) {
    itr = smi.itr;
    return *this; 
  }

  bool operator==(const SparseMatrixIterator& smi) {
    return itr == smi.itr;
  }
  bool operator!=(const SparseMatrixIterator& smi) {
    return itr != smi.itr;
  }

  IteratorType itr;
};

struct SparseMatrixConstIterator
{
  typedef std::list<Cell>::const_iterator ConstIteratorType;

  SparseMatrixConstIterator() {}

  SparseMatrixConstIterator(const SparseMatrixConstIterator& smi) : itr(smi.itr) {}
  SparseMatrixConstIterator(const ConstIteratorType& smi_itr) : itr(smi_itr) {}

  SparseMatrixConstIterator operator++() { 
    ++itr; 
    return *this;
  }
  SparseMatrixConstIterator operator++(int n) { 
    SparseMatrixConstIterator tmp = *this;
    itr++; 
    return tmp;
  }

  ConstIteratorType operator->() {
    return itr;
  }

  SparseMatrixConstIterator operator=(const SparseMatrixConstIterator& smi) {
    itr = smi.itr;
    return *this; 
  }

  bool operator==(const SparseMatrixConstIterator& smi) {
    return itr == smi.itr;
  }
  bool operator!=(const SparseMatrixConstIterator& smi) {
    return itr != smi.itr;
  }

  ConstIteratorType itr;
};

class SparseMatrix 
{
public:
  typedef Cell::ValueType ValueType;

protected:
  SparseMatrix() : Rows_(0), Nrow_(0), Ncol_(0) {}
  virtual ~SparseMatrix() {}

public:
  void reset_matrix(uint n, uint m) {
    Nrow_ = n;
    Ncol_ = m;
    Rows_.clear();
    Rows_.resize(n, std::list<Cell>());
  }

  void set_row(const std::vector<ValueType>& r, uint ridx) {
    assert(r.size()==Ncol_ && ridx < Nrow_);
    Rows_[ridx].clear();
    for (uint i=0; i!=r.size(); ++i) 
      if (r[i] != DefVal_) 
	Rows_[ridx].push_back(Cell(i, r[i]));
  }

  void get_row(std::vector<ValueType>& r, uint ridx) const {
    assert(ridx < Nrow_);
    r.clear();
    r.resize(Ncol_, DefVal_);
    for (SparseMatrixConstIterator itr=begin(ridx); itr!=end(ridx); ++itr)
      r[itr->col] = itr->val;
  }

  bool is_null(uint ridx) const {
    return Rows_[ridx].size() == 0;
  }

  SparseMatrixConstIterator begin(uint i) const {
    return SparseMatrixConstIterator(Rows_[i].begin());
  }
  SparseMatrixConstIterator end(uint i) const {
    return SparseMatrixConstIterator(Rows_[i].end());
  }

  ValueType operator()(uint i, uint j) const { 
    assert(i < Nrow_ && j < Ncol_);
    for (SparseMatrixConstIterator itr=begin(i); itr!=end(i); ++itr) {
      if (itr->col == j) { return itr->val; }
      else if (itr->col > j) { return DefVal_; }
    }
    return DefVal_;
  }

  uint nrow() const {
    return Nrow_;
  }
  uint ncol() const {
    return Ncol_;
  }

protected:
  std::vector<std::list<Cell> > Rows_;
  uint Nrow_, Ncol_;

protected:
  static const ValueType DefVal_ = 0.0;
};

#endif

