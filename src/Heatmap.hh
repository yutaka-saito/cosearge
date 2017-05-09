#ifndef __INC_HEATMAP_HH__
#define __INC_HEATMAP_HH__

#include <cassert>
#include <string>
#include <vector>

#include "Utility.hh"
#include "SparseMatrix.hh"

class Heatmap : public SparseMatrix
{
public:
  Heatmap() : SparseMatrix() {}
  ~Heatmap() {}

  bool load(const std::string& in_dir, const std::string& suffix);
  bool output_stat(const std::string& out_dir, const std::string& suffix) const;
  bool output_val(const std::string& out_dir, const std::string& suffix) const;
  bool output_val(const std::string& filename, 
		  uint st_row, uint sp_row, uint st_col, uint sp_col) const;
  bool output_val_chrsep(const std::string& out_dir, const std::string& suffix);

  void normalize(float aggr, bool chrwise, bool verbose);
  void logarithm();
  void scale();
  void scaletj(); // scaling described in [Tjong, Genome Res, 2012]
  void subtract(const Heatmap& hm);
  void divide(const Heatmap& hm);
  bool compatible(const Heatmap& hm) const;

  void print();

  std::string chr_label(uint i) const {
    assert(i < Nchr_);
    return ChrLabel_[i];
  }
  uint chr_start(uint i) const {
    assert(i < Nchr_);
    return ChrStart_[i];
  }
  uint chr_stop(uint i) const {
    assert(i < Nchr_);
    if (i == Nchr_-1) return Nbin_;
    else return ChrStart_[i+1];
  }
  uint nbin() const {
    return Nbin_;
  }
  uint nchr() const {
    return Nchr_;
  }
  uint resol() const {
    return Resol_;
  }

  bool pos2bin(const std::string label, const uint pos, uint& chr, uint& bin) const;
  bool pos2bin(const std::string label, const uint pos, uint& bin) const;
  bool bin2pos(const uint bin, std::string& label, uint& pos) const;

private:
  // IntraMean_[chri][i] and IntraStdev_[chri][i] are statistics for 
  // intra-chromosomal interactions in chromosome chri
  // whose distance = i * Resol_.
  // InterMean_[chri][chrj] and InterStdev_[chri][chrj] are statistics for 
  // inter-chromosomal interactions between chromosomes chri and chrj.
  // All*_ are statistics for 
  // all chromosomal interactions without chromosomes distinguished.
  std::vector<std::vector<ValueType> > IntraMean_;
  std::vector<std::vector<ValueType> > IntraStdev_; 
  std::vector<std::vector<ValueType> > InterMean_;
  std::vector<std::vector<ValueType> > InterStdev_;
  std::vector<ValueType> AllIntraMean_;
  std::vector<ValueType> AllIntraStdev_; 
  ValueType AllInterMean_;
  ValueType AllInterStdev_;
  std::vector<uint> ChrIdx_; // chromosome number for each bin
  std::vector<uint> PosIdx_; // nucleotide coordinate for each bin
  std::vector<std::string> ChrLabel_; // label (e.g. chrX) for each chromosome
  std::vector<uint> ChrStart_; // start bin coordinate for each chromosome
  uint Nbin_;
  uint Nchr_;
  uint Resol_;
};

#endif
