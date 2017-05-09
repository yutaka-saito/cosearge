#ifndef __INC_HYPOTHESIS_TEST_HH__
#define __INC_HYPOTHESIS_TEST_HH__

#include <string>
#include <vector>
#include <map>

#include "Utility.hh"
#include "Heatmap.hh"
#include "Bed.hh"

struct ChromosomeChunk
{
public:
  ChromosomeChunk(const std::vector<uint>& elm, uint start, uint stop);

  void shuffle(RandomGenerator& rng);

public:
  // Elm_ contains whole-genome bin coordinates of elements 
  // where Dist_[i] = Elm_[i+1] - Elm_[i]
  // shuffle()-ing Elm_ permutates the order of Dist_ members 
  // while preserving their values.
  // This chunk corresponds to [Start_, Stop_)-th bin region.
  std::vector<uint> Elm_; 
  std::vector<uint> Dist_;
  std::vector<uint> ShfElm_;
  std::vector<uint> ShfDist_;
  uint Nelm_;
  uint Dsum_;
  uint Start_;
  uint Stop_;
};

class HypothesisTest
{
public:
  HypothesisTest(uint ncnk, uint seed, bool verbose, 
		 const Heatmap& hm, const Bed& elm, const Bed& gap);
  ~HypothesisTest() {}

  bool test_stat(bool shuffled = false);
  bool null_stat();
  bool null_dist(uint nsmp, bool verbose);
  bool output_stat(const std::string& filename);
  bool output_dist(const std::string& filename);
  bool get_stat(double& stat, double& p, double& z);
  void print();

  void search_subset_single(Bed& subset, double& zmax, uint nsmp, bool verbose);
  void search_subset(Bed& subset, double& zmax, uint nsmp, bool verbose, uint nthr);
  void diff_search_subset(Bed& subset, double& zmax, uint nsmp, bool verbose, uint nthr, 
			  const Heatmap& ctl);

  bool look_interactive(uint nsmp, bool verbose, uint nthr);

private:
  void calc_pval();
  void calc_zval();

  bool parse_position(const std::string buf, 
		      std::string& label, uint& pos, uint& bin);
  bool format_position(std::string& buf, 
		       const std::string label, const uint pos, const uint bin);

private:
  // This class stores data in whole-genome bin coordinates, 
  // separately for each chromosome chunk. 
  // Heatmap class stores in whole-genome bin coordinates, while 
  // Bed class stores in chromosome-wise nucleotide coodinates.
  // LongArm_[i][j] and ShortArm_[i][j] are the j-th chunk from  
  // long and short arms of the i-th chromosome, respectively. 
  // Bin2Elm_ is a map from whole-genome bin coordinates to 
  // indexes of Elm_ contained there. 
  std::vector<std::vector<ChromosomeChunk> > ShortArm_;
  std::vector<std::vector<ChromosomeChunk> > LongArm_;
  std::map<uint,std::vector<uint> > Bin2Elm_;
  uint Nchr_;
  uint Ncnk_;
  std::vector<double> NullDist_;
  double TestStat_;
  double NullStat_;
  double Pval_;
  double Zval_;
  const Heatmap& Hm_;
  const Bed& Elm_;
  const Bed& Gap_;
  const uint seed_;
  RandomGenerator rng_;

private:
  static const double ZvalMax_ = INF;
};

#endif
