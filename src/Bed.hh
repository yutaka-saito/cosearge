#ifndef __INC_BED_HH__
#define __INC_BED_HH__

#include <cassert>
#include <string>
#include <vector>

#include "Utility.hh"

class Bed
{
public:
  Bed() : Chr_(), Start_(), Stop_(), 
	  Name_(), Value_(), Strand_(), Size_(0) {}
  ~Bed() {}

  bool load(const std::string& file);
  bool load_gap(const std::string& gap_file);
  void push_back(std::string ch, uint st, uint sp, 
		 std::string nm, double v, std::string d);
  void push_back(const Bed& elm, uint i);
  void pop_back();
  bool output(const std::string& file);
  void print();

  uint size() const {
    return Size_;
  }
  std::string chr(uint i) const {
    assert(i < Size_);
    return Chr_[i];
  }
  uint start(uint i) const {
    assert(i < Size_);
    return Start_[i];
  }
  uint stop(uint i) const {
    assert(i < Size_);
    return Stop_[i];
  }
  std::string name(uint i) const {
    assert(i < Size_);
    return Name_[i];
  }
  double value(uint i) const {
    assert(i < Size_);
    return Value_[i];
  }
  std::string strand(uint i) const {
    assert(i < Size_);
    return Strand_[i];
  }
  std::string line(uint i) const;

private:
  std::vector<std::string> Chr_;
  std::vector<uint> Start_;
  std::vector<uint> Stop_;
  std::vector<std::string> Name_;
  std::vector<double> Value_;
  std::vector<std::string> Strand_;
  uint Size_;
};

#endif
