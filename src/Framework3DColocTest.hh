#ifndef __INC_FRAMEWORK_3D_COLOC_TEST_HH__
#define __INC_FRAMEWORK_3D_COLOC_TEST_HH__

#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "Utility.hh"

struct Options
{
  // input
  std::string in_dir;
  std::string in_elm;
  std::string in_gap;
  // output
  std::string out_stat;
  std::string out_dist;
  // others
  std::string suffix;
  uint nsmp;
  uint ncnk;
  uint seed;
  bool verbose;

  Options() : in_dir(), in_elm(), in_gap(), out_stat(), out_dist() {}

  void add_options(boost::program_options::options_description& opts);
  void parse_extra_args(const std::vector<std::string>& extra_args);
};

class App
{
public:
  App(const Options& opts) : opts_(opts) {}

  bool execute();

private:
  bool test();

private:
  const Options& opts_;
};

#endif
