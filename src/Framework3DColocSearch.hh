#ifndef __INC_FRAMEWORK_3D_COLOC_SEARCH_HH__
#define __INC_FRAMEWORK_3D_COLOC_SEARCH_HH__

#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "Utility.hh"

struct Options
{
  // input
  std::string indir_smp;
  std::string indir_ctl;
  std::string in_elm;
  std::string in_gap;
  // output
  std::string out_subset;
  std::string out_stat;
  // others
  std::string sfx_smp;
  std::string sfx_ctl;
  uint nset;
  uint nsmp;
  uint ncnk;
  uint nthr;
  uint seed;
  bool verbose;

  Options() : indir_smp(), indir_ctl(), in_elm(), in_gap(), 
	      out_subset(), out_stat(), sfx_smp(), sfx_ctl() {}

  void add_options(boost::program_options::options_description& opts);
  void parse_extra_args(const std::vector<std::string>& extra_args);
};

class App
{
public:
  App(const Options& opts) : opts_(opts) {}

  bool execute();

private:
  bool search();
  bool diff_search();
  
private:
  const Options& opts_;
};

#endif
