#ifndef __INC_FRAMEWORK_3D_COLOC_DIFF_HH__
#define __INC_FRAMEWORK_3D_COLOC_DIFF_HH__

#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "Utility.hh"

struct Options
{
  // input
  std::string indir_smp;
  std::string indir_ctl;
  // others
  std::string sfx_smp;
  std::string sfx_ctl;
  std::string proc;
  float aggr;
  bool chrwise;
  bool chrsep;
  bool verbose;

  Options() : indir_smp(), indir_ctl(), sfx_smp(), sfx_ctl(), proc() {}

  void add_options(boost::program_options::options_description& opts);
  void parse_extra_args(const std::vector<std::string>& extra_args);
};

class App
{
public:
  App(const Options& opts) : opts_(opts) {}

  bool execute();

private:
  bool diff();

private:
  const Options& opts_;
};

#endif
