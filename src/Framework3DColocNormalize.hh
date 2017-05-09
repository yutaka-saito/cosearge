#ifndef __INC_FRAMEWORK_3D_COLOC_NORMALIZE_HH__
#define __INC_FRAMEWORK_3D_COLOC_NORMALIZE_HH__

#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "Utility.hh"

struct Options
{
  // input
  std::string in_dir;
  // others
  std::string suffix;
  std::string proc;
  float aggr;
  bool chrwise;
  bool chrsep;
  bool verbose;

  Options() : in_dir(), suffix(), proc() {}

  void add_options(boost::program_options::options_description& opts);
  void parse_extra_args(const std::vector<std::string>& extra_args);
};

class App
{
public:
  App(const Options& opts) : opts_(opts) {}

  bool execute();

private:
  bool normalize();

private:
  const Options& opts_;
};

#endif	
