#include <cassert>
#include <fstream>
#include <sstream>

#include "Framework3DColocLook.hh"
#include "Heatmap.hh"
#include "Bed.hh"
#include "HypothesisTest.hh"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

void Options::  
add_options(po::options_description& opts)
{
  opts.add_options()
    ("suffix",
     po::value<string>(&suffix)->default_value(""),
     "set the file name suffix")
    ("sample",
     po::value<uint>(&nsmp)->default_value(100),
     "set the number of Monte Carlo samples")
    ("chunk",
     po::value<uint>(&ncnk)->default_value(6),
     "set the number of chromosome chunks")
    ("thread",
     po::value<uint>(&nthr)->default_value(1),
     "set the number of threads")
    ("random", 
     po::value<uint>(&seed)->default_value(1985), 
     "set the seed for random generator")
    ("verbose",
     po::value<bool>(&verbose)->zero_tokens()->default_value(false),
     "make verbose reports to stdout")
    ;
}

void Options::
parse_extra_args(const vector<string>& extra_args)
{
  assert(extra_args.size() >= 3);
  in_dir = extra_args[0];
  in_elm = extra_args[1];
  in_gap = extra_args[2];
}

bool App::
execute()
{
  return look();
}

bool App::
look()
{
  bool res = false;

  progress("3DColocLook - Look at a heatmap and test statistics interactively");
  Heatmap hm;
  Bed elm, gap;

  progress("loading from a heatmap directory");
  res = hm.load(opts_.in_dir, opts_.suffix);
  if (!res) return false;
  
  progress("loading from an element file");
  res = elm.load(opts_.in_elm);
  if (!res) return false;
  
  progress("loading from a gap file");
  res = gap.load_gap(opts_.in_gap);
  if (!res) return false;

  progress("preparing a hypothesis test");
  HypothesisTest htest(opts_.ncnk, opts_.seed, opts_.verbose, hm, elm, gap);

  progress("looking at interactively");
  while (htest.look_interactive(opts_.nsmp, opts_.verbose, opts_.nthr));

  progress("3DColocLook - finished");
    
  return true;
}
