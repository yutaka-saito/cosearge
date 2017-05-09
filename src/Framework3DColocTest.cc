#include <cassert>

#include "Framework3DColocTest.hh"
#include "HypothesisTest.hh"
#include "Heatmap.hh"
#include "Bed.hh"

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
     po::value<uint>(&nsmp)->default_value(10000),
     "set the number of Monte Carlo samples")
    ("chunk",
     po::value<uint>(&ncnk)->default_value(6),
     "set the number of chromosome chunks")
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
  assert(extra_args.size() >= 5);
  in_dir = extra_args[0];
  in_elm = extra_args[1];
  in_gap = extra_args[2];
  out_stat = extra_args[3];
  out_dist = extra_args[4];
}

bool App::
execute()
{
  return test();
}

bool App::
test()
{
  bool res = false;

  progress("3DColocTest - test colocalization hypothesis on a heatmap");
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

  progress("calculating a test statistics");
  res = htest.test_stat();
  if (!res) return false;
    
  progress("calculating a null distribution");
  res = htest.null_dist(opts_.nsmp, opts_.verbose);
  if (!res) return false;
    
  progress("writing a test statistics");
  res = htest.output_stat(opts_.out_stat);
  if (!res) return false;
    
  progress("writing a null distribution");
  res = htest.output_dist(opts_.out_dist);
  if (!res) return false;

  progress("3DColocTest - finished");

  return true;
}
