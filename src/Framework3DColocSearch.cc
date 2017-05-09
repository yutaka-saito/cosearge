#include <cassert>
#include <fstream>
#include <sstream>

#include "Framework3DColocSearch.hh"
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
    ("suffix_smp",
     po::value<string>(&sfx_smp)->default_value(""),
     "set the file name suffix")
    ("heatdir_ctl",
     po::value<string>(&indir_ctl)->default_value(""),
     "use a heatmap in this directory as control")
    ("suffix_ctl",
     po::value<string>(&sfx_ctl)->default_value(""),
     "set the file name suffix for control")
    ("subset",
     po::value<uint>(&nset)->default_value(10),
     "set the number of subsets")
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
  assert(extra_args.size() >= 5);
  indir_smp = extra_args[0];
  in_elm = extra_args[1];
  in_gap = extra_args[2];
  out_subset = extra_args[3];
  out_stat = extra_args[4];
}

bool App::
execute()
{
  if (opts_.indir_ctl == "") return search(); 
  else return diff_search();
}

bool App::
search()
{
  bool res = false;

  progress("3DColocSearch - search a heatmap for colocalized subsets of elements");
  Heatmap hm;
  Bed elm, gap;

  progress("loading from a heatmap directory");
  res = hm.load(opts_.indir_smp, opts_.sfx_smp);
  if (!res) return false;
  
  progress("loading from an element file");
  res = elm.load(opts_.in_elm);
  if (!res) return false;
  
  progress("loading from a gap file");
  res = gap.load_gap(opts_.in_gap);
  if (!res) return false;

  progress("preparing a hypothesis test");
  HypothesisTest htest(opts_.ncnk, opts_.seed, opts_.verbose, hm, elm, gap);

  for (uint i=0; i!=opts_.nset; ++i) {
    progress("searching a new subset");
    Bed subset;
    double zmax = - INF;
    //htest.search_subset_single(subset, zmax, opts_.nsmp, opts_.verbose);
    htest.search_subset(subset, zmax, opts_.nsmp, opts_.verbose, opts_.nthr);
    if (subset.size() == 0) break;

    progress("writing the optimal subset");
    string out_subset_name;
    string out_stat_name;
    stringstream out_subset_ss;
    stringstream out_stat_ss;
    out_subset_ss << opts_.out_subset << "." << i;
    out_stat_ss << opts_.out_stat << "." << i;
    out_subset_ss >> out_subset_name;
    out_stat_ss >> out_stat_name;
    res = subset.output(out_subset_name);
    if (!res) return false;
    res = output_val(out_stat_name, zmax);
    if (!res) return false;
  }  

  progress("3DColocSearch - finished");
    
  return true;
}

bool App::
diff_search()
{
  bool res = false;

  progress("3DColocSearch - search a heatmap for colocalized subsets of elements");
  Heatmap hm_smp, hm_ctl;
  Bed elm, gap;

  progress("loading from a heatmap directory for sample");
  res = hm_smp.load(opts_.indir_smp, opts_.sfx_smp);
  if (!res) return false;

  progress("loading from a heatmap directory for control");
  res = hm_ctl.load(opts_.indir_ctl, opts_.sfx_ctl);
  if (!res) return false;

  progress("loading from an element file");
  res = elm.load(opts_.in_elm);
  if (!res) return false;
  
  progress("loading from a gap file");
  res = gap.load_gap(opts_.in_gap);
  if (!res) return false;

  progress("preparing a hypothesis test");
  HypothesisTest htest(opts_.ncnk, opts_.seed, opts_.verbose, hm_smp, elm, gap);

  for (uint i=0; i!=opts_.nset; ++i) {
    progress("searching a new subset");
    Bed subset;
    double zmax = - INF;
    htest.diff_search_subset(subset, zmax, opts_.nsmp, opts_.verbose, opts_.nthr, hm_ctl);
    if (subset.size() == 0) break;

    progress("writing the optimal subset");
    string out_subset_name;
    string out_stat_name;
    stringstream out_subset_ss;
    stringstream out_stat_ss;
    out_subset_ss << opts_.out_subset << "." << i;
    out_stat_ss << opts_.out_stat << "." << i;
    out_subset_ss >> out_subset_name;
    out_stat_ss >> out_stat_name;
    res = subset.output(out_subset_name);
    if (!res) return false;
    res = output_val(out_stat_name, zmax);
    if (!res) return false;
  }  

  progress("3DColocSearch - finished");
    
  return true;
}
