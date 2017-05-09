#include <cassert>

#include "Framework3DColocDiff.hh"
#include "Heatmap.hh"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

void Options::  
add_options(po::options_description& opts)
{
  opts.add_options()
    ("suffix_smp",
     po::value<string>(&sfx_smp)->default_value(""),
     "set the file name suffix for heatmap_dir1")
    ("suffix_ctl",
     po::value<string>(&sfx_ctl)->default_value(""),
     "set the file name suffix for heatmap_dir2")
    ("proc",
     po::value<string>(&proc)->default_value("s"),
     "specify the procedure")
    ("aggregate",
     po::value<float>(&aggr)->default_value(1e-2),
     "aggregate neibor bins if sample size is smaller than this fraction")
    ("chrwise",
     po::value<bool>(&chrwise)->zero_tokens()->default_value(false),
     "normalize with chromosome-wise statistics")
    ("chrsep",
     po::value<bool>(&chrsep)->zero_tokens()->default_value(false),
     "output heatmap values with chromosomes separated")
    ("verbose",
     po::value<bool>(&verbose)->zero_tokens()->default_value(false),
     "make verbose reports to stdout")
    ;
}

void Options::
parse_extra_args(const vector<string>& extra_args)
{
  assert(extra_args.size() >= 2);
  indir_smp = extra_args[0];
  indir_ctl = extra_args[1];
}

bool App::
execute()
{
  return diff();
}

bool App::
diff()
{
  bool res = false;
  
  progress("3DColocDiff - perform various operations for two heatmaps");
  Heatmap hm_smp, hm_ctl;
  
  progress("loading from a heatmap directory for sample");
  res = hm_smp.load(opts_.indir_smp, opts_.sfx_smp);
  if (!res) return false;
  progress("loading from a heatmap directory for control");
  res = hm_ctl.load(opts_.indir_ctl, opts_.sfx_ctl);
  if (!res) return false;

  progress("performing the procedure");
  bool oneside = false;
  bool normalized = false;

  for (uint i=0; i!=opts_.proc.size(); ++i) {
    if (opts_.proc[i] == 's') {
      cout << i+1 << ". subtraction, heatmap_smp - heatmap_ctl" << endl;
      hm_smp.subtract(hm_ctl);
      oneside = true; 
    }
    else if (opts_.proc[i] == 'd') {
      cout << i+1 << ". division, heatmap_smp / heatmap_ctl" << endl;
      hm_smp.divide(hm_ctl);
      oneside = true;
    }
    else if (opts_.proc[i] == 'n') {
      cout << i+1 << ". normalization, bar(heatmap_smp)" << endl;
      hm_smp.normalize(opts_.aggr, opts_.chrwise, opts_.verbose);
      normalized = true;
      if (oneside) continue;
      cout << i+1 << ". normalization, bar(heatmap_ctl)" << endl;
      hm_ctl.normalize(opts_.aggr, opts_.chrwise, opts_.verbose);
    }
    else if (opts_.proc[i] == 'l') {
      cout << i+1 << ". logarithm, log(heatmap_smp)" << endl;
      hm_smp.logarithm();
      if (oneside) continue;
      cout << i+1 << ". logarithm, log(heatmap_ctl)" << endl;
      hm_ctl.logarithm();
    }
    else if (opts_.proc[i] == 'm') {
      cout << i+1 << ". scaling, scale(heatmap_smp)" << endl;
      hm_ctl.scale();
      if (oneside) continue;
      cout << i+1 << ". scaling, scale(heatmap_ctl)" << endl;
      hm_smp.scale();
    }
    else {
      cout << i+1 << ". do nothing" << endl;
    }
    cout << flush;
  }

  if (normalized) {
    progress("writing heatmap_smp statistics");
    res = hm_smp.output_stat(opts_.indir_smp, opts_.sfx_smp + opts_.proc);
    if (!res) return false;
  }
  
  progress("writing heatmap_smp values");
  if (opts_.chrsep) res = hm_smp.output_val_chrsep(opts_.indir_smp, opts_.sfx_smp + opts_.proc);
  else res = hm_smp.output_val(opts_.indir_smp, opts_.sfx_smp + opts_.proc);
  if (!res) return false;
  
  progress("3DColocDiff - finished");
  
  return true; 
}

