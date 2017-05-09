#include <cassert>

#include "Framework3DColocNormalize.hh"
#include "Heatmap.hh"

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
    ("proc",
     po::value<string>(&proc)->default_value("n"),
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
  assert(extra_args.size() >= 1);
  in_dir = extra_args[0];
}

bool App::
execute()
{
  return normalize();
}

bool App::
normalize()
{
  bool res = false;
  
  progress("3DColocNormalize - normalize a heatmap produced by hiclib");
  Heatmap hm;
  
  progress("loading from a heatmap directory");
  res = hm.load(opts_.in_dir, opts_.suffix);
  if (!res) return false;

  progress("performing the procedure");
  bool normalized = false;

  for (uint i=0; i!=opts_.proc.size(); ++i) {
    if (opts_.proc[i] == 'n') {
      cout << i+1 << ". normalization" << endl;
      hm.normalize(opts_.aggr, opts_.chrwise, opts_.verbose);
      normalized = true;
    }
    else if (opts_.proc[i] == 'l') {
      cout << i+1 << ". logarithm" << endl;
      hm.logarithm();
    }
    else if (opts_.proc[i] == 'm') {
      cout << i+1 << ". scaling" << endl;
      hm.scale();
    }
    else if (opts_.proc[i] == 't') {
      cout << i+1 << ". Tjong's scaling" << endl;
      hm.scaletj();
    }
    else {
      cout << i+1 << ". do nothing" << endl;
    }
    cout << flush;
  }

  if (normalized) {
    progress("writing heatmap statistics");
    res = hm.output_stat(opts_.in_dir, opts_.suffix + opts_.proc);
    if (!res) return false;
  }
  
  progress("writing heatmap values");
  if (opts_.chrsep) res = hm.output_val_chrsep(opts_.in_dir, opts_.suffix + opts_.proc);
  else res = hm.output_val(opts_.in_dir, opts_.suffix + opts_.proc);
  if (!res) return false;
  
  progress("3DColocNormalize - finished");
  
  return true; 
}

