#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>

#include "Framework3DColocDiff.hh"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

int 
main(int argc, char** argv)
{
  Options opts;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
  opts.add_options(desc);

  po::variables_map vm;
  po::parsed_options parsed = 
    po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
  vector<string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  vector<po::option>::iterator new_end =
    remove_if(parsed.options.begin(), parsed.options.end(), 
	      bind(&po::option::unregistered, _1) );
  parsed.options.erase(new_end, parsed.options.end());
  po::store(parsed, vm);
  po::notify(vm);
  
  if (vm.count("help") || extra_args.size()<2) {
    cout << "3DColocDiff - perform various operations for two heatmaps" << endl
	 << endl
	 << "Usage:" << endl
	 << "  " << argv[0] << " [options] heatdir_smp heatdir_ctl" << endl
	 << endl
	 << desc << endl;
    return 1;
  }

  opts.parse_extra_args(extra_args);

  bool res = false;
  try {
    App app(opts);
    res = app.execute();
  }
  catch (const char* str) {
    cout << str << endl;
  }

  return !res;
}
