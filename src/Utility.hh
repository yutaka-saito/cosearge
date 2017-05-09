#ifndef __INC_UTILITY_HH__
#define __INC_UTILITY_HH__

#include <string>
#include <sstream>
#include <boost/random.hpp>

#define INF 2e20

typedef boost::mt19937 RandomGenerator;

struct RandomShuffler {
  RandomShuffler(RandomGenerator& rng) : rng_(rng) {}
  uint operator()(uint i);

  RandomGenerator& rng_;
};


void
progress(const char* message);
void 
progress(const std::string& message);

uint
uabs(uint a, uint b);

// for scientific notations e.g. 1.500000e+1
double 
ss2a2f(std::stringstream& ss);
int
ss2a2i(std::stringstream& ss);

void
clear_ss(std::stringstream& ss);

bool
output_val(const std::string& file, double val);

#endif
