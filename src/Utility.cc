#include <cstdlib>
#include <iostream>
#include <fstream>

#include "Utility.hh"

using namespace std;
using namespace boost;

uint RandomShuffler::
operator()(uint i)
{
  if (i==1) return 0; // necessary since uniform_int(a,b) does not support a==b 

  uniform_int<> unif_dist(0,i-1);
  variate_generator<RandomGenerator&, uniform_int<> > unif_samp(rng_, unif_dist);

  return unif_samp();
}

void
progress(const char* message) {
  cout << message << endl << flush;
}

void
progress(const string& message) {
  cout << message << endl << flush;
}

uint
uabs(uint a, uint b) {
  if (a > b) return a - b;
  else return b - a;
}

double 
ss2a2f(stringstream& ss) {
  string buf;
  ss >> buf;
  return atof(buf.c_str());
}

int
ss2a2i(stringstream& ss) {
  string buf;
  ss >> buf;
  return (int) atof(buf.c_str()); 
}

void
clear_ss(stringstream& ss) {
  ss.str("");
  ss.clear(stringstream::goodbit);
}

bool
output_val(const string& file, double val) {
  ofstream ofs(file.c_str());
  if (!ofs) return false;
  ofs << val << endl;
  ofs.close();
  return true;
}
