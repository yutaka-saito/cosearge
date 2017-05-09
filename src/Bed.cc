#include <iostream>
#include <fstream>
#include <sstream>

#include "Bed.hh"

using namespace std;

bool Bed::
load(const string& file) 
{
  ifstream ifs(file.c_str());
  string line;

  if (!ifs) return false;

  while (getline(ifs, line) != NULL) {
    string ch, nm, d;
    uint st, sp;
    double v;
    stringstream ss;
    ss << line;
    ss >> ch >> st >> sp >> nm >> v >> d;
    push_back(ch, st, sp, nm, v, d);
  }
  ifs.close();

  assert(Size_ > 1);

  //print();
  return true;
}

bool Bed::
load_gap(const string& gap_file) 
{
  ifstream ifs(gap_file.c_str());
  string line;

  if (!ifs) return false;

  while (getline(ifs, line) != NULL) {
    string ch, nm, buf;
    uint st, sp;
    stringstream ss;
    ss << line;
    ss >> buf >> ch >> st >> sp >> buf >> buf >> buf >> nm;
    if (nm == "centromere") 
      push_back(ch, st, sp, nm, 0, "+");
  }
  ifs.close();

  //print();
  return true;
}

void Bed::
push_back(string ch, uint st, uint sp, 
	  string nm, double v, string d)
{
  Chr_.push_back(ch);
  Start_.push_back(st);
  Stop_.push_back(sp);
  Name_.push_back(nm);
  Value_.push_back(v);
  Strand_.push_back(d);
  Size_++;
}

void Bed::
push_back(const Bed& elm, uint i)
{
  push_back(elm.chr(i), elm.start(i), elm.stop(i), 
	    elm.name(i), elm.value(i), elm.strand(i));
}

void Bed::
pop_back()
{
  assert(Size_ > 0);
  Chr_.pop_back();
  Start_.pop_back();
  Stop_.pop_back();
  Name_.pop_back();
  Strand_.pop_back();
  Value_.pop_back();
  Size_--;
}

bool Bed::
output(const string& file) 
{
  ofstream ofs(file.c_str());
  if (!ofs) return false;

  for (uint i=0; i!=Size_; ++i) 
    ofs << line(i);
  ofs.close();

  return true;
}

string Bed::
line(uint i) const
{
  stringstream ss;
  ss << Chr_[i] << "\t" << Start_[i] << "\t" << Stop_[i] << "\t" 
     << Name_[i] << "\t" << Value_[i] << "\t" << Strand_[i] << endl;

  return ss.str();
}

void Bed::
print()
{
  cout << "print Bed class:" << endl << endl << flush;

  for (uint i=0; i!=Size_; ++i) 
    cout << Chr_[i] << "\t" << Start_[i] << "\t" << Stop_[i] << "\t" 
	 << Name_[i] << "\t" << Value_[i] << "\t" << Strand_[i] << endl << flush;
}
