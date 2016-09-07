#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace vg{
  using namespace std;

  class AltLoci {
  public:
    string name;
    string chrom;
    int len;
    int start_pos;
    int end_pos;
    AltLoci(string name, string chrom, int start, int end);
    
  };

  class Chromosome {
  public:
    string name;
    int len;
    Chromosome(string name, int len);
  };
  
  Chromosome::Chromosome(string name, int len): name(name), len(len){}
  AltLoci::AltLoci(string name, string chrom, int start, int end): name(name), chrom(chrom), start_pos(start), end_pos(end){}

  void parse_alt_loci(vector<AltLoci> &altLoci, string filename) {
    string line;
    vector<string> parts;
    ifstream file1(filename);
    while ( getline(file1, line) ) {
      if (line.empty()) continue;
      if (line[0] == '#') continue;
      boost::split(parts, line, boost::is_any_of("\t "));
      if (! boost::starts_with(parts[7], "ALT")) continue;
      altLoci.push_back(AltLoci(parts[5], "chr" + parts[1],  stoi(parts[2]), stoi(parts[3])));
    }
    file1.close();
  }

  void parse_chrom_sizes(vector<AltLoci> &altLoci, vector<Chromosome> &chromosomes, string filename) {
    ifstream file2(filename);
    vector<string> parts;
    string line;
    while (getline(file2, line))  {
	boost::split(parts, line, boost::is_any_of("\t "));
	if (boost::ends_with(parts[0], "alt")) {
	  vector<string> subnames;
	  boost::split(subnames, parts[0], boost::is_any_of("_"));
	  string name  = subnames[1];
	  size_t vpos = name.find("v");
	  name.replace(vpos, 1, ".");
	  cout << name << endl;
	  for (auto alt : altLoci) {
	    if (alt.name == name) {
	      alt.len = stoi(parts[1]);
	      cout << alt.len << endl;
	      break;
	    }
	  }
	} else if (boost::ends_with(parts[0], "random")) {
	  continue;
	} else if (boost::starts_with(parts[0], "chrUn")) {
	  continue;
	} else {
	  chromosomes.push_back(Chromosome(parts[0], stoi(parts[1])));
	}
      }
    file2.close();
  }

  void make_blocks(const vector<AltLoci> &altLoci, const vector<Chromosome> &chromosomes)
  {
    map<string, vector<int>> events;
    for (auto chromosome : chromosomes)
      events[chromosome.name] = vector<int>();
    for (auto alt : altLoci) {
      events[alt.chrom].push_back(alt.start_pos);
      events[alt.chrom].push_back(alt.end_pos);
    }
    vector<AltLoci> nodes;

    for (auto &iter : events) {
      set<int> u_events(iter.second.begin(), iter.second.end());
      int old_int = 0;
      cout << iter.first << endl; 
      for (auto event : u_events){
	cout << old_int << "," << event << endl;
	old_int = event;
      }
    }
  }



void parse_data_files(string filename1, string filename2) {
  
  vector<AltLoci> altLoci;
  vector<Chromosome> chromosomes;
  parse_alt_loci(altLoci, filename1);
  parse_chrom_sizes(altLoci, chromosomes, filename2);
  //for (auto chr : chromosomes) cout << chr.name << endl;
  make_blocks(altLoci, chromosomes);
}

}
