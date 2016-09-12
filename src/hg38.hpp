#include <iostream>
#include <fstream>
#include <tuple>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include "vg.pb.h"
#include "vg.hpp"

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

    bool operator<(AltLoci other) const {
      if (chrom!=other.chrom) return chrom <other.chrom;
      return chrom<other.chrom;
    }

  };

  map<string, string> get_chr_map() {
    map<string, string> chr_map;
    ifstream mapfile("chrmap.txt");
    string line;
    vector<string> parts;
    while (getline(mapfile, line)) {
      boost::split(parts, line, boost::is_any_of("\t "));
      chr_map[parts[0]] = parts[1];
    }
    mapfile.close();
    return chr_map;
  }

  class AltReference {
  public:
    FastaReference main_ref, alt_ref;
    map<string, string> chr_map;
    AltReference(FastaReference &mainref, FastaReference &altref) {
      main_ref = mainref;
      alt_ref = altref;
      chr_map = get_chr_map();
    }
    void get_sequence(string chromosome, int start, int end, string &sequence) {
      string chr_id = chr_map[chromosome];
      cout << chromosome << endl;
      if (boost::ends_with(chromosome, "alt")){
	cout << "Find: "<< chr_id << ":" << start << "," <<end << endl;
	//for (auto seqname : alt_ref.index->sequenceNames) cout << seqname << endl;
	sequence = alt_ref.getSubSequence(chr_id, start, end);
      }
      else {
	cout << "Find: "<< chr_id << ":" << start << "," <<end << endl;
	//for (auto seqname : main_ref.index->sequenceNames) cout << seqname << endl;
	sequence = main_ref.getSubSequence(chr_id, start, end);
      }
    }
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

  // Read through chrom sizes file and
  // a) create chromosome objects for each chromosome
  // b) update altloci lenght for all alt loci
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
	  for (auto alt : altLoci) {
	    if (alt.name == name) {
	      alt.len = stoi(parts[1]);
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


  void make_events_map(map<string, vector<int> > &events, const vector<AltLoci> &altLoci, const vector<Chromosome> &chromosomes)
  {
    for (auto alt : altLoci) {
      events[alt.chrom].push_back(alt.start_pos);
      events[alt.chrom].push_back(alt.end_pos);
    }
    for (auto chr : chromosomes)
      events[chr.name].push_back(chr.len);
  }
  
  Node create_node(string chrom, int start, int stop, int id, AltReference &reference)
  {
    Node n;
    string seq;
    reference.get_sequence(chrom, start, stop, seq);
    n.set_sequence(seq);
    n.set_id(id);
    n.set_name(chrom+boost::lexical_cast<string>(start)); 
    return n;
  }

  Edge create_edge(Node from, Node to) {
    Edge e;
    e.set_from(from.id());
    e.set_to(to.id());
    return e;
  }
  
  vector<Node> get_main_nodes(const map<string, vector<int> > &events,
			      map<tuple<string,int>, Node> &starts,
			      map<tuple<string,int>, Node> &ends,
			      AltReference &reference
			      ) {
    vector<Node> nodes;
    int id = 1;
    for (auto &iter : events) {
      set<int> u_events(iter.second.begin(), iter.second.end());
      int start = 0;
      string chrom  = iter.first;
      for (auto stop : u_events){
	int len = stop-start;
	int steps = len/1000;
	int i=start;
	Node n = create_node(chrom, i, min(i+1000, stop), id++, reference);
	nodes.push_back(n);
	starts[make_tuple(chrom, start)] = n;
	i = min(i+1000, stop);
	for (; i<len+steps*1000; i+=1000) {
	  n = create_node(chrom, i, i+1000, id++, reference);
	  nodes.push_back(n);
	}
	if (i!=stop) {
	  n = create_node(chrom, i, stop, id++, reference);
	  nodes.push_back(n);
	}
	ends[make_tuple(chrom, stop)] = n;
	start = stop;
      }
    }
    return nodes;
  }

  map<AltLoci, Node> get_alt_nodes(const vector<AltLoci> &altLoci, int id, AltReference &reference) {
    map<AltLoci, Node> nodes;
    for (auto alt : altLoci) {
      nodes[alt] = create_node(alt.name, 0, alt.len, id++, reference);
    }
    return nodes;
  }

  
  vector<Edge> get_edges(const map<AltLoci, Node> & altNodes,
			 map<tuple<string, int>, Node> &starts,
			 map<tuple<string, int>, Node> &ends
			 ) {
    vector<Edge> edges;
    for (auto &iter : altNodes) {
      AltLoci alt = iter.first;
      Node n = iter.second;
      Node pre_node = ends[make_tuple(alt.chrom, alt.start_pos)];
      Node post_node = starts[make_tuple(alt.chrom, alt.end_pos)];
      if (post_node.id() == 0) cout << alt.chrom << ", " << alt.end_pos << endl;
	
      edges.push_back(create_edge(pre_node, n));
      edges.push_back(create_edge(n, post_node));
    }
    for (auto edge : edges) cout << edge.from() <<  ": " << edge.to() << endl;
    return edges;
  }

  void test_graph_gen() {
    Node a, b;
    a.set_sequence("cgt");
    a.set_name("A");
    a.set_id(0);

    b.set_sequence("ttg");
    b.set_name("B");
    b.set_id(1);

    Edge e;
    e.set_from(0);
    e.set_to(1);
    set<Node*> nodes;
    nodes.insert(&a);
    nodes.insert(&b);
    set<Edge*> edges;
    edges.insert(&e);
    VG graph(nodes, edges);
  }

  void parse_data_files(string filename1,
			string filename2,
			AltReference &reference
			) {
    vector<Node> nodes;
    vector<Edge> edges;
    vector<AltLoci> altLoci, altLociTmp;
    vector<Chromosome> chromosomes, chromosomesTmp;
    parse_alt_loci(altLociTmp, filename1);
    parse_chrom_sizes(altLociTmp, chromosomesTmp, filename2);
    for (auto alt : altLociTmp) 
      if (alt.chrom=="chr1") altLoci.push_back(alt);
    for (auto chr : chromosomesTmp)
      if (chr.name == "chr1") chromosomes.push_back(chr);
    
    map<string, vector<int> > events;
    make_events_map(events,  altLoci, chromosomes);
    map<tuple<string,int>, Node> starts;
    map<tuple<string,int>, Node> ends;
    nodes = get_main_nodes(events, starts, ends, reference);
    auto alt_nodes = get_alt_nodes(altLoci, nodes.size()+1, reference);
    edges = get_edges(alt_nodes, starts, ends);
    cout << "Got edges" << endl;
    for (auto node : nodes)
      cout << node.id() << "," << node.name() << endl;
    for (auto edge : edges)
      cout << edge.from() << "-" << edge.to() << endl;
      
    set<Node*> pnodes;
    set<Edge*> pedges;
    cout << "Creating nodes" << endl;


    for (auto node : nodes) {
      Node* n = new Node(node);
      pnodes.insert(n);
    }
    for (auto &node : alt_nodes) {
      Node* n = new Node(node.second);
      pnodes.insert(n);
    }
    cout << "Creating edges" << endl;
    for (auto edge : edges) {
      Edge* e = new Edge(edge);
      pedges.insert(e);
    }
    cout << "Creating graph" << endl;
    VG graph = VG(pnodes, pedges);
    cout << "Made graph" << endl;
  }
}
