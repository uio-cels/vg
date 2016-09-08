#include <iostream>
#include <fstream>
#include <tuple>
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

  // Should return the sequence on a given chromosome from start to end
  void get_sequence(string chromosome, int start, int end, string &sequence) {
    sequence = "tgc";
  }


  void make_events_map(map<string, vector<int> > &events, const vector<AltLoci> &altLoci)
  {
    for (auto alt : altLoci) {
      events[alt.chrom].push_back(alt.start_pos);
      events[alt.chrom].push_back(alt.end_pos);
    }
  }
  
  Node create_node(string chrom, int start, int stop, int id)
  {
    Node n;
    string seq;
    get_sequence(chrom, start, stop, seq);
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

  void make_blocks(const vector<AltLoci> &altLoci,
		   const vector<Chromosome> &chromosomes,
		   vector<Node> &nodes,
		   vector<Edge> &edges)
  {
    map<int, Node> node_ids;
    map<string, vector<int> > events;
    make_events_map(events, altLoci);

    map<tuple<string,int>, Node> starts;
    map<tuple<string,int>, Node> ends;
    int id = 1;

    for (auto &iter : events) {
      set<int> u_events(iter.second.begin(), iter.second.end());
      int start = 0;
      string chrom  = iter.first;
      for (auto stop : u_events){
	Node n = create_node(chrom, start, stop, id++);
	nodes.push_back(n);
	starts[make_tuple(chrom, start)] = n;
	ends[make_tuple(chrom, stop)] = n;
	node_ids[n.id()] = n;
	start = stop;
      }

      int stop = start;
      for (auto chr : chromosomes) {
	if (chr.name != chrom) continue;
	stop = chr.len;
	break;
      }

      Node n = create_node(chrom, start, stop, id++);
      nodes.push_back(n);
      starts[make_tuple(chrom, start)] = n;
      ends[make_tuple(chrom, stop)] = n;
      node_ids[n.id()] = n;
    }

    for (auto alt : altLoci) {
      Node n = create_node(alt.name, 0, alt.len, id++);
      nodes.push_back(n);
      
      Node pre_node = ends[make_tuple(alt.chrom, alt.start_pos)];
      Node post_node = starts[make_tuple(alt.chrom, alt.end_pos)];
      if (post_node.id() == 0) cout << alt.chrom << ", " << alt.end_pos << endl;
	
      edges.push_back(create_edge(pre_node, n));
      edges.push_back(create_edge(n, post_node));
    }
    for (auto edge : edges) cout << edge.from() <<  ": " << edge.to() << endl;
  }
  
  vector<Node> get_main_nodes(const map<string, vector<int> > &events,
			      map<tuple<string,int>, Node> &starts;
			      map<tuple<string,int>, Node> &ends;
			      ) {
    vector<Node> nodes;
    for (auto &iter : events) {
      set<int> u_events(iter.second.begin(), iter.second.end());
      int start = 0;
      string chrom  = iter.first;
      for (auto stop : u_events){
	Node n = create_node(chrom, start, stop, id++);
	nodes.push_back(n);
	starts[make_tuple(chrom, start)] = n;
	ends[make_tuple(chrom, stop)] = n;
	node_ids[n.id()] = n;
	start = stop;
      }

      int stop = start;
      for (auto chr : chromosomes) {
	if (chr.name != chrom) continue;
	stop = chr.len;
	break;
      }

      Node n = create_node(chrom, start, stop, id++);
      nodes.push_back(n);
      starts[make_tuple(chrom, start)] = n;
      ends[make_tuple(chrom, stop)] = n;
      node_ids[n.id()] = n;
    }
    return nodes;
  }
    
    

  map<AltLoci, Node> get_alt_nodes(const vector<AltLoci> &altLoci) {
    map<AltLoci, Node> nodes;
    for (auto alt : altLoci) {
      nodes[alt] = create_node(alt.name, 0, alt.len, id++);
    }
    return nodes;
  }

  
  vector<Edge> get_edges(const map<AltLoci, Node> & altNodes,
			 const map<tuple<string, int>, Node> &starts,
			 const map<tuple<string, int>, Node> &ends
			 ) {
    vector<Edge> edges;
    for (auto iter : altNodes) {
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


  void parse_data_files(string filename1, string filename2) {
    vector<Node> nodes;
    vector<Edge> edges;
    vector<AltLoci> altLoci;
    vector<Chromosome> chromosomes;
    parse_alt_loci(altLoci, filename1);
    parse_chrom_sizes(altLoci, chromosomes, filename2);
    make_blocks(altLoci, chromosomes, nodes, edges);
    set<Node*> pnodes;
    set<Edge*> pedges;
    for (auto node : nodes)
      pnodes.insert(&node);
    for (auto edge : edges)
      pedges.insert(&edge);
    VG graph = VG(pnodes, pedges);
  }
}
