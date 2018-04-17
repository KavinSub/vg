/**
* \file sample_main.cpp: sampling haplotypes from a sequence graph
*/

#include <stdlib.h>

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>
#include <string>
#include <set>
#include <random>

#include "subcommand.hpp"

#include "../haplotype_sampler.hpp"
#include "../snarls.hpp"
#include "../handle.hpp"
#include "../multipath_alignment.hpp"
#include "../vg.hpp"
#include "../phased_genome.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

SnarlManager* snarl_manager;
map<int, set<Snarl*> > id_map;
VG* graph;

// typedef struct node_t {
// 	int id;
// 	int forwardPaths;
// 	int reversePaths;
// 	int totalPaths;
// 	// vector<struct node_t> forwardNeighbors;
// 	// vector<struct node_t> reverseNeighbors;
// } RegularNode;

typedef struct count_graph_t {
	int head;
	int tail;
	map <int, int> count;
} CountGraph;

typedef struct directed_graph_t {
	int head;
	int tail;
	map<int, vector<int> > edges;
} DirectedGraph;

DirectedGraph reverse_graph(DirectedGraph &graph){
	DirectedGraph rgraph;
	rgraph.head = graph.tail;
	rgraph.tail = graph.head;

	for(map<int, vector<int> >::iterator iter = graph.edges.begin(); iter != graph.edges.end(); iter++){
		int node = iter->first;
		vector<int> neighbors = iter->second;
		for(int j = 0; j < neighbors.size(); j++){
			rgraph.edges[neighbors[j]].push_back(node);
		}
	}

	return rgraph;
}

bool is_head(DirectedGraph graph, int node){
	DirectedGraph rgraph = reverse_graph(graph);
	return rgraph.edges[node].size() == 0;
}

vector<int> topological_sort_graph(DirectedGraph graph){
	vector<int> L;
	vector<int> S;
	S.push_back(graph.head);

	while(S.size() > 0){
		int n = S.back();
		S.pop_back();
		L.push_back(n);
		vector<int> neighbors = graph.edges[n];
		graph.edges[n].clear();
		for(int i = 0; i < neighbors.size(); i++){
			int m = neighbors[i];
			if(is_head(graph, m)){
				S.push_back(m);
			}
		}
	}

	return L;
}

vector<int> get_nodes(DirectedGraph &graph){
	unordered_set<int> visited;
	vector<int> nodes;

	for(map<int, vector<int> >::iterator iter = graph.edges.begin(); iter != graph.edges.end(); iter++){
		int n = iter->first;
		vector<int> neighbors = iter->second;
		if(visited.find(n) == visited.end()){
			nodes.push_back(n);
			visited.insert(n);
		}
		for(int j = 0; j < neighbors.size(); j++){
			int m = neighbors[j];
			if(visited.find(m) == visited.end()){
				nodes.push_back(m);
				visited.insert(m);
			}
		}
	}

	return nodes;
}

void print_graph(DirectedGraph &graph){
	for(map<int, vector<int> >::iterator iter = graph.edges.begin(); iter != graph.edges.end(); iter++){
		int node = iter->first;
		vector<int> neighbors = iter->second;
		cout << node << ": ";
		for(int j = 0; j < neighbors.size(); j++){
			cout << neighbors[j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void print_count_graph(CountGraph &graph){
	for(map<int, int>::iterator iter = graph.count.begin(); iter != graph.count.end(); iter++){
		int node = iter->first;
		int count = iter->second;
		cout << node << ": " << count << endl;
	}
	cout << endl;
}

DirectedGraph build_graph(const Snarl &snarl, const HandleGraph* graph){

	DirectedGraph digraph;

	Visit start = snarl.start();
	Visit end = snarl.end();
	digraph.head = start.node_id();
	digraph.tail = end.node_id();

	handle_t begin = graph->get_handle(start.node_id());
	handle_t terminal = graph->get_handle(end.node_id());

	unordered_set<handle_t> visited;
	queue<handle_t> visit_queue;
	visit_queue.push(begin);

	while(visit_queue.size() > 0){
		handle_t current = visit_queue.front();
		visited.insert(current);
		visit_queue.pop();

		auto handle_other = [&](const handle_t& other){

			if(!visited.count(other) && graph->get_id(other) != graph->get_id(terminal)){
				visit_queue.push(other);
				visited.insert(other);
			}

			int current_id = graph->get_id(current);
			int other_id = graph->get_id(other);
			digraph.edges[current_id].push_back(other_id);
		};

		graph->follow_edges(current, false, handle_other);
	}

	return digraph;
}

CountGraph count_forward_paths(DirectedGraph &graph){
	CountGraph cgraph;
	cgraph.head = graph.head;
	cgraph.tail = graph.tail;

	cgraph.count[cgraph.head] = 1;
	vector<int> nodes = topological_sort_graph(graph);

	for(int i = 0; i < nodes.size(); i++){
		int n = nodes[i];
		vector<int> neighbors = graph.edges[n];
		for(int j = 0; j < neighbors.size(); j++){
			int m = neighbors[j];
			cgraph.count[m] += 1;
		}
	}

	return cgraph;
}

CountGraph count_reverse_paths(DirectedGraph &graph){
	DirectedGraph rgraph = reverse_graph(graph);
	CountGraph cgraph = count_forward_paths(rgraph);
	return cgraph;
}

CountGraph count_paths(DirectedGraph &graph){
	CountGraph cgraph;
	cgraph.head = graph.head;
	cgraph.tail = graph.tail;
	CountGraph fgraph = count_forward_paths(graph);
	CountGraph rgraph = count_reverse_paths(graph);

	vector<int> nodes = get_nodes(graph);

	for(int i = 0; i < nodes.size(); i++){
		int n = nodes[i];
		cgraph.count[n] = fgraph.count[n]*rgraph.count[n];
	}

	return cgraph;
}

void multipath_node_ids(const MultipathAlignment &mp, set<int> &nodes){
	for(int i = 0; i < mp.subpath().size(); i++){
		for(int j = 0; j < mp.subpath()[i].path().mapping().size(); j++){
			int node = mp.subpath()[i].path().mapping()[j].position().node_id();
			nodes.insert(node);
		}
	}
}

void print_vector(const vector<int> &vec){
	for(vector<int>::const_iterator it = vec.begin(); it != vec.end(); it++){
		cout << *it << endl;
	}
}

// Returns the probability distribution of paths conditioned on a node
static vector<double> get_distribution(CountGraph& cgraph, DirectedGraph& digraph, int node){
	double total = 0.0;
	for(int neighbor: digraph.edges[node])
		total += cgraph.count[neighbor];

	vector<double> dist;
	for(int neighbor: digraph.edges[node])
		dist.push_back(cgraph.count[neighbor]/total);

	return dist;
}

static list<Visit> sample_path(Snarl& snarl, HandleGraph* graph){
	DirectedGraph digraph = build_graph(snarl, graph);
	CountGraph cgraph = count_paths(digraph);

	int current = digraph.head;
	list<Visit> path;
	path.push_back(graph->to_visit(graph->get_handle(digraph.head)));

	while(current != digraph.tail){
		vector<double> dist = get_distribution(cgraph, digraph, current);
		default_random_engine generator;
		uniform_real_distribution<double> distribution(0.0, 1.0);

		double p = distribution(generator);
		double rp = 0.0;

		for(int i = 0; i < digraph.edges[current].size(); i++){
			int neighbor = digraph.edges[current][i];
			rp += dist[i];

			if(p <= rp){
				path.push_back(graph->to_visit(graph->get_handle(neighbor)));
				current = neighbor;
				break;
			}
		}
	}

	return path;
}

void help_sample(char** argv){
	cerr << "usage: " << argv[0] << " sample [options] .snarls .gamp .vg" << endl
	<< "Sample haplotypes from a sequence graph." << endl
	<< "options:" << endl
	<< "N/A" << endl;
	// << "	-s, --snarls FILE (required)" << endl;
}

const void snarl_inner(const Snarl* s){
	pair<unordered_set<Node*>, unordered_set<Edge*> > p = snarl_manager -> deep_contents(s, *graph, true);
	unordered_set<Node*> nodes = p.first;
	for(unordered_set<Node*>::const_iterator it = nodes.begin(); it != nodes.end(); it++){
		int node_id = (*it) -> id();
		id_map[node_id].insert((Snarl*) s);
	}
}

list<Visit> sequence;
bool follow_edges_callback(const handle_t& handle){
	// cout << "this is fun" << endl;
	Visit visit = graph -> to_visit(handle);
	sequence.push_back(visit);
	return true;
}

int main_sample(int argc, char** argv){

	if(argc < 4){
		help_sample(argv);
		return 1;
	}

	string snarls_name = argv[2];
	string gamp_name = argv[3];
	string graph_name = argv[4];

	cout << snarls_name << endl;
	cout << gamp_name << endl;
	cout << graph_name << endl;

	// int c;
	// optind = 2;
	// while(true){
	// 	static struct option long_options[] = {
	// 		{"snarls", required_argument, 0, 's'},
	// 		{0, 0, 0, 0}
	// 	};

	// 	int option_index = 0;
	// 	c = getopt_long (argc, argv, "s",
	//                      long_options, &option_index);

	// 	/* Detect the end of the options. */
	//     if (c == -1)
	//         break;

	//     switch(c){
	//     	case 's': 
	//     		// snarls_name = optarg;
	//       //       if (snarls_name.empty()) {
	//       //           cerr << "error:[vg sample] Must provide snarl file with -s." << endl;
	//       //           exit(1);
	//       //       }
	//             break;

	//         default:
	//             abort ();
	//     }
	// }

	// ********************************************************************************
    // 									Data Structures
    // ********************************************************************************

	// Snarl Manager
	snarl_manager = nullptr;
    if (!snarls_name.empty()) {
        ifstream snarl_stream(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg sample] Cannot open .snarls file " << snarls_name << endl;
            exit(1);
        }
        snarl_manager = new SnarlManager(snarl_stream);
    }
    cout << "Constructed Snarl Manager." << endl;

 //    // Multipath Alignments
 //    vector<MultipathAlignment> mp;
 //    function<void(MultipathAlignment&)> lambda = [&mp](MultipathAlignment& mp_aln) {
 //        mp.push_back(mp_aln);
 //    };
 //    get_input_file(gamp_name, [&](istream& in) {
 //        stream::for_each(in, lambda);
 //    });
 //    cout << "Constructed Multipath Alignments." << endl;
 //    cout << "Number of multipath alignments: " << mp.size() << endl;

    // Graph 
    get_input_file(graph_name, [&](istream& in){
    	graph = new VG(in);
    });
    cout << "Constructed Graph." << endl;

 //    // Phased Genome
 //    PhasedGenome p(*snarl_manager);
 //    cout << "Constructed Phased Genome." << endl;

    // ********************************************************************************
    // 									Graph testing
    // ********************************************************************************

	// *** Testing countgraph ***
    // DirectedGraph graph;
    // graph.head = 1;
    // graph.tail = 4;
    // graph.edges[1].push_back(2);
    // graph.edges[1].push_back(3);
    // graph.edges[2].push_back(4);
    // graph.edges[3].push_back(4);

    // CountGraph cgraph = count_paths(graph);
    // print_count_graph(cgraph);

	// *** Testing forward, reverse path counting ***
    // DirectedGraph rgraph = reverse_graph(graph);
    // print_graph(graph);
    // print_graph(rgraph);

    // CountGraph cgraph = count_forward_paths(graph);
    // print_count_graph(cgraph);
    // CountGraph crgraph = count_reverse_paths(graph);
    // print_count_graph(crgraph);

	// *** Building the directed graph from a snarl ***
	vector<const Snarl*> top_snarls = snarl_manager -> children_of(nullptr);
	Snarl snarl = *top_snarls[11];

	DirectedGraph digraph = build_graph(snarl, graph);
	CountGraph cgraph = count_paths(digraph);

	list<Visit> path = sample_path(snarl, graph);
	for(Visit v: path){
		cout << v.node_id() << " ";
	}
	cout << endl;

    // ********************************************************************************
    // 									Path counting
    // ********************************************************************************
    // vector<const Snarl*> snarls = snarl_manager -> children_of(nullptr);
    // build_graph(*snarls[0], *graph);


    // ********************************************************************************
    // 									Snarl to Multipath map
    // ********************************************************************************
    // // map: MultipathAlignments --> nodes
    // map<MultipathAlignment*, set<int> > mp_map;
    // for(int i = 0; i < mp.size(); i++){
    // 	set<int> nodes;
    // 	multipath_node_ids(mp[i], nodes);
    // 	mp_map[&mp[i]] = nodes;
    // }

    // // map: nodes --> snarls
    // snarl_manager -> for_each_snarl_preorder(snarl_inner);
    // cout << "Constructed nodes to snarls map." << endl;

    // // map: snarls --> MultipathAlignments
    // map<Snarl*, set<MultipathAlignment*> > snarl_map;
    // for(int i = 0; i < mp.size(); i++){
    // 	set<int> nodes = mp_map[&mp[i]];
    // 	for(set<int>::const_iterator node = nodes.begin(); node != nodes.end(); node++){
    // 		for(set<Snarl*>::const_iterator snarl = id_map[*node].begin(); snarl != id_map[*node].end(); snarl++){
    // 			snarl_map[*snarl].insert(&mp[i]);
    // 		} 
    // 	}
    // }
    // cout << "Constructed snarls to multipath map." << endl;

    // vector<const Snarl*> top_snarls = snarl_manager -> children_of(nullptr);
    // typedef vector<const Snarl*>::const_iterator iter;

    // // Snarl snarl = *top_snarls[0];
    // // Visit start = snarl.start();
    // // Visit end = snarl.end();
    // // cout << "Snarl ends: " << start.node_id() << " " << end.node_id() << endl;
    // // handle_t handle = graph -> get_handle(start.node_id(), false);
    // // graph -> follow_edges(handle, false, follow_edges_callback);
    // // cout << "Size of sequence: " << sequence.size() << endl;

    // // Initial Sequence
    // map<int, list<Visit>::iterator> seqmap;
    // for(iter it = top_snarls.begin(); it != top_snarls.end(); it++){
    // 	Snarl snarl = **it;
    // 	Visit start = snarl.start();

    // 	sequence.push_back(start);
    // 	handle_t handle = graph -> get_handle(start.node_id(), false);
    // 	graph -> follow_edges(handle, false, follow_edges_callback);

    // }
    // Snarl last_snarl = *top_snarls[top_snarls.size() - 1];
    // Visit end = last_snarl.end();
    // sequence.push_back(end);

    // // To verify sequence contains no overlaps
    // // unordered_set<Visit*> unique_visits;
    // // for(list<Visit>::iterator it = sequence.begin(); it != sequence.end(); it++){
    // // 	Visit* visit = &(*it);
    // // 	unique_visits.insert(visit);
    // // }
    // // cout << "Size of sequence: " << sequence.size() << endl;
    // // cout << "Size of unique visits: " << unique_visits.size() << endl;


	return 0;
}

static Subcommand vg_add("sample", "Sample haplotypes from a sequence graph", main_sample);