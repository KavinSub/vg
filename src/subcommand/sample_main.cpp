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

void multipath_node_ids(const MultipathAlignment &mp, set<int> &nodes){
	for(int i = 0; i < mp.subpath().size(); i++){
		for(int j = 0; j < mp.subpath()[i].path().mapping().size(); j++){
			int node = mp.subpath()[i].path().mapping()[j].position().node_id();
			nodes.insert(node);
		}
	}
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

	// *** Building the directed graph from a snarl ***
	vector<const Snarl*> top_snarls = snarl_manager -> children_of(nullptr);
	Snarl snarl = *top_snarls[0];

	list<Visit> path = sample_path(snarl, graph);
	for(Visit v: path){
		cout << v.node_id() << " ";
	}
	cout << endl;


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