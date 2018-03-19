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
#include <iterator>
#include <string>
#include <set>

#include "subcommand.hpp"

#include "../haplotype_sampler.hpp"
#include "../snarls.hpp"
#include "../multipath_alignment.hpp"
#include "../vg.hpp"
#include "../phased_genome.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

SnarlManager* snarl_manager;
map<int, vector<Snarl*> > id_map;
VG* graph;

void multipath_node_ids(const MultipathAlignment &mp, vector<int> &node_ids){
	set<int> encountered_ids;
	for(int i = 0; i < mp.subpath().size(); i++){
		for(int j = 0; j < mp.subpath()[i].path().mapping().size(); j++){
			int node_id = mp.subpath()[i].path().mapping()[j].position().node_id();
			if(encountered_ids.find(node_id) == encountered_ids.end()){
				node_ids.push_back(node_id);
				encountered_ids.insert(node_id);
			}
		}
	}
}

void print_vector(const vector<int> &vec){
	for(vector<int>::const_iterator it = vec.begin(); it != vec.end(); it++){
		cout << *it << endl;
	}
}

void help_sample(char** argv){
	cerr << "usage: " << argv[0] << " sample [options] .snarls .gamp" << endl
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
		if(id_map.count(node_id) == 0){
			vector<Snarl*> snarl_set;
			snarl_set.push_back((Snarl*) s);
			id_map[node_id] = snarl_set;
		}else{
			id_map[node_id].push_back((Snarl*) s);
		}
	}
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

	snarl_manager = nullptr;
    if (!snarls_name.empty()) {
        ifstream snarl_stream(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg sample] Cannot open .snarls file " << snarls_name << endl;
            exit(1);
        }
        snarl_manager = new SnarlManager(snarl_stream);
    }

    cout << "Snarl Manager created." << endl;

    vector<MultipathAlignment> mp;
    function<void(MultipathAlignment&)> lambda = [&mp](MultipathAlignment& mp_aln) {
        mp.push_back(mp_aln);
    };
    get_input_file(gamp_name, [&](istream& in) {
        stream::for_each(in, lambda);
    });

    cout << "Multipath alignments read." << endl;
    cout << "Number of multipath alignments: " << mp.size() << endl;

    get_input_file(graph_name, [&](istream& in){
    	graph = new VG(in);
    });

    cout << "Graph read." << endl;

    PhasedGenome p(*snarl_manager);
    cout << "Phased genome constructed." << endl;

    // Map from Multpath alignments to contained node ids
    map<MultipathAlignment*, vector<int> > mp_map;
    for(int i = 0; i < mp.size(); i++){
    	vector<int> node_ids;
    	multipath_node_ids(mp[i], node_ids);
    	mp_map[&mp[i]] = node_ids;
    }

    // Build map from node ids to snarls
    snarl_manager -> for_each_snarl_preorder(snarl_inner);

    cout << "Constructed id to snarl map." << endl;

    // Construct map from snarls to multipath alignment
    map<Snarl*, vector<MultipathAlignment*> > snarl_map;
    for(int i = 0; i < mp.size(); i++){
    	vector<int> node_ids = mp_map[&mp[i]];
    	unordered_set<Snarl*> encountered_snarls;
    	for(vector<int>::const_iterator it = node_ids.begin(); it != node_ids.end(); it++){
    		vector<Snarl*> snarls = id_map[*it];
    		for(vector<Snarl*>::const_iterator st = snarls.begin(); st != snarls.end(); st++){
    			if(encountered_snarls.find(*st) == encountered_snarls.end()){
	    			if(snarl_map.find(*st) == snarl_map.end()){
	    				vector<MultipathAlignment*> vmp;
	    				vmp.push_back(&mp[i]);
	    				snarl_map[*st] = vmp;
	    			}else{
	    				snarl_map[*st].push_back(&mp[i]);
	    			}
	    			encountered_snarls.insert(*st);
    			}
    		}
    	}
    }

    cout << "Constructed snarls to multipath map." << endl;

    vector<const Snarl*> snarls = snarl_manager -> children_of(nullptr);
    int rindex = rand() % snarls.size();
    Snarl *s = (Snarl*) snarls[rindex];

    cout << snarl_map[s].size() << endl;

	return 0;
}

static Subcommand vg_add("sample", "Sample haplotypes from a sequence graph", main_sample);