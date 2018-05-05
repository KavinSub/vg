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


void help_sample(char** argv){
	cerr << "usage: " << argv[0] << " sample [options] .snarls .gamp .vg" << endl
	<< "Sample haplotypes from a sequence graph." << endl
	<< "options:" << endl
	<< "N/A" << endl;
	// << "	-s, --snarls FILE (required)" << endl;
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
	SnarlManager* snarl_manager;
    if (!snarls_name.empty()) {
        ifstream snarl_stream(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg sample] Cannot open .snarls file " << snarls_name << endl;
            exit(1);
        }
        snarl_manager = new SnarlManager(snarl_stream);
    }
    cout << "Constructed Snarl Manager." << endl;

    // Multipath Alignments
    vector<MultipathAlignment> mp;
    function<void(MultipathAlignment&)> lambda = [&mp](MultipathAlignment& mp_aln) {
        mp.push_back(mp_aln);
    };
    get_input_file(gamp_name, [&](istream& in) {
        stream::for_each(in, lambda);
    });
    cout << "Constructed Multipath Alignments." << endl;
    cout << "Number of multipath alignments: " << mp.size() << endl;

    // Graph 
    VG* graph;
    get_input_file(graph_name, [&](istream& in){
    	graph = new VG(in);
    });
    cout << "Constructed Graph." << endl;

    // Phased Genome
    PhasedGenome p(*snarl_manager);
    cout << "Constructed Phased Genome." << endl;

	return 0;
}

static Subcommand vg_add("sample", "Sample haplotypes from a sequence graph", main_sample);