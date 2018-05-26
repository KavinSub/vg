#ifndef haplotype_sampler_hpp
#define haplotype_sampler_hpp


#include "haplotype_sampler.hpp"
#include "handle.hpp"
#include "snarls.hpp"
#include "multipath_alignment.hpp"
#include "vg.hpp"
#include "phased_genome.hpp"

#include <list>

using namespace vg;
using namespace std;

void initialize_sampler();
list<Visit> sample_path(const Snarl&, HandleGraph*);
list<Visit> random_sequence(SnarlManager*, HandleGraph*);
map<int, list<Visit>::iterator> visit_map(list<Visit>);
void replace_subsequence(list<Visit>&, list<Visit>);
const Snarl* random_snarl(SnarlManager*);
const Snarl* random_snarl(vector<const Snarl*>&);
map<const Snarl*, set<MultipathAlignment*>> construct_snarl_map(SnarlManager*, vector<MultipathAlignment>&, VG*);
void naive_sampling(int, string, SnarlManager*, vector<MultipathAlignment>&, VG*, PhasedGenome&);

#endif