#ifndef haplotype_sampler_hpp
#define haplotype_sampler_hpp


#include "haplotype_sampler.hpp"
#include "handle.hpp"
#include "snarls.hpp"
#include "multipath_alignment.hpp"
#include "vg.hpp"

#include <list>

using namespace vg;
using namespace std;

list<Visit> sample_path(Snarl&, HandleGraph*);
list<Visit> random_sequence(SnarlManager*, HandleGraph*);
map<int, list<Visit>::iterator> visit_map(list<Visit>);
void replace_subsequence(list<Visit>&, list<Visit>);
Snarl random_snarl(SnarlManager*);
map<Snarl*, set<MultipathAlignment*>> construct_snarl_map(SnarlManager*, vector<MultipathAlignment>&, VG*);

#endif