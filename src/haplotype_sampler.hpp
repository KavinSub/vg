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

#endif