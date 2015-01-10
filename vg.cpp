#include "vg.hpp"

namespace vg {

using namespace std;

// construct from a stream of protobufs
VG::VG(istream& in) {
    init();

    uint64_t count = 0;

    ::google::protobuf::io::ZeroCopyInputStream *raw_in =
          new ::google::protobuf::io::IstreamInputStream(&in);
    ::google::protobuf::io::CodedInputStream *coded_in =
          new ::google::protobuf::io::CodedInputStream(raw_in);

    coded_in->ReadVarint64(&count);
    delete coded_in;

    if (show_progress) {
        progress_message = "loading graph";
        progress = new ProgressBar(count, progress_message.c_str());
    }

    std::string s;

    for (uint64_t i = 0; i < count; ++i) {

        ::google::protobuf::io::CodedInputStream *coded_in =
          new ::google::protobuf::io::CodedInputStream(raw_in);

        if (progress) progress->Progressed(i);

        uint32_t msgSize = 0;
        coded_in->ReadVarint32(&msgSize);

        if ((msgSize > 0) &&
            (coded_in->ReadString(&s, msgSize))) {
            Graph o;
            o.ParseFromString(s);
            extend(o);
        }

        delete coded_in;
    }

    delete raw_in;

    if (progress) {
        delete progress;
        progress = NULL;
        cerr << endl;
    }

    //topologically_sort_graph();
    //build_indexes();
}

void VG::serialize_to_ostream(ostream& out, int64_t chunk_size) {

    // for chunks of something
    ::google::protobuf::io::ZeroCopyOutputStream *raw_out =
          new ::google::protobuf::io::OstreamOutputStream(&out);
    ::google::protobuf::io::CodedOutputStream *coded_out =
          new ::google::protobuf::io::CodedOutputStream(raw_out);

    // save the number of the messages to be serialized into the output file
    int64_t count = graph.node_size() / chunk_size + 1;
    coded_out->WriteVarint64(count);

    std::string s;
    uint64_t written = 0;
    for (int64_t n = 0; n < count; ++n) {
        VG g;
        for (int64_t j = n * chunk_size;
             j < (n+1)*chunk_size && j < graph.node_size();
             ++j) {
            Node* node = graph.mutable_node(j);
            node_context(node, g);
        }
        g.graph.SerializeToString(&s);
        coded_out->WriteVarint32(s.size());
        coded_out->WriteRaw(s.data(), s.size()); // ->WriteString(s)
        ++written;
    }

    delete coded_out;
    delete raw_out;

}

VG::~VG(void) {
    destroy_alignable_graph();
}

VG::VG(void) {
    init();
}

void VG::init(void) {
    gssw_aligner = NULL;
    current_id = 1;
    show_progress = false;
    progress_message = "progress";
    progress = NULL;
}

VG::VG(set<Node*>& nodes, set<Edge*>& edges) {
    init();
    add_nodes(nodes);
    add_edges(edges);
    topologically_sort_graph();
}

// check for conflict (duplicate nodes and edges) occurs within add_* functions

void VG::add_nodes(set<Node*>& nodes) {
    for (auto node : nodes) {
        add_node(*node);
    }
}

void VG::add_edges(set<Edge*>& edges) {
    for (auto edge : edges) {
        add_edge(*edge);
    }
}

void VG::add_nodes(vector<Node>& nodes) {
    for (auto& node : nodes) {
        add_node(node);
    }
}

void VG::add_edges(vector<Edge>& edges) {
    for (auto& edge : edges) {
        add_edge(edge);
    }
}


void VG::add_node(Node& node) {
    if (!has_node(node)) {
        Node* new_node = graph.add_node(); // add it to the graph
        new_node->set_sequence(node.sequence());
        new_node->set_id(node.id());
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = graph.node_size()-1;
    }
}

void VG::add_edge(Edge& edge) {
    if (!has_edge(edge)) {
        Edge* new_edge = graph.add_edge(); // add it to the graph
        new_edge->set_from(edge.from());
        new_edge->set_to(edge.to());
        set_edge(edge.from(), edge.to(), new_edge);
        edge_index[new_edge] = graph.edge_size()-1;
    }
}

int64_t VG::node_count(void) {
    return graph.node_size();
}

int64_t VG::edge_count(void) {
    return graph.edge_size();
}

int VG::in_degree(Node* node) {
    auto in = edges_to_from.find(node->id());
    if (in == edges_to_from.end()) {
        return 0;
    } else {
        return in->second.size();
    }
}

int VG::out_degree(Node* node) {
    auto out = edges_from_to.find(node->id());
    if (out == edges_from_to.end()) {
        return 0;
    } else {
        return out->second.size();
    }
}

void VG::edges_of_node(Node* node, vector<Edge*>& edges) {
    hash_map<int64_t, vector<int64_t> >::iterator in = edges_to_from.find(node->id());
    if (in != edges_to_from.end()) {
        vector<int64_t>::iterator e = in->second.begin();
        for (; e != in->second.end(); ++e) {
            Edge* edge = edge_by_id[make_pair(*e, node->id())];
            if (!edge) {
                cerr << "error:[VG::edges_of_node] nonexistent edge " << *e << "->" << node->id() << endl;
                exit(1);
            }
            edges.push_back(edge);
        }
    }
    hash_map<int64_t, vector<int64_t> >::iterator out = edges_from_to.find(node->id());
    if (out != edges_from_to.end()) {
        vector<int64_t>::iterator e = out->second.begin();
        for (; e != out->second.end(); ++e) {
            Edge* edge = edge_by_id[make_pair(node->id(), *e)];
            if (!edge) {
                cerr << "error:[VG::edges_of_node] nonexistent edge " << node->id() << "->" << *e << endl;
                exit(1);
            }
            edges.push_back(edge);
        }
    }
}

void VG::edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges) {
    for (set<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        vector<Edge*> ev;
        edges_of_node(*n, ev);
        for (vector<Edge*>::iterator e = ev.begin(); e != ev.end(); ++e) {
            edges.insert(*e);
        }
    }
}

int64_t VG::total_length_of_nodes(void) {
    int64_t length = 0;
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        length += n->sequence().size();
    }
    return length;
}

void VG::build_indexes(void) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        edge_index[e] = i;
        set_edge(e->from(), e->to(), e);
    }
}

void VG::clear_indexes(void) {
    node_index.clear();
    node_by_id.clear();
    edge_by_id.clear();
    edge_index.clear();
    edges_from_to.clear();
    edges_to_from.clear();
}

void VG::clear_indexes_no_resize(void) {
#ifdef USE_DENSE_HASH
    node_index.clear_no_resize();
    node_by_id.clear_no_resize();
    edge_by_id.clear_no_resize();
    edge_index.clear_no_resize();
    edges_from_to.clear_no_resize();
    edges_to_from.clear_no_resize();
#else
    clear_indexes();
#endif
}

void VG::resize_indexes(void) {
    node_index.resize(graph.node_size());
    node_by_id.resize(graph.node_size());
    edge_by_id.resize(graph.edge_size());
    edge_index.resize(graph.edge_size());
    edges_from_to.resize(graph.edge_size());
    edges_to_from.resize(graph.edge_size());
}

void VG::rebuild_indexes(void) {
    //clear_indexes();
    //resize_indexes();
    clear_indexes_no_resize();
    build_indexes();
}

bool VG::empty(void) {
    return graph.node_size() == graph.edge_size() == 0;
}

bool VG::has_node(Node* node) {
    return node && has_node(node->id());
}

bool VG::has_node(Node& node) {
    return has_node(node.id());
}

bool VG::has_node(int64_t id) {
    return node_by_id.find(id) != node_by_id.end();
}

bool VG::has_edge(Edge* edge) {
    return edge && has_edge(edge->from(), edge->to());
}

bool VG::has_edge(Edge& edge) {
    return has_edge(edge.from(), edge.to());
}

bool VG::has_edge(int64_t from, int64_t to) {
    return edge_by_id.find(make_pair(from, to)) != edge_by_id.end();
}

// remove duplicated nodes and edges that would occur if we merged the graphs
void VG::remove_duplicated_in(VG& g) {
    vector<Node*> nodes_to_destroy;
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (g.has_node(n)) {
            nodes_to_destroy.push_back(n);
        }
    }
    vector<Edge*> edges_to_destroy;
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (g.has_edge(e)) {
            edges_to_destroy.push_back(e);
        }
    }
    for (vector<Node*>::iterator n = nodes_to_destroy.begin();
         n != nodes_to_destroy.end(); ++n) {
        g.destroy_node(g.get_node((*n)->id()));
    }
    for (vector<Edge*>::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        g.destroy_edge(g.get_edge((*e)->from(), (*e)->to()));
    }
}

void VG::merge_union(VG& g) {
    // remove duplicates, then merge
    remove_duplicated_in(g);
    if (g.graph.node_size() > 0) {
        merge(g.graph);
    }
}

void VG::merge(VG& g) {
    merge(g.graph);
}

// this merges without any validity checks
// this could be rather expensive if the graphs to merge are largely overlapping
void VG::merge(Graph& g) {
    graph.mutable_node()->MergeFrom(g.node());
    graph.mutable_edge()->MergeFrom(g.edge());
    rebuild_indexes();
}

// iterates over nodes and edges, adding them in when they don't already exist
void VG::extend(VG& g) {
    for (int64_t i = 0; i < g.graph.node_size(); ++i) {
        Node* n = g.graph.mutable_node(i);
        if (!has_node(n)) {
            add_node(*n);
        }
    }
    for (int64_t i = 0; i < g.graph.edge_size(); ++i) {
        Edge* e = g.graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        }
    }
}

void VG::extend(Graph& graph) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (!has_node(n)) {
            add_node(*n);
        }
    }
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        if (!has_edge(e)) {
            add_edge(*e);
        }
    }
}

// extend this graph by g, connecting the tails of this graph to the heads of the other
// the ids of the second graph are modified for compact representation
void VG::append(VG& g) {

    // compact and increment the ids of g out of range of this graph
    //g.compact_ids();

    // assume we've already compacted the other, or that id compaction doesn't matter
    // just get out of the way
    g.increment_node_ids(max_node_id());

    // get the heads of the other graph, now that we've compacted the ids
    vector<Node*> heads = g.head_nodes();
    // collect ids as node*'s may change
    vector<int64_t> heads_ids;
    for (vector<Node*>::iterator n = heads.begin(); n != heads.end(); ++n) {
        heads_ids.push_back((*n)->id());
    }

    // get the current tails of this graph
    vector<Node*> tails = tail_nodes();
    vector<int64_t> tails_ids;
    for (vector<Node*>::iterator n = tails.begin(); n != tails.end(); ++n) {
        tails_ids.push_back((*n)->id());
    }

    // add in the other graph
    // note that we don't use merge_union because we are ensured non-overlapping ids
    merge(g);

    /*
    cerr << "this graph size " << node_count() << " nodes " << edge_count() << " edges" << endl;
    cerr << "in append with " << heads.size() << " heads and " << tails.size() << " tails" << endl;
    */

    // now join the tails to heads
    for (vector<int64_t>::iterator t = tails_ids.begin(); t != tails_ids.end(); ++t) {
        for (vector<int64_t>::iterator h = heads_ids.begin(); h != heads_ids.end(); ++h) {
            create_edge(*t, *h);
        }
    }
}

void VG::combine(VG& g) {
    // compact and increment the ids of g out of range of this graph
    //g.compact_ids();
    g.increment_node_ids(max_node_id());
    // now add it into the current graph, without connecting any nodes
    extend(g);
}

int64_t VG::max_node_id(void) {
    int64_t max_id = 0;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (n->id() > max_id) {
            max_id = n->id();
        }
    }
    return max_id;
}

void VG::compact_ids(void) {
    hash_map<int64_t, int64_t> new_id;
    int64_t id = 1; // start at 1
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        new_id[n->id()] = id++;
    }
//#pragma omp parallel for
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        n->set_id(new_id[n->id()]);
    }
//#pragma omp parallel for
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        e->set_from(new_id[e->from()]);
        e->set_to(new_id[e->to()]);
    }
    rebuild_indexes();
}

void VG::increment_node_ids(int64_t increment) {
    for_each_node_parallel([increment](Node* n) {
            n->set_id(n->id()+increment);
        });
    for_each_edge_parallel([increment](Edge* e) {
            e->set_from(e->from()+increment);
            e->set_to(e->to()+increment);
        });
    rebuild_indexes();
}

void VG::decrement_node_ids(int64_t decrement) {
    increment_node_ids(-decrement);
}

void VG::swap_node_id(int64_t node_id, int64_t new_id) {
    swap_node_id(node_by_id[node_id], new_id);
}

void VG::swap_node_id(Node* node, int64_t new_id) {

    //cerr << "swapping " << node->id() << " for new id " << new_id << endl;
    int edge_n = edge_count();
    int64_t old_id = node->id();
    node->set_id(new_id);
    node_by_id.erase(old_id);

    // we check if the old node exists, and bail out if we're not doing what we expect
    assert(node_by_id.find(new_id) == node_by_id.end());

    // otherwise move to a new id
    node_by_id[new_id] = node;

    set<pair<int64_t, int64_t> > edges_to_destroy;
    set<pair<int64_t, int64_t> > edges_to_create;

    vector<int64_t>& to = edges_to(old_id);
    for (vector<int64_t>::iterator e = to.begin(); e != to.end(); ++e) {
        edges_to_create.insert(make_pair(*e, new_id));
        edges_to_destroy.insert(make_pair(*e, old_id));
    }

    vector<int64_t>& from = edges_from(old_id);
    for (vector<int64_t>::iterator e = from.begin(); e != from.end(); ++e) {
        edges_to_create.insert(make_pair(new_id, *e));
        edges_to_destroy.insert(make_pair(old_id, *e));
    }

    assert(edges_to_destroy.size() == edges_to_create.size());

    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        destroy_edge(e->first, e->second);
    }

    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_create.begin();
         e != edges_to_create.end(); ++e) {
        create_edge(e->first, e->second);
    }

    assert(edge_n == edge_count());

    // we maintain a valid graph
    // this an expensive check but should work (for testing only)
    //assert(is_valid());

}



// construct from VCF records
// --------------------------
// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)
//
VG::VG(vector<vcf::Variant>& records, string seq, string chrom, int offset) {
    init();
    from_vcf_records(&records, seq, chrom, offset);
}

void VG::from_vcf_records(vector<vcf::Variant>* r, string seq, string chrom, int offset) {

    //init();
    vector<vcf::Variant>& records = *r;

/*
    int tid = omp_get_thread_num();
#pragma omp critical
    {
        cerr << tid << ": in from_vcf_records" << endl;
        cerr << tid << ": with " << records.size() << " vars" << endl;
        cerr << tid << ": and " << seq.size() << "bp" << endl;
    }
*/


    map<long, Node*> reference_path;
    // track the last nodes so that we can connect everything
    // completely when variants occur in succession
    map<long, set<Node*> > nodes_by_end_position;
    map<long, set<Node*> > nodes_by_start_position;

    // parsed alternates
    map<long, set<vcf::VariantAllele> > altp;

    Node* ref_node = create_node(seq);
    reference_path[0] = ref_node;

    for (auto& var : records) {

        // adjust the variant position relative to the sequence
        var.position -= offset;

//#pragma omp critical
//        cerr << omp_get_thread_num() << ": " << var << endl;

        // decompose the alt
        bool flat_input_vcf = false; // hack
        map<string, vector<vcf::VariantAllele> > alternates
            = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());

        for (auto& alleles : alternates) {
            for (auto& allele : alleles.second) {
                altp[allele.position].insert(allele);
            }
        }
    }

    for (auto& va : altp) {

        set<vcf::VariantAllele>& alleles = va.second;

        for (auto allele : alleles) {

            // reference alleles are provided naturally by the reference itself
            if (allele.ref == allele.alt) {
                continue;
            }

            // 0/1 based conversion happens in offset
            long allele_start_pos = allele.position;
            long allele_end_pos = allele_start_pos + allele.ref.size();
            // for ordering, set insertion start position at +1
            // otherwise insertions at the same position will loop infinitely
            //if (allele_start_pos == allele_end_pos) allele_end_pos++;

            if (allele_start_pos == 0) {
                // ensures that we can handle variation at first position
                // (important when aligning)
                Node* root = create_node("");
                reference_path[-1] = root;
                nodes_by_start_position[-1].insert(root);
                nodes_by_end_position[0].insert(root);
            }

            Node* left_ref_node = NULL;
            Node* middle_ref_node = NULL;
            Node* right_ref_node = NULL;

            divide_path(reference_path,
                        allele_start_pos,
                        left_ref_node,
                        right_ref_node);

            /*
              cerr << "nodes: left: " << left_ref_node->id() << ":" << left_ref_node->sequence()
              << " right: " << right_ref_node->id() << ":" << right_ref_node->sequence() << endl;
            */

            // if the ref portion of the allele is not empty, then we need to make another cut
            if (!allele.ref.empty()) {
                divide_path(reference_path,
                            allele_end_pos,
                            middle_ref_node,
                            right_ref_node);
            }

            Node* alt_node = NULL;
            // create a new alt node and connect the pieces from before
            if (!allele.alt.empty() && !allele.ref.empty()) {
                //cerr << "both alt and ref have sequence" << endl;

                alt_node = create_node(allele.alt);
                create_edge(left_ref_node, alt_node);
                create_edge(alt_node, right_ref_node);

                nodes_by_end_position[allele_end_pos].insert(alt_node);
                nodes_by_end_position[allele_end_pos].insert(middle_ref_node);
                //nodes_by_end_position[allele_start_pos].insert(left_ref_node);
                nodes_by_start_position[allele_start_pos].insert(alt_node);
                nodes_by_start_position[allele_start_pos].insert(middle_ref_node);

            } else if (!allele.alt.empty()) { // insertion

                //cerr << "alt has sequence" << endl;
                alt_node = create_node(allele.alt);
                create_edge(left_ref_node, alt_node);
                create_edge(alt_node, right_ref_node);
                nodes_by_end_position[allele_end_pos].insert(alt_node);
                nodes_by_end_position[allele_end_pos].insert(left_ref_node);
                nodes_by_start_position[allele_start_pos].insert(alt_node);

            } else {// otherwise, we have a deletion

                    //cerr << "ref has sequence" << endl;
                create_edge(left_ref_node, right_ref_node);
                nodes_by_end_position[allele_end_pos].insert(left_ref_node);
                nodes_by_start_position[allele_start_pos].insert(left_ref_node);

            }

            /*
            cerr << allele << endl;
            if (left_ref_node) cerr << "left_ref " << left_ref_node->id()
                                    << " " << left_ref_node->sequence() << endl;
            if (middle_ref_node) cerr << "middle_ref " << middle_ref_node->id()
                                      << " " << middle_ref_node->sequence() << endl;
            if (alt_node) cerr << "alt_node " << alt_node->id()
                               << " " << alt_node->sequence() << endl;
            if (right_ref_node) cerr << "right_ref " << right_ref_node->id()
                                     << " " << right_ref_node->sequence() << endl;
            */


            if (allele_end_pos == seq.size()) {
                // ensures that we can handle variation at last position (important when aligning)
                Node* end = create_node("");
                reference_path[allele_end_pos] = end;
                // for consistency, this should be handled below in the start/end connections
                if (alt_node) {
                    create_edge(alt_node, end);
                }
                if (middle_ref_node) {
                    create_edge(middle_ref_node, end);
                }
            }

            //print_edges();
            /*
            if (!is_valid()) {
                cerr << "graph is invalid after variant " << *a << endl;
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
            */

        }

        map<long, set<Node*> >::iterator ep
            = nodes_by_end_position.find(va.first);
        map<long, set<Node*> >::iterator sp
            = nodes_by_start_position.find(va.first);
        if (ep != nodes_by_end_position.end()
            && sp != nodes_by_start_position.end()) {
            set<Node*>& previous_nodes = ep->second;
            set<Node*>& current_nodes = sp->second;
            for (set<Node*>::iterator n = previous_nodes.begin();
                 n != previous_nodes.end(); ++n) {
                for (set<Node*>::iterator m = current_nodes.begin();
                     m != current_nodes.end(); ++m) {
                    if (node_index.find(*n) != node_index.end()
                        && node_index.find(*m) != node_index.end()
                        && !(previous_nodes.count(*n) && current_nodes.count(*n)
                             && previous_nodes.count(*m) && current_nodes.count(*m))
                        ) {
                        /*
                        cerr << "connecting previous "
                             << (*n)->id() << " @end=" << ep->first << " to current "
                             << (*m)->id() << " @start=" << sp->first << endl;
                        */
                        create_edge(*n, *m);
                    }
                }
            }
        }

        // clean up previous
        while (!nodes_by_end_position.empty() && nodes_by_end_position.begin()->first < va.first) {
            nodes_by_end_position.erase(nodes_by_end_position.begin()->first);
        }

        while (!nodes_by_start_position.empty() && nodes_by_start_position.begin()->first < va.first) {
            nodes_by_start_position.erase(nodes_by_start_position.begin()->first);
        }

    }

    topologically_sort_graph();
    compact_ids();

}

void VG::print_edges(void) {
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        int64_t f = e->from();
        int64_t t = e->to();
        cerr << f << "->" << t << " ";
    }
    cerr << endl;
}

bool allATGC(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
}

VG::VG(vcf::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target_region,
       int vars_per_region,
       bool showprog) {

    init();

    omp_set_dynamic(1); // use dynamic scheduling

    show_progress = showprog;

    map<string, VG*> refseq_graph;

    vector<string> targets;
    if (!target_region.empty()) {
        targets.push_back(target_region);
    } else {
        for (vector<string>::iterator r = reference.index->sequenceNames.begin();
             r != reference.index->sequenceNames.end(); ++r) {
            targets.push_back(*r);
        }
    }

    // to scale up, we have to avoid big string memcpys
    // this could be accomplished by some deep surgery on the construction routines
    // however, that could be a silly thing to do,
    // because why break something that's conceptually clear
    // and anyway, we want to break the works into chunks
    //
    // there is a development that could be important
    // our chunk size isn't going to reach into the range where we'll have issues (>several megs)
    // so we'll run this for regions of moderate size, scaling up in the case that we run into a big deletion
    // 
    // 

    for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {

        //string& seq_name = *t;
        string seq_name;
        string target = *t;
        int start_pos = 0, stop_pos = 0;
        // nasty hack for handling single regions
        parse_region(target,
                     seq_name,
                     start_pos,
                     stop_pos);
        if (stop_pos > 0) {
            variantCallFile.setRegion(seq_name, start_pos, stop_pos);
        } else {
            variantCallFile.setRegion(seq_name);
            stop_pos = reference.sequenceLength(seq_name);
        }
        vcf::Variant var(variantCallFile);

        vector<vcf::Variant>* region = NULL;

        // convert from 1-based input to 0-based internal format
        // and handle the case where we are already doing the whole chromosome
        int64_t start = start_pos ? start_pos - 1 : 0;
        int64_t end = start;
        bool done_with_chrom = false;
        // track if the variant we are looking at has the 3'-most reference extent of any variant in the bunch
        bool var_is_at_end = false;
        struct Plan {
            VG* graph;
            vector<vcf::Variant>* vars;
            string seq;
            string name;
            int offset;
        };
        deque<Plan*> construction;

        // our graph for this refseq
        refseq_graph[target] = new VG;
        // and the construction queue
        list<VG*> graphq;
        omp_lock_t appending_graphs;
        omp_init_lock(&appending_graphs);

        // for tracking progress through the chromosome
        map<VG*, unsigned long> graph_end;
        string message = "constructing graph for " + seq_name;
        ProgressBar* progress = NULL;
        if (show_progress) {
            progress = new ProgressBar(stop_pos-start_pos, message.c_str());
        }
        if (progress) progress->Progressed(0);

        set<VG*> graph_completed;

        // omp pragma here to define parallel section

        // this system is not entirely general
        // there will be a problem when the regions of overlapping deletions become too large
        // then the inter-dependence of each region will make parallel construction in this way difficult
        // because the chunks will get too large

#pragma omp parallel default(none)                                      \
    shared(vars_per_region, region, target, stop_pos,                   \
           variantCallFile, done_with_chrom, reference,                 \
           seq_name, start, start_pos, var_is_at_end, progress,         \
           refseq_graph, graphq, appending_graphs,                      \
           end, var, cerr, construction, graph_completed, graph_end)

        while (!done_with_chrom || !construction.empty() || !graphq.empty()) {

            int tid = omp_get_thread_num();

/*
#pragma omp critical (debug)
            cerr << tid << ": " << "at top of loop "
                 << (done_with_chrom ? "done" : "more") << " "
                 << construction.size() << " " << graphq.size() << endl;
*/

// each thread should check if there is a new item in the queue to remove
// and append

            VG* g = refseq_graph[target];

            if (!graphq.empty()
                && graph_completed.count(graphq.front())
                && omp_test_lock(&appending_graphs)) {
                // this thread has the lock
                while (!graphq.empty() && graph_completed.count(graphq.front())) {
                    VG* o = NULL;
#pragma omp critical (graphq)
                    {
                        o = graphq.front();
                        graphq.pop_front();
                        graph_completed.erase(o);
                    }
                    g->append(*o);
                    delete o;
#pragma omp critical (progress)
                    if (progress) progress->Progressed(graph_end[o]-start_pos);
                }
                // reset lock
                omp_unset_lock(&appending_graphs);
            }

            // processing of VCF file should only be handled by one thread at a time
#pragma omp critical (vcf_input)
            if (!done_with_chrom) {

                if (!region) {
                    region = new vector<vcf::Variant>;
                }

                done_with_chrom = !variantCallFile.getNextVariant(var);
                var.position -= 1; // convert to 0-based

                vcf::Variant* lvar = (region->empty() ? NULL : &region->back());
//#pragma omp critical (debug)
//                cerr << tid << ": got variant = " << var << endl;

                // skip non-DNA sequences, such as SVs
                bool isDNA = allATGC(var.ref);
                for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                    if (!allATGC(*a)) isDNA = false;
                }

                if (!done_with_chrom && isDNA) {
                    if (lvar) {
                        end = max(end, (int64_t) (lvar->position + lvar->ref.size()));
                        if (var.position <= end) {
                            var_is_at_end = false;
                        } else {
                            var_is_at_end = true;
                        }
                    }
                    region->push_back(var);
                }

                if (done_with_chrom) {
                    var_is_at_end = true;
                }

                if (region
                    && (region->size() - 1 >= vars_per_region
                        && var_is_at_end
                        || done_with_chrom)) {

                    vector<vcf::Variant>* vars = region;
                    region = NULL;
                
                    // find the end of the region
                    if (done_with_chrom) {
                        end = stop_pos;
                    }


                    if (!done_with_chrom) {
                        // push an empty list back to handle the next region
                        region = new vector<vcf::Variant>;
                        region->push_back(vars->back());
                        vars->pop_back();
                    }

                    // makes a new construction plan
                    Plan* plan = new Plan;
                    plan->graph = new VG;

                    // retain reference to graph
#pragma omp critical (graphq)
                    graphq.push_back(plan->graph);
                    // variants
                    plan->vars = vars;
                    // reference sequence
                    plan->seq = reference.getSubSequence(seq_name, start, end - start);
                    // chromosome name
                    plan->name = seq_name;
                    // offset
                    plan->offset = start;
                    // store the plan in our construction queue
                    construction.push_front(plan);
                    // record our end position (for progress logging only)
                    if (progress) graph_end[plan->graph] = end;

                    // this start is the next end
                    start = end;
                    // and reset end
                    end = 0;
                    // reset var_is_at_end
                    var_is_at_end = false;

                }
            }

            // execute this as a parallel block
            // it can be done in an entirely separate thread
            if (!construction.empty()) {
                Plan* plan = NULL;
#pragma omp critical (construction)
                {
                    if (!construction.empty()) {
                        plan = construction.back();
                        construction.pop_back();
                    }
                }
                
                if (plan) {

/*
#pragma omp critical (debug)
                    cerr << tid << ": " << "constructing graph " << plan->graph << " over "
                         << plan->vars->size() << " variants in " <<plan->seq.size() << "bp "
                         << plan->name << ":" << plan->offset
                         << "-" << plan->offset + plan->seq.size() << endl;
*/

                    plan->graph->from_vcf_records(plan->vars,
                                                  plan->seq,
                                                  plan->name,
                                                  plan->offset);

#pragma omp critical (graphq)
                    {
                        graph_completed.insert(plan->graph);
//#pragma omp critical (debug)
//                        cerr << tid << ": " << "completed graph " << plan->graph << endl;
                    }
                    // TODO plan should probably be a class, so that this can be managed
                    delete plan->vars;
                    delete plan;
                }
            }

            // if we are waiting on the completion of graph appending or construction, slow down
            if (done_with_chrom && (!construction.empty() || !graphq.empty())) {
                usleep(10000); // 10ms
            }
        }
        // parallel end

        // clean up "null" nodes that are used for maintaining structure between temporary subgraphs
        refseq_graph[target]->remove_null_nodes_forwarding_edges();
        // then use topological sorting and re-compression of the id space to make sure that
        refseq_graph[target]->topologically_sort_graph();
        // we get identical graphs no matter what the region size is
        refseq_graph[target]->compact_ids();

        if (progress) {
            delete progress;
            progress = NULL;
            cerr << endl;
        }

        omp_destroy_lock(&appending_graphs);

    }

    // small hack for efficiency when constructing over a single chromosome
    if (refseq_graph.size() == 1) {
        *this = *refseq_graph[targets.front()];
    } else {
        // where we have multiple targets
        for (vector<string>::iterator t = targets.begin(); t != targets.end(); ++t) {
            // merge the variants into one graph
            VG& g = *refseq_graph[*t];
            combine(g);
        }
    }
}

void VG::topologically_sort_graph(void) {
    deque<Node*> sorted_nodes;
    topological_sort(sorted_nodes);
    deque<Node*>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < graph.node_size() && n != sorted_nodes.end();
          ++i, ++n) {
        swap_nodes(graph.mutable_node(i), *n);
    }
}

void VG::swap_nodes(Node* a, Node* b) {
    int aidx = node_index[a];
    int bidx = node_index[b];
    graph.mutable_node()->SwapElements(aidx, bidx);
    node_index[a] = bidx;
    node_index[b] = aidx;
}

Edge* VG::create_edge(Node* from, Node* to) {
    return create_edge(from->id(), to->id());
}

Edge* VG::create_edge(int64_t from, int64_t to) {
    //cerr << "creating edge " << from << "->" << to << endl;
    // prevent self-linking (violates DAG/partial ordering property)
    if (to == from) return NULL;
    // ensure the edge does not already exist
    Edge* edge = get_edge(from, to);
    if (edge) return edge;
    // if not, create it
    edge = graph.add_edge();
    edge->set_from(from);
    edge->set_to(to);
    set_edge(from, to, edge);
    edge_index[edge] = graph.edge_size()-1;
    //cerr << "created edge " << edge->from() << "->" << edge->to() << endl;
    return edge;
}

Edge* VG::get_edge(int64_t from, int64_t to) {
    pair_hash_map<pair<int64_t, int64_t>, Edge*>::iterator e = edge_by_id.find(make_pair(from, to));
    if (e != edge_by_id.end()) {
        return e->second;
    } else {
        return NULL;
    }
}

void VG::set_edge(int64_t from, int64_t to, Edge* edge) {
    if (!has_edge(edge)) {
        edge_by_id[make_pair(from, to)] = edge;
        edges_from_to[from].push_back(to);
        edges_to_from[to].push_back(from);
    }
}

void VG::for_each_edge_parallel(function<void(Edge*)> lambda) {
#pragma omp parallel for
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        lambda(graph.mutable_edge(i));
    }
}

void VG::for_each_edge(function<void(Edge*)> lambda) {
    for (int64_t i = 0; i < graph.edge_size(); ++i) {
        lambda(graph.mutable_edge(i));
    }
}

vector<int64_t>& VG::edges_from(int64_t id) {
    hash_map<int64_t, vector<int64_t> >::iterator e = edges_from_to.find(id);
    if (e == edges_from_to.end()) {
        return empty_ids;
    } else {
        return e->second;
    }
}

vector<int64_t>& VG::edges_to(int64_t id) {
    hash_map<int64_t, vector<int64_t> >::iterator e = edges_to_from.find(id);
    if (e == edges_to_from.end()) {
        return empty_ids;
    } else {
        return e->second;
    }
}

vector<int64_t>& VG::edges_from(Node* node) {
    return edges_from(node->id());
}

vector<int64_t>& VG::edges_to(Node* node) {
    return edges_to(node->id());
}

void VG::remove_edge_fti(int64_t from, int64_t to) {
    vector<int64_t>& f = edges_from(from);
    swap_remove(f, to);
}

void VG::remove_edge_tfi(int64_t from, int64_t to) {
    vector<int64_t>& t = edges_to(to);
    swap_remove(t, from);
}

void VG::destroy_edge(int64_t from, int64_t to) {
    destroy_edge(get_edge(from, to));
}

void VG::destroy_edge(Edge* edge) {
    //cerr << "destroying edge " << edge->from() << "->" << edge->to() << endl;

    // noop on NULL pointer or non-existent edge
    if (!has_edge(edge)) { return; }
    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }
    // erase from indexes

    //cerr << "erasing from indexes" << endl;

    remove_edge_fti(edge->from(), edge->to());
    remove_edge_tfi(edge->from(), edge->to());

    //assert(edges_from_to[edge->from()].find(edge->to()) == edges_from_to[edge->from()].end());
    //assert(edges_to_from[edge->to()].find(edge->from()) == edges_to_from[edge->to()].end());

    // removing the sub-indexes if they are now empty
    // we must do this to maintain a valid structure
    if (edges_from_to[edge->from()].empty()) edges_from_to.erase(edge->from());
    if (edges_to_from[edge->to()].empty()) edges_to_from.erase(edge->to());

    // erase from edges by moving to end and dropping

    // get the last edge index (lei) and this edge index (tei)
    int lei = graph.edge_size()-1;
    int tei = edge_index[edge];

    // erase this edge from the index
    // we'll fix up below
    edge_index.erase(edge);

    // Why do we check that lei != tei?
    //
    // It seems, after an inordinate amount of testing and probing,
    // that if we call erase twice on the same entry, we'll end up corrupting the hash_map
    //
    // So, if the element is already at the end of the table,
    // take a little break and just remove the last edge in graph

    // if we need to move the element to the last position in the array...
    if (lei != tei) {

        // get a pointer to the last element
        Edge* last = graph.mutable_edge(lei);

        // erase from our index
        edge_index.erase(last);

        // swap
        graph.mutable_edge()->SwapElements(tei, lei);

        // point to new position
        Edge* nlast = graph.mutable_edge(tei);

        // insert the new edge index position
        edge_index[nlast] = tei;

        // and fix edge indexes for moved edge object
        set_edge(nlast->from(), nlast->to(), nlast);

    }

    // drop the last position, erasing the node
    graph.mutable_edge()->RemoveLast();

    //if (!is_valid()) { cerr << "graph ain't valid" << endl; }

}

Node* VG::get_node(int64_t id) {
    hash_map<int64_t, Node*>::iterator n = node_by_id.find(id);
    if (n != node_by_id.end()) {
        return n->second;
    } else {
        // again... should this throw an error?
        return NULL;
    }
}

// use the VG class to generate ids
Node* VG::create_node(string seq) {
    // create the node
    Node* node = graph.add_node();
    node->set_sequence(seq);
    node->set_id(current_id++);
    //cerr << "creating node " << node->id() << endl;
    // copy it into the graph
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = graph.node_size()-1;
    //if (!is_valid()) cerr << "graph invalid" << endl;
    return node;
}

void VG::for_each_node_parallel(function<void(Node*)> lambda) {
    if (show_progress) {
        progress = new ProgressBar(graph.node_size(), progress_message.c_str());
    }
    int64_t completed = 0;
#pragma omp parallel for shared(completed)
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        lambda(graph.mutable_node(i));
        if (progress && completed++ % 1000 == 0) {
#pragma omp critical (progress_bar)
            progress->Progressed(completed);
        }
    }
    if (show_progress) {
        progress->Progressed(completed);
        delete progress;
        progress = NULL;
        cerr << endl;
    }
}

void VG::for_each_node(function<void(Node*)> lambda) {
    for (int64_t i = 0; i < graph.node_size(); ++i) {
        lambda(graph.mutable_node(i));
    }
}

// a graph composed of this node and its edges
void VG::node_context(Node* node, VG& g) {
    // add the node
    g.add_node(*node);
    // and its edges
    vector<int64_t>& to = edges_to(node->id());
    for (vector<int64_t>::iterator e = to.begin(); e != to.end(); ++e) {
        g.add_edge(*get_edge(*e, node->id()));
    }
    vector<int64_t>& from = edges_from(node->id());
    for (vector<int64_t>::iterator e = from.begin(); e != from.end(); ++e) {
        g.add_edge(*get_edge(node->id(), *e));
    }
}

void VG::destroy_node(int64_t id) {
    destroy_node(get_node(id));
}

void VG::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
    //cerr << "destroying node " << node->id() << endl;
    // noop on NULL/nonexistent node
    if (!has_node(node)) { return; }
    // remove edges associated with node
    set<pair<int64_t, int64_t> > edges_to_destroy;
    hash_map<int64_t, vector<int64_t> >::iterator e = edges_from_to.find(node->id());
    if (e != edges_from_to.end()) {
        for (vector<int64_t>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(make_pair(node->id(), *f));
        }
    }
    e = edges_to_from.find(node->id());
    if (e != edges_to_from.end()) {
        for (vector<int64_t>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(make_pair(*f, node->id()));
        }
    }
    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        destroy_edge(e->first, e->second);
    }
    // assert cleanup
    edges_to_from.erase(node->id());
    edges_from_to.erase(node->id());

    // swap node with the last in nodes
    // call RemoveLast() to drop the node
    int lni = graph.node_size()-1;
    int tni = node_index[node];

    if (lni != tni) {
        // swap this node with the last one
        Node* last = graph.mutable_node(lni);
        graph.mutable_node()->SwapElements(tni, lni);
        Node* nlast = graph.mutable_node(tni);

        // and fix up the indexes
        node_by_id[last->id()] = nlast;
        node_index.erase(last);
        node_index[nlast] = tni;
    }

    // remove this node (which is now the last one) and remove references from the indexes
    node_by_id.erase(node->id());
    node_index.erase(node);
    graph.mutable_node()->RemoveLast();
    //if (!is_valid()) { cerr << "graph is invalid after destroy_node" << endl; exit(1); }
}

void VG::remove_null_nodes(void) {
    vector<int64_t> to_remove;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        if (node->sequence().size() == 0) {
            to_remove.push_back(node->id());
        }
    }
    for (vector<int64_t>::iterator n = to_remove.begin(); n != to_remove.end(); ++n) {
        destroy_node(*n);
    }
}

void VG::remove_null_nodes_forwarding_edges(void) {
    vector<Node*> to_remove;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        if (node->sequence().size() == 0) {
            to_remove.push_back(node);
        }
    }
    for (vector<Node*>::iterator n = to_remove.begin(); n != to_remove.end(); ++n) {
        remove_node_forwarding_edges(*n);
    }
}

void VG::remove_node_forwarding_edges(Node* node) {
    vector<int64_t>& to = edges_to(node);
    vector<int64_t>& from = edges_from(node);
    // for edge to
    set<pair<int64_t, int64_t> > edges_to_create;
    for (vector<int64_t>::iterator f = to.begin(); f != to.end(); ++f) {
        for (vector<int64_t>::iterator t = from.begin(); t != from.end(); ++t) {
            // connect
            edges_to_create.insert(make_pair(*f, *t));
        }
    }
    for (set<pair<int64_t, int64_t> >::iterator e = edges_to_create.begin();
         e != edges_to_create.end(); ++e) {
        create_edge(e->first, e->second);
    }
    destroy_node(node);
}

// utilities
void VG::divide_node(Node* node, int pos, Node*& left, Node*& right) {

    //cerr << "dividing node " << node->id() << endl;

#pragma omp critical
    if (pos < 0) {
        cerr << omp_get_thread_num() << ": in divide_node "
             << "position (" << pos << ") is less than 0!" << endl
             << "cannot divide node " << node->id() << ":" << node->sequence() << endl;
        exit(1);
    }
/*
#pragma omp critical
    cerr << omp_get_thread_num() << ": in divide_node " << pos << " of " << node->sequence().size() << endl;
*/

    // make our left node
    left = create_node(node->sequence().substr(0,pos));

    hash_map<int64_t, vector<int64_t> >::const_iterator e;
    set<pair<int64_t, int64_t> > edges_to_create;

    // replace node connections to prev (left)
    e = edges_to_from.find(node->id());
    if (e != edges_to_from.end()) {
        for (vector<int64_t>::const_iterator p = e->second.begin();
             p != e->second.end(); ++p) {
            edges_to_create.insert(make_pair(*p, left->id()));
        }
    }

    // make our right node
    right = create_node(node->sequence().substr(pos,node->sequence().size()-1));

    // replace node connections to next (right)
    e = edges_from_to.find(node->id());
    if (e != edges_from_to.end()) {
        for (vector<int64_t>::const_iterator n = e->second.begin();
             n != e->second.end(); ++n) {
            edges_to_create.insert(make_pair(right->id(), *n));
        }
    }

    // create the edges here as otherwise we will invalidate the iterators
    for (set<pair<int64_t, int64_t> >::iterator c = edges_to_create.begin();
         c != edges_to_create.end(); ++c) {
        create_edge(c->first, c->second);
    }

    // connect left to right
    create_edge(left, right);

    destroy_node(node);

}

// for dividing a path of nodes with an underlying coordinate system
void VG::divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {

    map<long, Node*>::iterator target = path.upper_bound(pos);
    --target; // we should now be pointing to the target ref node

    long node_pos = target->first;
    Node* old = target->second;
    
    // nothing to do
    if (node_pos == pos) {

        map<long, Node*>::iterator n = target; --n;
        left = n->second;
        right = target->second;

    } else {

        // divide the target node at our pos
        int diff = pos - node_pos;

        divide_node(old, diff, left, right);

        // left
        path[node_pos] = left;

        // right
        path[pos] = right;
    }
}

void VG::nodes_prev(Node* node, vector<Node*>& nodes) {
    vector<int64_t>& from = edges_to(node);
    for (vector<int64_t>::iterator f = from.begin(); f != from.end(); ++f) {
        nodes.push_back(node_by_id[*f]);
    }
}

void VG::nodes_next(Node* node, vector<Node*>& nodes) {
    vector<int64_t>& to = edges_from(node);
    for (vector<int64_t>::iterator t = to.begin(); t != to.end(); ++t) {
        nodes.push_back(node_by_id[*t]);
    }
}

void VG::prev_kpaths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths) {
    if (length == 0) { return; }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    postfix.push_front(node);
    vector<Node*> prev_nodes;
    nodes_prev(node, prev_nodes);
    if (prev_nodes.empty()) {
        list<Node*> new_path = postfix;
        paths.insert(new_path);
    } // implicit else
    for (vector<Node*>::iterator p = prev_nodes.begin(); p != prev_nodes.end(); ++p) {
        if ((*p)->sequence().size() < length) {
            prev_kpaths_from_node(*p, length - (*p)->sequence().size(), postfix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = postfix;
            new_path.push_front(*p);
            paths.insert(new_path);
        }
    }
}

void VG::next_kpaths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths) {
    if (length == 0) { return; }
    // start at node
    // do a leftward DFS up to length limit to establish paths from the left of the node
    prefix.push_back(node);
    vector<Node*> next_nodes;
    nodes_next(node, next_nodes);
    if (next_nodes.empty()) {
        list<Node*> new_path = prefix;
        paths.insert(new_path);
    } // implicit else
    for (vector<Node*>::iterator n = next_nodes.begin(); n != next_nodes.end(); ++n) {
        if ((*n)->sequence().size() < length) {
            next_kpaths_from_node(*n, length - (*n)->sequence().size(), prefix, paths);
        } else {
            // create a path for this node
            list<Node*> new_path = prefix;
            new_path.push_back(*n);
            paths.insert(new_path);
        }
    }
}

// iterate over the kpaths in the graph, doing something

void VG::for_each_kpath(int k, function<void(list<Node*>&)> lambda) {
    auto by_node = [k, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, lambda);
    };
    for_each_node(by_node);
}

void VG::for_each_kpath(int k, function<void(Path&)> lambda) {
    auto by_node = [k, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, lambda);
    };
    for_each_node(by_node);
}

// parallel versions of above
// this isn't by default because the lambda may have side effects
// that need to be guarded explicitly

void VG::for_each_kpath_parallel(int k, function<void(list<Node*>&)> lambda) {
    auto by_node = [k, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, lambda);
    };
    for_each_node_parallel(by_node);
}

void VG::for_each_kpath_parallel(int k, function<void(Path&)> lambda) {
    auto by_node = [k, &lambda, this](Node* node) {
        for_each_kpath_of_node(node, k, lambda);
    };
    for_each_node_parallel(by_node);
}

// per-node kpaths

void VG::for_each_kpath_of_node(Node* n, int k,
                                function<void(Path&)> lambda) {
    auto apply_to_path = [&lambda, this](list<Node*>& p) {
        Path path = create_path(p);
        lambda(path);
    };
    for_each_kpath_of_node(n, k, apply_to_path);
}

void VG::for_each_kpath_of_node(Node* node, int k,
                                function<void(list<Node*>&)> lambda) {
    // get left, then right
    set<list<Node*> > prev_paths;
    set<list<Node*> > next_paths;
    list<Node*> empty_list;
    prev_kpaths_from_node(node, k, empty_list, prev_paths);
    next_kpaths_from_node(node, k, empty_list, next_paths);
    // now take the cross and give to the callback
    for (set<list<Node*> >::iterator p = prev_paths.begin(); p != prev_paths.end(); ++p) {
        for (set<list<Node*> >::iterator n = next_paths.begin(); n != next_paths.end(); ++n) {
            list<Node*> path = *p;
            list<Node*>::const_iterator m = n->begin(); ++m; // skips current node, which is included in *p
            while (m != n->end()) {
                path.push_back(*m);
                ++m;
            }
            lambda(path);
        }
    }
}

void VG::kpaths_of_node(Node* node, set<list<Node*> >& paths, int length) {
    auto collect_path = [&paths](list<Node*>& path) {
        paths.insert(path);
    };
    for_each_kpath_of_node(node, length, collect_path);
}

void VG::kpaths_of_node(Node* node, vector<Path>& paths, int length) {
    set<list<Node*> > unique_paths;
    kpaths_of_node(node, unique_paths, length);
    for (set<list<Node*> >::iterator p = unique_paths.begin(); p != unique_paths.end(); ++p) {
        Path path = create_path(*p);
        paths.push_back(path);
    }
}

// aggregators, when a callback won't work

void VG::kpaths(set<list<Node*> >& paths, int length) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* node = graph.mutable_node(i);
        kpaths_of_node(node, paths, length);
    }
}

void VG::kpaths(vector<Path>& paths, int length) {
    set<list<Node*> > unique_paths;
    kpaths(unique_paths, length);
    for (set<list<Node*> >::iterator p = unique_paths.begin(); p != unique_paths.end(); ++p) {
        Path path = create_path(*p);
        paths.push_back(path);
    }
}

// path utilities
// these are in this class because attributes of the path (such as its sequence) are a property of the graph

Path VG::create_path(const list<Node*>& nodes) {
    Path path;
    for (list<Node*>::const_iterator n = nodes.begin(); n != nodes.end(); ++n) {
        Mapping* mapping = path.add_mapping();
        mapping->set_node_id((*n)->id());
    }
    return path;
}

string VG::path_string(const list<Node*>& nodes) {
    string seq;
    for (list<Node*>::const_iterator n = nodes.begin(); n != nodes.end(); ++n) {
        seq.append((*n)->sequence());
    }
    return seq;
}

string VG::path_string(Path& path) {
    string seq;
    for (int i = 0; i < path.mapping_size(); ++i) {
        Mapping* m = path.mutable_mapping(i);
        Node* n = node_by_id[m->node_id()];
        seq.append(n->sequence());
    }
    return seq;
}

void VG::expand_path(const list<Node*>& path, vector<Node*>& expanded) {
    for (list<Node*>::const_iterator n = path.begin(); n != path.end(); ++n) {
        Node* node = *n;
        int s = node->sequence().size();
        for (int i = 0; i < s; ++i) {
            expanded.push_back(node);
        }
    }
}

void VG::node_starts_in_path(const list<Node*>& path,
                             map<Node*, int>& node_start) {
    int i = 0;
    for (list<Node*>::const_iterator n = path.begin(); n != path.end(); ++n) {
        node_start[*n] = i;
        int l = (*n)->sequence().size();
        i += l;
    }
}

void VG::kpaths_of_node(int64_t node_id, vector<Path>& paths, int length) {
    hash_map<int64_t, Node*>::iterator n = node_by_id.find(node_id);
    if (n != node_by_id.end()) {
        Node* node = n->second;
        kpaths_of_node(node, paths, length);
    }
}

string VG::path_sequence(Path& path) {
    string sequence;
    for (int i = 0; i < path.mapping_size(); ++i) {
        sequence.append(node_by_id[path.mapping(i).node_id()]->sequence());
    }
    return sequence;
}

bool VG::is_valid(void) {

    if (node_by_id.size() != graph.node_size()) {
        cerr << "graph invalid: node count is not equal to that found in node by-id index" << endl;
        return false;
    }

    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (node_by_id.find(n->id()) == node_by_id.end()) {
            cerr << "graph invalid: node " << n->id() << " missing from by-id index" << endl;
            return false;
        }
    }

    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        int64_t f = e->from();
        int64_t t = e->to();

        //cerr << "edge " << e << " " << e->from() << "->" << e->to() << endl;

        if (node_by_id.find(f) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " (" << f << "->" << t << ") cannot find node (from) " << f << endl;
            return false;
        }
        if (node_by_id.find(t) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " (" << f << "->" << t << ") cannot find node (to) " << t << endl;
            return false;
        }

        if (edges_from_to.find(f) == edges_from_to.end()) {
            // todo check if it's in the vector
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_from_to for node " << f << endl;
            return false;
        }

        if (edges_to_from.find(t) == edges_to_from.end()) {
            // todo check if it's in the vector
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_to_from for node " << t << endl;
            return false;
        }

    }

    //cerr << "there are " << edges_from_to.size() << " items in the edges_from_to index" << endl;
    for (hash_map<int64_t, vector<int64_t> >::iterator f = edges_from_to.begin();
         f != edges_from_to.end(); ++f) {
        vector<int64_t>& to = f->second;
        for (vector<int64_t>::iterator t = to.begin(); t != to.end(); ++t) {
            Edge* e = get_edge(f->first, *t);
            //cerr << "edges_from_to " << e << " " << e->from() << "->" << e->to() << endl;
            if (!e) {
                cerr << "graph invalid, edge is null" << endl;
                return false;
            }
            if (*t != e->to() || f->first != e->from()) {
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " stored in to_from index under " << f->first << "->" << *t << endl;
                return false;
            }
            if (!has_node(e->from())) {
                cerr << "graph invalid: edge from a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
            if (!has_node(e->to())) {
                cerr << "graph invalid: edge to a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
        }
    }

    //cerr << "there are " << edges_to_from.size() << " items in the edges_to_from index" << endl;
    for (hash_map<int64_t, vector<int64_t> >::iterator t = edges_to_from.begin();
         t != edges_to_from.end(); ++t) {
        vector<int64_t>& from = t->second;
        for (vector<int64_t>::iterator f = from.begin(); f != from.end(); ++f) {
            Edge* e = get_edge(*f, t->first);
            //cerr << "edges_to_from " << e << " " << e->from() << "->" << e->to() << endl;
            if (!e) {
                cerr << "graph invalid, edge is null" << endl;
                return false;
            }
            if (t->first != e->to() || *f != e->from()) {
                cerr << "graph invalid: edge " << e->from() << "->" << e->to()
                     << " stored in to_from index under " << *f << "->" << t->first << endl;
                return false;
            }
            if (!has_node(e->from())) {
                cerr << "graph invalid: edge from a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
            if (!has_node(e->to())) {
                cerr << "graph invalid: edge to a non-existent node " << e->from() << "->" << e->to() << endl;
                return false;
            }
        }
    }

    if (head_nodes().empty()) {
        cerr << "graph invalid: no head nodes" << endl;
        return false;
    }

    if (tail_nodes().empty()) {
        cerr << "graph invalid: no tail nodes" << endl;
        return false;
    }

    //cerr << "all is well" << endl;

    return true;
}

void VG::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\"];" << endl;
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "    " << p->id() << " -> " << n->id() << ";" << endl;
    }
    out << "}" << endl;
}

void VG::to_gfa(ostream& out) {
    out << "H" << "\t" << "HVN:Z:1.0" << endl;
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        out << "S" << "\t" << n->id() << "\t" << n->sequence() << endl;
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        Edge* e = graph.mutable_edge(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "L" << "\t" << p->id() << "\t" << "-" << "\t" << n->id() << "\t" << "+" << "\t" << "0M" << endl;
    }
}

void VG::destroy_alignable_graph(void) {
    if (gssw_aligner != NULL) {
        delete gssw_aligner;
    }
}

void VG::connect_node_to_nodes(Node* node, vector<Node*>& nodes) {
    for (vector<Node*>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        create_edge(node->id(), (*n)->id());
    }
}

Node* VG::join_heads(void) {
    Node* root = create_node("N");
    vector<Node*> heads;
    head_nodes(heads);
    connect_node_to_nodes(root, heads);
    return root;
}

Alignment& VG::align(Alignment& alignment) {

    // to be completely aligned, the graph's head nodes need to be fully-connected to a common root
    Node* root = join_heads();

    gssw_aligner = new GSSWAligner(graph);
    gssw_aligner->align(alignment);
    delete gssw_aligner;
    gssw_aligner = NULL;

    // adjust the alignment position to respect the trim of the first node
    Path* path = alignment.mutable_path();
    Node* first_node = node_by_id[path->mutable_mapping(0)->node_id()];
    if (first_node->has_trim()) {
        int32_t t = path->target_position();
        path->set_target_position(t - first_node->trim());
    }

    // remove root
    destroy_node(root);
    
    return alignment;
}

Alignment VG::align(string& sequence) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment);
}

void VG::for_each_kmer_parallel(int kmer_size,
                                function<void(string&, Node*, int)> lambda,
                                int stride) {
    _for_each_kmer(kmer_size, lambda, true, stride);
}

void VG::for_each_kmer(int kmer_size,
                       function<void(string&, Node*, int)> lambda,
                       int stride) {
    _for_each_kmer(kmer_size, lambda, false, stride);
}

void VG::_for_each_kmer(int kmer_size,
                        function<void(string&, Node*, int)> lambda,
                        bool parallel,
                        int stride) {

    // use an LRU cache to clean up duplicates over the last 1mb
    // use one per thread so as to avoid contention
    map<int, string_hash_map<string, bool> > lru;
#pragma omp parallel
    {
#pragma omp single
        for (int i = 0; i < (parallel ? omp_get_num_threads() : 1); ++i) {
            lru[i];// = new string_hash_map<string, bool>();
        }
    }
    // constructs the cache key
    // experiment -- use a struct here
    auto make_cache_key = [](string& kmer, Node* node, int start) -> string {
        string cache_key = kmer;
        cache_key.resize(kmer.size() + sizeof(Node*) + sizeof(int));
        memcpy((char*)cache_key.c_str()+kmer.size(), &node, sizeof(Node*));
        memcpy((char*)cache_key.c_str()+kmer.size()+sizeof(Node*), &start, sizeof(Node*));
        return cache_key;
    };

    auto handle_path = [this,
                        lambda,
                        kmer_size,
                        stride,
                        &lru,
                        &make_cache_key](list<Node*>& path) {

        // expand the path into a vector :: 1,1,1,2,2,2,2,3,3 ... etc.
        // this makes it much easier to quickly get all the node matches of each kmer
        vector<Node*> node_by_path_position;
        expand_path(path, node_by_path_position);

        auto& cache = lru[omp_get_thread_num()];

        map<Node*, int> node_start;
        node_starts_in_path(path, node_start);
        
        // now process the kmers of this sequence
        // by first getting the sequence
        string seq = path_string(path);

        // and then stepping across the path, finding the kmers, and then implied node overlaps
        for (int i = 0; i < seq.size() - kmer_size; i+=stride) {

            // get the kmer
            string kmer = seq.substr(i, kmer_size);
            // record when we get a kmer match
            set<Node*> nodes;

            // execute our callback on each kmer/node/position
            int j = 0;
            while (j < kmer_size) {
                Node* node = node_by_path_position[i+j];
                if (nodes.find(node) == nodes.end()) {
                    nodes.insert(node);
                    int node_position = node_start[node];
                    int kmer_relative_start = i - node_position;
                    string cache_key = make_cache_key(kmer, node, kmer_relative_start);
                    if (cache.size() > 1000000) {
                        cache.clear();
                    }
                    if (cache.find(cache_key) == cache.end()) {
                        cache[cache_key] = true;
                        lambda(kmer, node, kmer_relative_start);
                    }
                }
                ++j;
            }
        }
    };

    if (parallel) {
        for_each_kpath_parallel(kmer_size, handle_path);
    } else {
        for_each_kpath(kmer_size, handle_path);
    }

}

void VG::collect_subgraph(Node* node, set<Node*>& subgraph) {

    // add node to subgraph
    subgraph.insert(node);

    // for each predecessor of node
    vector<Node*> prev;
    nodes_prev(node, prev);
    for (vector<Node*>::iterator p = prev.begin(); p != prev.end(); ++p) {
        // if it's not already been examined, collect its neighborhood
        if (!subgraph.count(*p)) {
            collect_subgraph(*p, subgraph);
        }
    }

    // for each successor of node
    vector<Node*> next;
    nodes_next(node, next);
    for (vector<Node*>::iterator n = next.begin(); n != next.end(); ++n) {
        if (!subgraph.count(*n)) {
            collect_subgraph(*n, subgraph);
        }
    }

}

void VG::disjoint_subgraphs(list<VG>& subgraphs) {
    vector<Node*> heads;
    head_nodes(heads);
    map<Node*, set<Node*> > subgraph_by_head;
    map<Node*, set<Node*>* > subgraph_membership;
    // start at the heads, but keep in mind that we need to explore fully
    for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
        if (subgraph_membership.find(*h) == subgraph_membership.end()) {
            set<Node*>& subgraph = subgraph_by_head[*h];
            collect_subgraph(*h, subgraph);
            for (set<Node*>::iterator n = subgraph.begin(); n != subgraph.end(); ++n) {
                subgraph_membership[*n] = &subgraph;
            }
        }
    }
    for (map<Node*, set<Node*> >::iterator g = subgraph_by_head.begin();
         g != subgraph_by_head.end(); ++ g) {
        set<Node*>& nodes = g->second;
        set<Edge*> edges;
        edges_of_nodes(nodes, edges);
        subgraphs.push_back(VG(nodes, edges));
    }
}

void VG::head_nodes(vector<Node*>& nodes) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (edges_to_from.find(n->id()) == edges_to_from.end()) {
            nodes.push_back(n);
        }
    }
}

vector<Node*> VG::head_nodes(void) {
    vector<Node*> heads;
    head_nodes(heads);
    return heads;
}

void VG::tail_nodes(vector<Node*>& nodes) {
    for (int i = 0; i < graph.node_size(); ++i) {
        Node* n = graph.mutable_node(i);
        if (edges_from_to.find(n->id()) == edges_from_to.end()) {
            nodes.push_back(n);
        }
    }
}

vector<Node*> VG::tail_nodes(void) {
    vector<Node*> tails;
    tail_nodes(tails);
    return tails;
}

void VG::wrap_with_null_nodes(void) {
    vector<Node*> heads;
    head_nodes(heads);
    Node* head = create_node("");
    for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
        create_edge(head, *h);
    }

    vector<Node*> tails;
    tail_nodes(tails);
    Node* tail = create_node("");
    for (vector<Node*>::iterator t = tails.begin(); t != tails.end(); ++t) {
        create_edge(*t, tail);
    }
}


    /*
Kahn's topological sort (1962)

L ← Empty list that will contain the sorted elements
S ← Set of all nodes with no incoming edges
while S is non-empty do
    remove a node n from S
    add n to tail of L
    for each node m with an edge e from n to m do
        remove edge e from the graph
        if m has no other incoming edges then
            insert m into S
if graph has edges then
    return error (graph has at least one cycle)
else 
    return L (a topologically sorted order)
    */

void VG::topological_sort(deque<Node*>& l) {
    //assert(is_valid());

    // using a map instead of a set ensures a stable sort across different systems
    map<int64_t, Node*> s;
    vector<Node*> heads;
    head_nodes(heads);
    for (vector<Node*>::iterator n = heads.begin(); n != heads.end(); ++n) {
        s[(*n)->id()] = *n;
    }

    // check that we have heads of the graph
    if (heads.empty() && graph.node_size() > 0) {
        cerr << "error:[VG::topological_sort] No heads of graph given, but graph not empty. "
             << "In-memory indexes of nodes and edges may be out of sync." << endl;
        exit(1);
    }

    while (!s.empty()) {
        Node* n = s.begin()->second;
        s.erase(n->id());
        l.push_back(n);
        vector<int64_t>& from = edges_from(n);
        vector<int64_t> to_erase;
        for (vector<int64_t>::iterator f = from.begin(); f != from.end(); ++f) {
            Node* m = node_by_id[*f];
            to_erase.push_back(m->id());
            remove_edge_tfi(n->id(), m->id());
            if (edges_to(m->id()).empty()) {
                s[m->id()] = m;
            }
        }
        for (vector<int64_t>::iterator t = to_erase.begin(); t != to_erase.end(); ++t) {
            remove_edge_fti(n->id(), *t);
        }
    }
    // if we have a cycle, signal an error, as we are not guaranteed an order
    for (hash_map<int64_t, vector<int64_t> >::iterator f = edges_from_to.begin();
         f != edges_from_to.end(); ++f) {
        if (!f->second.empty()) {
            cerr << "error:[VG::topological_sort] graph has a cycle from " << f->first
                 << " to " << f->second.front() << endl
                 << "thread " << omp_get_thread_num() << endl;
#pragma omp critical
            {
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
        }
    }
    for (hash_map<int64_t, vector<int64_t> >::iterator t = edges_to_from.begin();
         t != edges_to_from.end(); ++t) {
        if (!t->second.empty()) {
            cerr << "error:[VG::topological_sort] graph has a cycle to " << t->first
                 << " to " << t->second.front() << endl
                 << "thread " << omp_get_thread_num() << endl;
#pragma omp critical
            {
                std::ofstream out("fail.vg");
                serialize_to_ostream(out);
                out.close();
                exit(1);
            }
        }
    }
    // we have destroyed the graph's index to ensure its order
    // rebuild the indexes
    rebuild_indexes();
}

} // end namespace
