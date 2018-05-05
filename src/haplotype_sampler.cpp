#include "haplotype_sampler.hpp"
#include "handle.hpp"
#include "snarls.hpp"
#include "multipath_alignment.hpp"
#include "vg.hpp"

using namespace std;
using namespace vg;

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

DirectedGraph build_graph(const Snarl&, const HandleGraph*);
DirectedGraph reverse_graph(DirectedGraph&);
bool is_head(DirectedGraph, int);
vector<int> topological_sort_graph(DirectedGraph);
vector<int> get_nodes(DirectedGraph&);
void print_graph(DirectedGraph&);
void print_count_graph(CountGraph&);
CountGraph count_forward_paths(DirectedGraph&);
CountGraph count_reverse_paths(DirectedGraph&);
CountGraph count_paths(DirectedGraph&);
vector<double> get_distribution(CountGraph&, DirectedGraph&, int);
void multipath_node_ids(const MultipathAlignment&, set<int>&);
map<MultipathAlignment*, set<int>> construct_node_map(vector<MultipathAlignment>&);
map<int, set<Snarl*>> construct_node_snarl_map(SnarlManager*, VG*);
map<Snarl*, set<MultipathAlignment*>> construct_snarl_alignment_map(vector<MultipathAlignment>&, map<int, set<Snarl*>>&, map<MultipathAlignment*, set<int>>&);


list<Visit> sample_path(Snarl& snarl, HandleGraph* graph){
	DirectedGraph digraph = build_graph(snarl, graph);
	CountGraph cgraph = count_paths(digraph);

	int current = digraph.head;
	list<Visit> path;
	path.push_back(graph->to_visit(graph->get_handle(digraph.head)));

	while(current != digraph.tail){
		vector<double> dist = get_distribution(cgraph, digraph, current);

		double p = ((double) rand() / (RAND_MAX));
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

list<Visit> random_sequence(SnarlManager* snarl_manager, HandleGraph* graph){
	list<Visit> sequence;
	vector<const Snarl*> top_snarls = snarl_manager -> children_of(nullptr);

	for(int i = 0; i < top_snarls.size(); i++){
		Snarl snarl = *top_snarls[i];

		list<Visit> subsequence = sample_path(snarl, graph);

		DirectedGraph digraph = build_graph(snarl, graph);

		if(i == 0){
			sequence.insert(sequence.end(), subsequence.begin(), subsequence.end());
		}else{
			if(sequence.back() == subsequence.front())
				sequence.insert(sequence.end(), ++subsequence.begin(), subsequence.end());
			else
				sequence.insert(sequence.end(), subsequence.begin(), subsequence.end());
		}
	}

	return sequence;
}

// [NAIVE]
map<int, list<Visit>::iterator> visit_map(list<Visit> sequence){
	map<int, list<Visit>::iterator> vmap;
	
	for(list<Visit>::iterator it = sequence.begin(); it != sequence.end(); it++){
		vmap[(*it).node_id()] = it;
	}

	return vmap;
}

// [NAIVE] Assumes subsequence refers to snarl represented in genome
// subsequence should contain both endpoints of snarl
// construct a new visit map after calling
void replace_subsequence(list<Visit> &sequence, list<Visit> subsequence){
	Visit start = subsequence.front();
	Visit end = subsequence.back();

	typedef list<Visit>::iterator iter;
	iter first = sequence.begin();

	// Remove the original
	while(first->node_id() != start.node_id())
		first++;

	iter last = first;
	while(last->node_id() != end.node_id())
		last++;
	last++;

	sequence.erase(first, last);

	// Insert the new
	iter position = last;
	sequence.insert(last, subsequence.begin(), subsequence.end());
}

// [NAIVE]
bool in_sequence(Snarl snarl, list<Visit> sequence){
	list<Visit>::iterator it = sequence.begin();

	while(it != sequence.end() || it->node_id() != snarl.start().node_id())
		it++;

	return it != sequence.end();
}

// [NAIVE]
// Selects a random top level snarl
Snarl random_snarl(SnarlManager* manager){
	vector<const Snarl*> snarls = manager -> children_of(nullptr);
	int index = rand() % snarls.size();
	return *snarls[index];
}

map<Snarl*, set<MultipathAlignment*>> construct_snarl_map(SnarlManager* snarl_manager, vector<MultipathAlignment>& mp, VG* graph){
	map<MultipathAlignment*, set<int>> mp_map = construct_node_map(mp);
	map<int, set<Snarl*>> id_map = construct_node_snarl_map(snarl_manager, graph);
	return construct_snarl_alignment_map(mp, id_map, mp_map);
}

// ********************************************************************************
// 							Private methods (for internal use)
// ********************************************************************************

// Constructs a directed graph representation of the given snarl
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

// Constructs the reversed directed graph
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

// Determines if a node is a head node (no incoming edges)
bool is_head(DirectedGraph graph, int node){
	DirectedGraph rgraph = reverse_graph(graph);
	return rgraph.edges[node].size() == 0;
}

// Performs topological sort on a digraph using Kahn's algorithm
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

// Returns a list of all nodes in the graph
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

// Debug method to print out directed graph
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

// Debug method to print out count graph
void print_count_graph(CountGraph &graph){
	for(map<int, int>::iterator iter = graph.count.begin(); iter != graph.count.end(); iter++){
		int node = iter->first;
		int count = iter->second;
		cout << node << ": " << count << endl;
	}
	cout << endl;
}

// Builds count graph where value at node indicates
// number of paths from the head to the node.
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

// Builds count graph where value at node indicates
// number of paths from the node to the tail
CountGraph count_reverse_paths(DirectedGraph &graph){
	DirectedGraph rgraph = reverse_graph(graph);
	CountGraph cgraph = count_forward_paths(rgraph);
	return cgraph;
}

// Builds count graph where value at node indicates
// number of paths from head to tail that include
// that node
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

vector<double> get_distribution(CountGraph& cgraph, DirectedGraph& digraph, int node){
	double total = 0.0;
	for(int neighbor: digraph.edges[node])
		total += cgraph.count[neighbor];

	vector<double> dist;
	for(int neighbor: digraph.edges[node])
		dist.push_back(cgraph.count[neighbor]/total);

	return dist;
}

void multipath_node_ids(const MultipathAlignment &mp, set<int> &nodes){
	for(int i = 0; i < mp.subpath().size(); i++){
		for(int j = 0; j < mp.subpath()[i].path().mapping().size(); j++){
			int node = mp.subpath()[i].path().mapping()[j].position().node_id();
			nodes.insert(node);
		}
	}
}

map<MultipathAlignment*, set<int>> construct_node_map(vector<MultipathAlignment>& mp){
	map<MultipathAlignment*, set<int>> mp_map;
    for(int i = 0; i < mp.size(); i++){
    	set<int> nodes;
    	multipath_node_ids(mp[i], nodes);
    	mp_map[&mp[i]] = nodes;
    }
    return mp_map;
}

map<int, set<Snarl*>> construct_node_snarl_map(SnarlManager* snarl_manager, VG* graph){
	map<int, set<Snarl*>> id_map;

	snarl_manager -> for_each_snarl_preorder([&](const Snarl* s){
		pair<unordered_set<Node*>, unordered_set<Edge*>> p = snarl_manager -> deep_contents(s, *graph, true);
		unordered_set<Node*> nodes = p.first;
		for(unordered_set<Node*>::const_iterator it = nodes.begin(); it != nodes.end(); it++){
			int node_id = (*it) -> id();
			id_map[node_id].insert((Snarl*) s);
		}
	});

	return id_map;
}

map<Snarl*, set<MultipathAlignment*>> construct_snarl_alignment_map(vector<MultipathAlignment>& mp, map<int, set<Snarl*>> &id_map, map<MultipathAlignment*, set<int>> &mp_map){
	map<Snarl*, set<MultipathAlignment*> > snarl_map;
	for(int i = 0; i < mp.size(); i++){
    	set<int> nodes = mp_map[&mp[i]];
    	for(set<int>::const_iterator node = nodes.begin(); node != nodes.end(); node++){
    		for(set<Snarl*>::const_iterator snarl = id_map[*node].begin(); snarl != id_map[*node].end(); snarl++){
    			snarl_map[*snarl].insert(&mp[i]);
    		} 
    	}
    }
	return snarl_map;
}