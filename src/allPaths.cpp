
/*
allPaths.cpp - This file is part of the Bayesembler (v1.2.0)


The MIT License (MIT)

Copyright (c) 2015 Lasse Maretty and Jonas Andreas Sibbesen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#include <allPaths.h>

#include <iostream>
#include <boost/regex.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/algorithm/string.hpp>
#include <tr1/unordered_set>



// Constructor
AllPaths::AllPaths (string dotfile_in) {
        
    parseGraphviz(dotfile_in);
    vertex_map = parseVertexProperties();
    junction_threshold = 1;
    
}

// Parses a graphviz formatted dot file
void AllPaths::parseGraphviz (string dot_file) {
    
	boost::dynamic_properties dp;

    // Define map between node name (dotfile) and vertex name (graph)
	name = get(boost::vertex_name, graph);
	dp.property("node_id",name); 
    
    // Define map between sequence (dotfile) and vertex color (graph)
	label = get(boost::vertex_color, graph);
	dp.property("label",label);
    
	coverage = get(boost::edge_weight, graph);
	dp.property("coverage",coverage);
	
    // Extract graph name using regular expressions "digraph <graph_name> {"
    boost::regex re("^digraph\\s{1}(\\S+)\\s{1}\\{$", boost::regex::extended);
		    
    boost::match_results<const char*> res;
    	    
    assert (boost::regex_search(dot_file.c_str(), res, re));
    graph_name = res[1];

    // Parse dot file into graph object
    bool status = read_graphviz(dot_file, graph, dp, "node_id");
        
    assert (status);
}



// Parses vertex properties (labels in dot file)
AllPaths::VertexMap AllPaths::parseVertexProperties () {

    VertexMap vertex_map_temp;
                
    pair<VertexIterator, VertexIterator> vp;
    
    for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {
        
        vector <string> vertex_prop_vec;
        boost::split(vertex_prop_vec, label[*vp.first], boost::is_any_of(";"));
    
        Vertex vertex_temp;
        vertex_temp.reference = vertex_prop_vec[0];
        vertex_temp.strand = vertex_prop_vec[1];
        vertex_temp.start = atoi(vertex_prop_vec[2].c_str());
        vertex_temp.end = atoi(vertex_prop_vec[3].c_str());
        
        vertex_map_temp[name[*vp.first]] = vertex_temp;     
    }

    VertexMap::iterator first_it = vertex_map_temp.begin();

    for (VertexMap::iterator it = vertex_map_temp.begin(); it != vertex_map_temp.end(); it++) {

        assert (first_it->second.reference == it->second.reference);
        assert (first_it->second.strand == it->second.strand);
        assert (first_it->second.end >= first_it->second.start);
    }
    
    return vertex_map_temp;
}        

// Finds all full paths in a directed acyclic graph
AllPaths::CandidateVector AllPaths::findFullPaths (bool no_pre_mrna, int max_candidate_number, stringstream& log_stream) {
	
    CandidateVector all_candidates;
    
    tr1::unordered_set <string> non_sinks;
    tr1::unordered_set <string> non_source;
    
    EdgeIterator ei, ei_end;
	for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
	
        non_source.insert(name[target(*ei, graph)]);
        non_sinks.insert(name[source(*ei, graph)]);
	}
 
    boost::graph_traits < SpliceGraph >::vertices_size_type vertex_num_graph = num_vertices(graph);

    log_stream << "\n[" << getLocalTime() << "] " << "Number of vertices in graph: " << vertex_num_graph << endl;   
    log_stream << "[" << getLocalTime() << "] " << "Number of sources in graph: " << vertex_num_graph - non_source.size() << endl;
    log_stream << "[" << getLocalTime() << "] " << "Number of sinks in graph: " << vertex_num_graph - non_sinks.size() << endl;      
    
	if ((non_source.size() == 0) and (non_sinks.size() == 0) and (!no_pre_mrna) and (vertex_num_graph == 2)) {
		
		if ((vertex_map["0"].start == vertex_map["pre_mrna"].start) and (vertex_map["0"].end == vertex_map["pre_mrna"].end)) {
			
			assert (vertex_map["0"].reference == vertex_map["pre_mrna"].reference);
			assert (vertex_map["0"].strand == vertex_map["pre_mrna"].strand);
			assert (num_edges(graph) == 0);
			
		    pair<VertexIterator, VertexIterator> vp;
            
            VertexDescriptor pre_mrna_vertex;
            bool found_pre_mRNA = false;
            			
		    log_stream << "[" << getLocalTime() << "] " << "WARNING: Removed pre-mRNA from candidates due to exact match with another candidate ..." << endl;
			
			for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {  
				
				if (name[*vp.first] == "pre_mrna") {
					
                    assert (!found_pre_mRNA);
                    pre_mrna_vertex = *vp.first;
                    found_pre_mRNA = true;
					break;
				}	
			}

            assert (found_pre_mRNA);
            
            remove_vertex(pre_mrna_vertex, graph);
		}
	}
    
    bool too_many_candidates = true;
    int candidate_idx;
    int path_count;
    
    while (too_many_candidates) {
        
        candidate_idx = 0;
        path_count = 0;
        all_candidates.clear();
        non_source.clear();
        
        vector<EdgeDescriptor> edges_to_remove;
        
        EdgeIterator ei, ei_end;
    	for (tie(ei, ei_end) = edges(graph); ei != ei_end; ei++) {
		
            if (junction_threshold > coverage[*ei]) {
            
                edges_to_remove.push_back(*ei);
            }
    	} 
        
        for (int i = 0; i < edges_to_remove.size(); i++) {
            
            remove_edge(edges_to_remove[i], graph);
        }
                                
    	for (tie(ei, ei_end) = edges(graph); ei != ei_end; ei++) {
		
            non_source.insert(name[target(*ei, graph)]);
    	}

        pair<VertexIterator, VertexIterator> vp;
    	for (vp = vertices(graph); vp.first != vp.second; vp.first++) {        
                 
            if (non_source.count(name[*vp.first])) {

                continue;
            }
		
            // Initilize vector of paths
    		Vertices path;
    		VerticesVector paths;
                
            // Find all paths between two vertices
    		too_many_candidates = findAllPathsSource(*vp.first, path, paths, path_count, max_candidate_number, log_stream);

            if (too_many_candidates) {
                
                junction_threshold++;
                assert((path_count - 1) == max_candidate_number);
                break;
            }

    		for (int i = 0; i < paths.size(); i++) {	
            		
    			int candidate_len = 0;
                VertexVector vertex_vector_temp;
                Candidate candidate_temp;
    
                stringstream seq_name;
                seq_name << graph_name;
			
    			// Print paths
    			for (int j = 0; j < paths[i].size(); j++) {
				
                    candidate_len += vertex_map[name[paths[i][j]]].end - vertex_map[name[paths[i][j]]].start + 1;
                                
                    if ((j > 0) and ((vertex_map[name[paths[i][j-1]]].end + 1) == vertex_map[name[paths[i][j]]].start)) {
                    
                        vertex_vector_temp.back().end = vertex_map[name[paths[i][j]]].end;
                                  
                    } else {
                    
                        vertex_vector_temp.push_back(vertex_map[name[paths[i][j]]]);
                    }              
    			}
						            
                if (name[paths[i][0]] == "pre_mrna") {
                
					assert(paths[i].size() == 1);
                    assert(vertex_vector_temp.size() == 1);
                
                    if (!no_pre_mrna) {
                                      
                        seq_name << "_pre";
                    
                        candidate_temp.idx = candidate_idx;
                        candidate_temp.name = seq_name.str();
                        candidate_temp.graph_name = graph_name;
                        candidate_temp.length = candidate_len;
                        candidate_temp.vertices = vertex_vector_temp;
                        candidate_temp.isPreMrna = true;
                        candidate_idx += 1;                    
                        all_candidates.push_back(candidate_temp);
                    }
                                    
                } else {
			
                    // Create sequence name
                    seq_name << "_seq" << candidate_idx;
                
                    candidate_temp.idx = candidate_idx;
                    candidate_temp.name = seq_name.str();
                    candidate_temp.graph_name = graph_name;
                    candidate_temp.length = candidate_len;
                    candidate_temp.vertices = vertex_vector_temp;
                    candidate_temp.isPreMrna = false;
                    candidate_idx += 1;
                    all_candidates.push_back(candidate_temp);
    			}
    		}
    	}
        
        if (num_edges(graph) == 0) {
            
            too_many_candidates = false;
        }
    }
	
	log_stream << "[" << getLocalTime() << "] " << "Number of paths in splice-graph (transcripts): " << candidate_idx << endl;

    return all_candidates;
}

// Finds all possible paths between two vertices
bool AllPaths::findAllPathsSource (VertexDescriptor s, Vertices path, VerticesVector& paths, int& path_count, int max_candidate_number, stringstream& log_stream) {
	
    // Initialise vector of paths
	VerticesVector temp;
    
    // Extend path
	path.push_back(s);
	
    // Check if source vertex has reached target
    AdjacencyIterator ai_init, ai_end;
    tie(ai_init, ai_end) = adjacent_vertices(s, graph);
    
	if (ai_init == ai_end) {
        
	 	paths = VerticesVector(1, path);
        path_count++;
	
        if (path_count > max_candidate_number) {
            
            log_stream << "[" << getLocalTime() << "] Too many candidates to continue - adjusting junction threshold from " << junction_threshold << " to " << junction_threshold + 1 << endl;
            return true;
        }
        
	} else {
	
		paths = temp;
        
        // Loop over all vertices adjacent to source
        AdjacencyIterator ai_init, ai_end;
		for (tie(ai_init, ai_end) = adjacent_vertices(s, graph); ai_init != ai_end; ++ai_init) {
            
            // Find all paths from adjacent vertex to target
            if (findAllPathsSource(*ai_init, path, paths, path_count, max_candidate_number, log_stream)) {
                
                return true;
            }
            			
            // Add new paths to paths array
            for (int i = 0; i < paths.size(); i++) {

                temp.push_back(paths[i]);			
            }

            paths = temp;			    
		}
	}
    
    return false;
}











