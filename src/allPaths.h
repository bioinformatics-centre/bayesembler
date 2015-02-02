
/*
allPaths.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef ALLPATHS_H
#define ALLPATHS_H


#include <string>
#include <tr1/unordered_map>
#include <utils.h>
#include <boost/graph/adjacency_list.hpp>
#include <vector>

using namespace std;



class AllPaths {

  private:
        
    // Define vertex and edge properties
	typedef boost::property < boost::vertex_name_t, string, boost::property < boost::vertex_color_t, string > > VertexProperty;
    typedef boost::property < boost::edge_weight_t, int > EdgeProperty;

	// Define graph structure
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS, VertexProperty, EdgeProperty > SpliceGraph;
	
	typedef boost::graph_traits<SpliceGraph>::vertex_descriptor VertexDescriptor;
    typedef boost::graph_traits<SpliceGraph>::edge_descriptor EdgeDescriptor;
    typedef vector <VertexDescriptor> Vertices;
    typedef vector <vector <VertexDescriptor> > VerticesVector;
	
    typedef boost::graph_traits<SpliceGraph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<SpliceGraph>::edge_iterator EdgeIterator;
    typedef boost::graph_traits<SpliceGraph>::adjacency_iterator AdjacencyIterator;

    // Define map between graph properties and its associated value
	boost::property_map <SpliceGraph, boost::vertex_name_t>::type name;
	boost::property_map <SpliceGraph, boost::vertex_color_t>::type label;
    boost::property_map <SpliceGraph, boost::edge_weight_t>::type coverage;

    string graph_name;
    
    // Initialise graph
	SpliceGraph graph;
    
    int junction_threshold;
    
	typedef vector <Vertex> VertexVector;
    typedef vector <Candidate> CandidateVector;
    typedef tr1::unordered_map <string, Vertex> VertexMap;
    
    // Graph transversal utility function
    bool findAllPathsSource(VertexDescriptor, Vertices, VerticesVector&, int&, int, stringstream&);
    
    VertexMap vertex_map;
        
    void parseGraphviz (string);
    VertexMap parseVertexProperties();

    
  public:


    // Constructor
    AllPaths(string);
    
    // Graph parsing and transversal runner functions
    CandidateVector findFullPaths(bool, int, stringstream&);
    		
};


#endif
