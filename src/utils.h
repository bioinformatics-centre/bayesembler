
/*
utils.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef UTILS_H
#define UTILS_H

#include <tr1/unordered_map>
#include <list>
#include <string>
#include <iomanip>
#include <vector>
#include <time.h>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>
#include <math.h>
#include <limits>

typedef unsigned int uint;

using namespace std;

const double double_underflow = numeric_limits<double>::min();
const double double_overflow = numeric_limits<double>::max();
const double double_precision = numeric_limits<double>::epsilon();
const double double_almost_one = 1-double_precision;
const double accumulated_precision_error = 100;

inline bool double_compare(double a, double b) {

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision * accumulated_precision_error));
}

typedef pair<uint, vector<BamTools::CigarOp> > ReadData;


struct FragmentLengthContainer {

    uint transcript_length; 
    uint map_count;
    tr1::unordered_map<uint,uint> fragment_lengths;
};

struct FragmentAlignment {

    pair<ReadData, ReadData> read_data;
    double probability;
};


struct GraphInfo {

    list<string> dot_strings;
    string reference;
    string strand;
    int left;
    int right;
    bool single_path_graph;
    uint read_count;
};

struct Ensemble {

    vector<int> indices;
    Eigen::VectorXd means;
    Eigen::VectorXd sd;

};

// Container for vertex properties
struct Vertex {
    
    string reference;
    string strand;
    int start;
    int end;
};

// Container for candidate properties
struct Candidate {
    
    int idx;
    string name;
    string graph_name;
    int length;
    vector<Vertex> vertices;
    bool isPreMrna;
};

struct AssemblyInfo {
    
    string name;
    string ensemble;
    string candidates;
    int return_code;
    string assembly_log;
    
};

struct CollapsedMap {

    Eigen::MatrixXd probability_matrix;
    vector<int> counts;
    tr1::unordered_map<int,pair<int,double> > new_idx_to_old_idx_effective_length;
    int return_code;
   
};


// Container for option variables (all parameters used in the code)
struct OptionsContainer {

    string bam_file;

    string samtools_path;
    string cem_processsam_path;

    int num_threads;
    string strand_specific;
    double marginal_threshold;
    double count_threshold;
    bool no_pre_mrna;

    int seed;
    string output_prefix;
    string output_mode;
    int total_num_fragments;
    bool keep_temp_files;

    int max_candidate_number; 
    double exon_base_coverage;

    double gamma;
    double frag_mean; 
    double frag_sd;
    bool estimate_fragment_length; 
    uint fragment_est_min_transcript_length;

    unsigned int gibbs_base_iterations;
    unsigned int gibbs_scale_iterations;
};
	
// Vector flush capability
template < typename T >
ostream& operator << (ostream& os, const vector<T>& v) {
    
    for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {
        
		os << scientific << setprecision(6) << *it << "  ";	
	}

    return os;
}

// 2D-vector flush capability
template < typename T >
ostream& operator << (ostream& os, const vector<vector<T> >& v) {
    
    for (typename vector<vector<T> >::const_iterator it = v.begin(); it != v.end(); it++) {
		
        os << *it << "\n";
	}
	
    return os;
}


inline string getLocalTime () {
    
    time_t now = time(NULL);
    struct tm * lt;
    lt = localtime (&now);
    
    stringstream output_string;
        
    if (lt->tm_mday < 10) {
        
        output_string << "0" << lt->tm_mday << "/";
        
    } else {
        
        output_string << lt->tm_mday << "/";
        
    }
    
    if ((1 + lt->tm_mon) < 10) {
    
        output_string << "0" << (1 + lt->tm_mon) << "/" << (1900 + lt->tm_year) << " ";
        
    } else {
            
        output_string << (1 + lt->tm_mon) << "/" << (1900 + lt->tm_year) << " ";
        
    }

    if (lt->tm_hour < 10) {
        
        output_string << "0" << lt->tm_hour << ":";
        
    } else {
        
        output_string << lt->tm_hour << ":";
        
    }
    
    if (lt->tm_min < 10) {
        
        output_string << "0" << lt->tm_min << ":";
        
    } else {
        
        output_string << lt->tm_min << ":";
        
    }

    if (lt->tm_sec < 10) {
        
        output_string << "0" << lt->tm_sec;
        
    } else {
        
        output_string << lt->tm_sec;
        
    }
    
    return output_string.str();
    
}
	
#endif