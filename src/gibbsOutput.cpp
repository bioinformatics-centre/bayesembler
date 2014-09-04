
/*
gibbsOutput.cpp - This file is part of the Bayesembler (v1.1.1)


The MIT License (MIT)

Copyright (c) 2014 Lasse Maretty and Jonas Andreas Sibbesen

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


#include <gibbsOutput.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <sstream>
#include <bitset>
#include <math.h>

using namespace std; 

GibbsOutput::GibbsOutput(int num_iterations_in, int num_transcripts_in) {
    
    num_iterations = num_iterations_in;
    num_transcripts = num_transcripts_in;
    
    marginal_occurence = Eigen::VectorXi::Zero(num_transcripts);
    marginal_mean_cumulant = Eigen::VectorXd::Zero(num_transcripts);
    marginal_var_cumulant = Eigen::VectorXd::Zero(num_transcripts);
    
}

void GibbsOutput::addSample(vector<int> plus, vector<double> values) {
    
    Eigen::VectorXd plus_values(plus.size());
     
	sort(plus.begin(), plus.end());
    
    for (int i=0; i < plus.size(); i++) {
        
        int idx = plus[i];        
        marginal_occurence(idx)++;
        
        // If first occurence
        if (marginal_occurence(idx) == 1) {
            
            marginal_mean_cumulant(idx) = values[idx];
                  
        } else {
            
            double delta = values[idx] - marginal_mean_cumulant(idx);
            marginal_mean_cumulant(idx) += delta/marginal_occurence(idx);
            marginal_var_cumulant(idx) += delta * (values[idx] - marginal_mean_cumulant(idx));       
        }
    }
    
}

Ensemble GibbsOutput::getTopMarginals(double marginal_threshold) {

    double count_threshold = marginal_threshold * num_iterations;
    vector<int> indices;

    for (int i=0; i < num_transcripts; i++) {
        
        if (marginal_occurence(i) >= count_threshold) {
            
            indices.push_back(i);
        }
    }
    
    Eigen::VectorXd means(indices.size());
    Eigen::VectorXd sds(indices.size());
    
    for (int i=0; i < indices.size(); i++) {
        
        means(i) = marginal_mean_cumulant(indices[i]);
        sds(i) = sqrt(marginal_var_cumulant(indices[i])/(marginal_occurence(indices[i])-1));
    }

    Ensemble top_marginals;
    
    top_marginals.indices = indices;
    top_marginals.means = means;
    top_marginals.sd = sds;
    
    return top_marginals;
}


string GibbsOutput::writeEnsembleGtf(Ensemble top_ensemble, vector<Candidate> candidates, tr1::unordered_map<int,pair<int,double> > new_idx_to_old_idx_effective_length, int total_num_fragments, double count_threshold) {

    stringstream output_string;

    for (int i=0; i < top_ensemble.indices.size(); i++) {
    
        int new_idx = top_ensemble.indices[i];
        int old_idx = new_idx_to_old_idx_effective_length[new_idx].first;
        double expected_count = top_ensemble.means(i) * new_idx_to_old_idx_effective_length[new_idx].second * total_num_fragments/1e9;

        if (expected_count < count_threshold) {

            continue;

        }

        Candidate current_candidate = candidates[old_idx];
        
        assert(current_candidate.idx == old_idx);
                
        if (!current_candidate.isPreMrna) {
        
            vector<Vertex> current_vertices = current_candidate.vertices;
        
            for (int j=0; j < current_vertices.size(); j++) {
                
                output_string << current_vertices[j].reference << "\tBayesembler\texon\t" << current_vertices[j].start << "\t" << current_vertices[j].end << "\t.\t" << current_vertices[j].strand << "\t.\tgene_id \"" << current_candidate.graph_name << "\"; transcript_id \"" << current_candidate.name << "\"; transcript_confidence \"" << ((double) marginal_occurence(new_idx)/num_iterations) << "\"; FPKM \"" << top_ensemble.means(i) << "\"; FPKM_sd \"" << top_ensemble.sd(i) << "\"; expected_count \"" << expected_count << "\";" << endl;            
            }
        }
    }

    return output_string.str();
}

string GibbsOutput::writeCandidateGtf(vector<Candidate> candidates, tr1::unordered_map<int,pair<int,double> > new_idx_to_old_idx_effective_length, int total_num_fragments) {
        
    stringstream output_string;

    for (int i=0; i < num_transcripts; i++) {
    
        int old_idx = new_idx_to_old_idx_effective_length[i].first;
        double expected_count = marginal_mean_cumulant(i) * new_idx_to_old_idx_effective_length[i].second * total_num_fragments/1e9;

        Candidate current_candidate = candidates[old_idx];
        
        assert(current_candidate.idx == old_idx);
                
        vector<Vertex> current_vertices = current_candidate.vertices;
    
        for (int j=0; j < current_vertices.size(); j++) {

            output_string << current_vertices[j].reference << "\tBayesembler\texon\t" << current_vertices[j].start << "\t" << current_vertices[j].end << "\t.\t" << current_vertices[j].strand << "\t.\tgene_id \"" << current_candidate.graph_name << "\"; transcript_id \"" << current_candidate.name << "\"; transcript_confidence \"" << ((double) marginal_occurence(i)/num_iterations) << "\"; FPKM \"" << marginal_mean_cumulant(i) << "\"; FPKM_sd \"" << sqrt(marginal_var_cumulant(i)/(marginal_occurence(i) - 1)) << "\"; expected_count \"" << expected_count << "\";" << endl;                         
        }
    }

    return output_string.str();
    
}



