
/*
assignmentGenerator.cpp - This file is part of the Bayesembler (v1.2.0)


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


#include <assignmentGenerator.h>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <algorithm>

typedef boost::random::uniform_int_distribution<> uniform_sampler_t;
typedef boost::random::binomial_distribution<> binomial_distribution_t;
typedef boost::random::variate_generator<boost::random::mt19937*, boost::random::binomial_distribution<> > binomial_sampler_t;
typedef boost::random::uniform_01<boost::random::mt19937*> uniform_01_sampler_t;


AssignmentGenerator::AssignmentGenerator(mt_rng_pt_t mt_rng_pt_in, Eigen::MatrixXd * collapsed_map_eigen_in, vector<int> * collapsed_count_vec_in) {

	mt_rng_pt = mt_rng_pt_in;
    collapsed_map_eigen = collapsed_map_eigen_in;
    collapsed_count_vec = collapsed_count_vec_in;
    num_transcripts = collapsed_map_eigen->cols();

}

void AssignmentGenerator::calcMinimumReadCover(stringstream& log_stream) {

    log_stream << "\n[" << getLocalTime() << "] " << "Calculating minimum read cover" << endl;
    
    tr1::unordered_set<int> minimum_set_temp;
    
    Eigen::Matrix< int,1,Eigen::Dynamic > logical_vector(collapsed_map_eigen->rows());
    logical_vector.setOnes();
    
    Eigen::Matrix< int,Eigen::Dynamic,Eigen::Dynamic > binary_collapsed_map_eigen((collapsed_map_eigen->array() >= double_underflow).matrix().cast<int>());
    
    while (logical_vector.sum() > 0) {
        
        Eigen::RowVectorXi covered_reads(logical_vector * binary_collapsed_map_eigen);
        int max_val = covered_reads.maxCoeff();
        assert (max_val > 0);
        
        vector<int> max_indices;
        max_indices.reserve(covered_reads.size());
        
        for (int i=0; i < covered_reads.size(); i++) {
            
            if (max_val == covered_reads[i]) {
                
                max_indices.push_back(i);
                
            }
        }
        
        uniform_sampler_t set_index_dist(0,max_indices.size()-1);
        int max_pos = max_indices[set_index_dist(*mt_rng_pt)];
        
        minimum_set_temp.insert(max_pos);
        logical_vector = (logical_vector.array() - logical_vector.array() * binary_collapsed_map_eigen.col(max_pos).transpose().array()).matrix();
    }
    
    minimum_set = minimum_set_temp;
    
    log_stream << "[" << getLocalTime() << "] " << "Number of candidates in minimum read cover: " << minimum_set.size() << "\n" << endl;
}

int AssignmentGenerator::getMinimumReadCoverSize() {
    
    return minimum_set.size();
}

CountValueContainer AssignmentGenerator::initAssignmentMinimum(stringstream& log_stream) {

    log_stream << "[" << getLocalTime() << "] " << "Initialising read assignments using minimum set cover" << endl;
        
    CountValueContainer map_counts(num_transcripts);
        
    for (int i=0; i < collapsed_count_vec->size(); i++) {
        
        vector<int> minimum_set_vec;
        double norm_const = 0;
		
		for (int j=0; j < collapsed_map_eigen->cols(); j++) {
            
			if ((minimum_set.count(j) == 1) and ((*collapsed_map_eigen)(i,j) >= double_underflow)) {
                
				minimum_set_vec.push_back(j);
                norm_const += (*collapsed_map_eigen)(i,j);
                
			}
		}

        assert(norm_const >= double_underflow);
        int num_samples = collapsed_count_vec->at(i);
        double divider = 1;
        
        for (int j = 0; j < minimum_set_vec.size(); j++) {
            
            double norm_probability = (*collapsed_map_eigen)(i,minimum_set_vec[j])/norm_const;
            
            if (norm_probability >= double_underflow) {
                
                binomial_distribution_t binom_dist(num_samples, norm_probability/divider);
                binomial_sampler_t binom_sampler(mt_rng_pt, binom_dist);
                
                int num_maps = binom_sampler();
                map_counts.addToCount(num_maps,minimum_set_vec[j]);
                
                if (num_maps > 0) {
                    
                    map_counts.addToPlus(minimum_set_vec[j]);
                    
                }
                
                num_samples -= num_maps;                
            }
            
            divider -= norm_probability;
        }
        
        assert (num_samples == 0);      
	}
    
	map_counts.fetchNullSet();
    vector <int> counts_out = map_counts.getCounts();
	// cout << "Initialised assignment with counts:\n" << counts_out << endl;
	
	return map_counts;
    
}

CountValueContainer AssignmentGenerator::initEnsembleAssignment(vector<int> indices) {
        
    CountValueContainer map_counts(num_transcripts);
                
    for (int i=0; i < collapsed_count_vec->size(); i++) {
                
        double norm_const = 0;
        
		for (int j=0; j < indices.size(); j++) {
            
            norm_const += (*collapsed_map_eigen)(i,indices[j]);
    	}
        
        assert(norm_const >= double_underflow);
        int num_samples = collapsed_count_vec->at(i);
        double divider = 1;
        
        for (int j = 0; j < indices.size(); j++) {
            
            int current_idx = indices[j];
            double norm_probability = (*collapsed_map_eigen)(i, current_idx)/norm_const;
            
            if (norm_probability >= double_underflow) {
                
                binomial_distribution_t binom_dist(num_samples, norm_probability/divider);
                binomial_sampler_t binom_sampler(mt_rng_pt, binom_dist);
                
                int num_maps = binom_sampler();
                map_counts.addToCount(num_maps, current_idx);
                
                if (num_maps > 0) {
                    
                    map_counts.addToPlus(current_idx);
                    
                }
                
                num_samples -= num_maps;
                
            }
            
            divider -= norm_probability;
            
        }
        
        assert (num_samples == 0);    
	}
    
	map_counts.fetchNullSet();
    vector <int> counts_out = map_counts.getCounts();
	// cout << "Initialised assignment with counts:\n" << counts_out << endl;
	
	return map_counts;
    
}

CountValueContainer AssignmentGenerator::initAssignment(stringstream& log_stream) {

    log_stream << "[" << getLocalTime() << "] " << "Initialising read assignments uniformly" << endl;
        
    CountValueContainer map_counts(num_transcripts);
        
    for (int i=0; i < collapsed_count_vec->size(); i++) {
        
        int num_samples = collapsed_count_vec->at(i);
        double divider = 1;
        
        for (int j = 0; j < num_transcripts; j++) {
            
            double norm_probability = (*collapsed_map_eigen)(i,j);
            
            if (norm_probability >= double_underflow) {
                
                binomial_distribution_t binom_dist(num_samples, norm_probability/divider);
                binomial_sampler_t binom_sampler(mt_rng_pt, binom_dist);
                
                int num_maps = binom_sampler();
                map_counts.addToCount(num_maps,j);
                
                if (num_maps > 0) {
                    
                    map_counts.addToPlus(j);
                    
                }
                
                num_samples -= num_maps;
                
            }
            
            divider -= norm_probability;
            
        }
        
        assert (num_samples == 0);
        
	}
    
	map_counts.fetchNullSet();
    vector <int> counts_out = map_counts.getCounts();
	// cout << "Initialised assignment with counts:\n" << counts_out << endl;
	
	return map_counts;
    
}


CountValueContainer AssignmentGenerator::generateAssignment(ExpressionValueContainer expression) {
	
    // Reset assignment map
	CountValueContainer map_counts(num_transcripts);
    
    // Create expression vector and calculate normalisation constants
	Eigen::RowVectorXd expression_values(num_transcripts);
	
	for (int i=0; i < num_transcripts; i++) {
	    
		expression_values(i) = expression.getValue(i);
	}
      
    uniform_01_sampler_t sample_uniform_01(mt_rng_pt);
    
    for (int i=0; i < collapsed_count_vec->size(); i++) {
        
        Eigen::RowVectorXd probabilities_unnormalised = collapsed_map_eigen->row(i).array() * expression_values.array();
                
        double normalisation_constant = probabilities_unnormalised.sum();
        
        assert (normalisation_constant >= double_underflow);
        assert (probabilities_unnormalised.size() == num_transcripts);
        
        int num_samples = collapsed_count_vec->at(i);
        
        if (num_samples < 20) {
        
            vector <double> uniform_samples(num_samples);
            
            for (int j = 0; j < num_samples; j++) {
                
                uniform_samples[j] = sample_uniform_01();            
            }
            
            sort (uniform_samples.begin(), uniform_samples.end());
            double cum_sum = 0;
            int index_num = 0;
            
            for (int j = 0; j < num_transcripts; j++) {
            
                cum_sum += probabilities_unnormalised[j]/normalisation_constant;
                
                if (uniform_samples[index_num] < cum_sum) {
                    
                    int samp_counter = 1;
                    index_num +=1;
                                        
                    while (index_num < uniform_samples.size()) {
                        
                        if (uniform_samples[index_num] >= cum_sum) {
                        
                            break;   
                        }
                                                
                        samp_counter += 1;
                        index_num +=1;
                        
                    }
                    
                    map_counts.addToCount(samp_counter, j);
                    map_counts.addToPlus(j);
                    
                }
                
                if (index_num == uniform_samples.size()) {
                    
                    break;
                    
                }
              
            }
                    
        } else {
        
            double divider = 1;
        
            for (int j = 0; j < num_transcripts; j++) {
                
                double norm_probability = probabilities_unnormalised[j]/normalisation_constant;
                
                if (norm_probability >= double_underflow) {
                
                    binomial_distribution_t binom_dist(num_samples, norm_probability/divider);
                    binomial_sampler_t binom_sampler(mt_rng_pt, binom_dist);
                    
                    int num_maps = binom_sampler();
                    
                    if (num_maps > 0) {

                        map_counts.addToCount(num_maps, j);                    
                        map_counts.addToPlus(j);
                    }
                    
                    num_samples -= num_maps;
                }
                
                divider -= norm_probability;            
            }
            
            assert (num_samples == 0);
        
        }
        
	}

	map_counts.fetchNullSet();

	return map_counts;

}



