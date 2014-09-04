
/*
expressionGenerator.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include <expressionGenerator.h>
#include <boost/random/uniform_01.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/discrete_distribution.hpp>

typedef boost::random::mt19937* mt_rng_pt_t;
typedef boost::random::uniform_01<boost::random::mt19937*> uniform_01_sampler_t;
typedef boost::random::gamma_distribution<> gamma_distribution_t;
typedef boost::random::variate_generator<boost::random::mt19937*, boost::random::gamma_distribution<> > gamma_sampler_t;
typedef boost::random::uniform_int_distribution<> uniform_sampler_t;

ExpressionGenerator::ExpressionGenerator(int num_fragments_in, int num_transcripts_in, mt_rng_pt_t mt_rng_pt_in) {
	
	num_fragments = num_fragments_in;
	num_transcripts = num_transcripts_in;
	mt_rng_pt = mt_rng_pt_in;
}

ExpressionValueContainer ExpressionGenerator::generateEnsembleExpression(vector<int> indices, double gamma, CountValueContainer map_counts) {
	
	// Init container
	ExpressionValueContainer expression(num_transcripts); 
	
	// Sample gammas for s-plus
	double norm_const_expression = 0;
	    
	for (int i = 0; i < indices.size(); i++) {
	        
	    int trans_id = indices[i];
	
	    gamma_distribution_t gamma_dist((map_counts.getCount(trans_id) + gamma),1);
        gamma_sampler_t sample_gamma(mt_rng_pt, gamma_dist);
		
        double gamma_sample = sample_gamma();
        
	    norm_const_expression += gamma_sample;
	    expression.setBinaryOn(trans_id);
	    expression.setValue(gamma_sample, trans_id);
        expression.addToPlus(trans_id);
	}
	
	expression.normalise(norm_const_expression);
	
	return expression;
}

ExpressionValueContainer ExpressionGenerator::generateExpression(int b, double gamma, CountValueContainer map_counts) {
	
	// Init container
	ExpressionValueContainer expression(num_transcripts); 
	
	// Sample gammas for s-plus
	double norm_const_expression = 0;
	
	for (int i = 0; i < map_counts.getPlusSize(); i++) {
	        
	    int trans_id = map_counts.getPlus(i);
	
	    gamma_distribution_t gamma_dist((map_counts.getCount(trans_id) + gamma),1);
        gamma_sampler_t sample_gamma(mt_rng_pt, gamma_dist);
		
        double gamma_sample = sample_gamma();
        
	    norm_const_expression += gamma_sample;
	    expression.setBinaryOn(trans_id);
	    expression.setValue(gamma_sample, trans_id);
        expression.addToPlus(trans_id);
	}
	
    gamma_distribution_t gamma_dist_base (gamma, 1);
    gamma_sampler_t sample_gamma_base (mt_rng_pt, gamma_dist_base);

	// Sample gammas for the expanded simplex
	for (int i = map_counts.getPlusSize(); i < b ; i++) {
	       
	    uniform_sampler_t sample_trans(0,(map_counts.getNullSize()-1));
	        
	    int s_null_idx = sample_trans(*mt_rng_pt);
	    int trans_id = map_counts.getNull(s_null_idx);
		
	    assert (map_counts.getCount(trans_id) == 0);
	        	
	    double gamma_base_sample = sample_gamma_base();
	    assert(gamma_base_sample >= double_underflow);

            
        // if (gamma_base_sample < almost_zero) {
            
        //     gamma_base_sample = almost_zero;    
        // }
        
	    norm_const_expression += gamma_base_sample;
	    expression.setBinaryOn(trans_id);
	    expression.setValue(gamma_base_sample, trans_id);
        expression.addToPlus(trans_id);
	    map_counts.eraseNull(s_null_idx);	
	}
	
	expression.normalise(norm_const_expression);
	
	return expression;
}


