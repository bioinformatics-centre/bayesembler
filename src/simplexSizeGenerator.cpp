
/*
simplexSizeGenerator.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/lexical_cast.hpp>
#include <simplexSizeGenerator.h>


typedef boost::random::gamma_distribution<> gamma_distribution_t;
typedef boost::random::variate_generator<boost::random::mt19937*, boost::random::gamma_distribution<> > gamma_sampler_t;
typedef boost::random::mt19937* mt_rng_pt_t;
typedef boost::random::uniform_01<boost::random::mt19937*> uniform_01_sampler_t;

int SimplexSizeGenerator::sampleSimplexSize(double pi, double gamma, int count_plus_size) {
	
	assert(count_plus_size <= num_transcripts);
	assert(pi >= double_underflow);
	assert(pi <= double_almost_one);
	assert(gamma > double_underflow);

	// Init binomial distribution over simplex sizes
	vector <double> simplex_prob_vector;
    simplex_prob_vector.reserve(num_transcripts - count_plus_size + 1);
	
	// Cardinality of equivalence class of size |s+| is one
	double cardinal_eq_z_log = 0;
	
	// Calculate probability of member of equivalence class
	double prob_z_log = count_plus_size*log(pi) + (num_transcripts-count_plus_size)*log(1-pi);
	
	// Probability of assignment given the binary vector
	double prob_t_log = boost::math::lgamma(count_plus_size*gamma) - boost::math::lgamma(num_fragments + count_plus_size*gamma);
	
	// Full probability of the binary vector		
	double prob_eq_z_log = cardinal_eq_z_log + prob_z_log + prob_t_log;
	double row_sum = prob_eq_z_log;
	
	simplex_prob_vector.push_back(row_sum);
    	
	for (int i = count_plus_size + 1; i < num_transcripts + 1; i++) {
		
		// Calculate cardinality of equivalence class
		cardinal_eq_z_log = boost::math::lgamma(num_transcripts-count_plus_size+1)-(boost::math::lgamma(i-count_plus_size+1)+boost::math::lgamma(num_transcripts - i + 1));
		
		// Calculate probability of member of equivalence class
		prob_z_log = i*log(pi) + (num_transcripts-i)*log(1-pi);
		
		// Probability of assignment given the binary vector
		prob_t_log = boost::math::lgamma(i*gamma) - boost::math::lgamma(num_fragments + i*gamma);
		
		// Full probability of the binary vector
		prob_eq_z_log = cardinal_eq_z_log + prob_z_log + prob_t_log;
		
		row_sum += log(1 + exp(prob_eq_z_log - row_sum));
		simplex_prob_vector.push_back(row_sum);
                
    	if (double_compare(simplex_prob_vector.back(), *(simplex_prob_vector.rbegin() + 1))) {
            
            break;     
        }
	}
	
	// Row-normalise and transform back from log-space
	for (int i = 0; i < simplex_prob_vector.size(); i++) {
		
	    simplex_prob_vector[i] = exp(simplex_prob_vector[i] - row_sum);
        
	}
    
    assert (simplex_prob_vector.back() > double_almost_one);
    
	uniform_01_sampler_t sample_uniform_01(mt_rng_pt);
	
	int b = int(upper_bound(simplex_prob_vector.begin(), simplex_prob_vector.end(), sample_uniform_01()) - simplex_prob_vector.begin()) + count_plus_size;
    
	return b;		
}


bool SimplexSizeGenerator::isFixed() {
	
	return is_fixed;
	
}

/* SIMPLEX-SIZE GENERATOR CLASS */
FixedBinomialFixedGammaSimplexSizeGenerator::FixedBinomialFixedGammaSimplexSizeGenerator(double pi_in, double gamma_in, int num_transcripts_in, mt_rng_pt_t mt_rng_pt_in, int num_fragments_in) {

    // Set members 
    pi = pi_in;
    gamma = gamma_in;
    num_transcripts = num_transcripts_in;
    is_fixed = true;
    mt_rng_pt = mt_rng_pt_in;
    num_fragments = num_fragments_in;
    
    simplex_prob_matrix = vector<vector<double> >(num_transcripts, vector<double>());
    
    // Loop over count plus sizes 
    for (int i=0; i < num_transcripts; i++) {
        
        int count_plus_size = i + 1;
        
    	// Init binomial distribution over simplex sizes
    	vector < double> simplex_prob_vector;
	
    	// Cardinality of equivalence class of size |s+| is zero
    	double cardinal_eq_z_log = 0;
	
    	// Calculate probability of member of equivalence class
    	double prob_z_log = count_plus_size*log(pi) + (num_transcripts-count_plus_size)*log(1-pi);
	
    	// Probability of assignment given the binary vector
    	double prob_t_log = boost::math::lgamma(count_plus_size*gamma) - boost::math::lgamma(num_fragments + count_plus_size*gamma);
	
    	// Full probability of the binary vector		
    	double prob_eq_z_log = cardinal_eq_z_log + prob_z_log + prob_t_log;
    	double row_sum = prob_eq_z_log;
	
    	simplex_prob_vector.push_back(row_sum);
	
    	for (int j = count_plus_size + 1; j < num_transcripts + 1; j++) {
		
    		// Calculate cardinality of equivalence class
    		cardinal_eq_z_log = boost::math::lgamma(num_transcripts-count_plus_size+1)-(boost::math::lgamma(j-count_plus_size+1)+boost::math::lgamma(num_transcripts - j + 1));
		
    		// Calculate probability of member of equivalence class
    		prob_z_log = j*log(pi) + (num_transcripts-j)*log(1-pi);
		
    		// Probability of assignment given the binary vector
    		prob_t_log = boost::math::lgamma(j*gamma) - boost::math::lgamma(num_fragments + j*gamma);
		
    		// Full probability of the binary std::vector<char> v;
    		prob_eq_z_log = cardinal_eq_z_log + prob_z_log + prob_t_log;
		    
            row_sum += log(1 + exp(prob_eq_z_log - row_sum));
    		simplex_prob_vector.push_back(row_sum);
            
            if (double_compare(simplex_prob_vector.back(), *(simplex_prob_vector.rbegin() + 1))) {

                break;     
            }
    	}
    
    	// Row-normalise and transform back from log-space
    	for (int j = 0; j < simplex_prob_vector.size(); j++) {
		
    	    simplex_prob_vector[j] = exp(simplex_prob_vector[j] - row_sum);
        
    	}

        assert (simplex_prob_vector.back() > double_almost_one);
        
        simplex_prob_matrix[i] = simplex_prob_vector;
        
    }
}

string FixedBinomialFixedGammaSimplexSizeGenerator::getParameterString(){
        
    string parameter_str;       
    parameter_str += "pi";
    parameter_str += boost::lexical_cast<string>(pi);
    parameter_str += "gamma";
    parameter_str += boost::lexical_cast<string>(gamma);
    
    return parameter_str;
}

// Initialise simplex probability matrix
pair<int, double> FixedBinomialFixedGammaSimplexSizeGenerator::initSimplexSize(int count_plus_size, double gamma) {
            
    int b = sampleSimplexSize(pi, gamma, count_plus_size);
    
    return pair<int, double>(b, pi);
}

// Samples a simplex size from the simplex size probability matrix
pair<int, double> FixedBinomialFixedGammaSimplexSizeGenerator::generateSimplexSize(int expression_plus_size, int count_plus_size, double gamma) {
           
    int b = sampleSimplexSize(pi, gamma, count_plus_size);
    
    return pair<int, double>(b, pi);
}

int FixedBinomialFixedGammaSimplexSizeGenerator::sampleSimplexSize(double pi, double gamma, int count_plus_size) {
	
    uniform_01_sampler_t sample_uniform_01(mt_rng_pt);
    
    int b = int(upper_bound(simplex_prob_matrix[count_plus_size-1].begin(), simplex_prob_matrix[count_plus_size-1].end(), sample_uniform_01()) - simplex_prob_matrix[count_plus_size-1].begin()) + count_plus_size;
        
    return b;        

}


FixedBinomialSimplexSizeGenerator::FixedBinomialSimplexSizeGenerator(double pi_in, int num_transcripts_in, mt_rng_pt_t mt_rng_pt_in, int num_fragments_in) {

	// Set members 
	pi = pi_in;
	num_transcripts = num_transcripts_in;
	is_fixed = true;
	mt_rng_pt = mt_rng_pt_in;
	num_fragments = num_fragments_in;	    
}

string FixedBinomialSimplexSizeGenerator::getParameterString(){
    	
    string parameter_str;       
	parameter_str += "pi";
    parameter_str += boost::lexical_cast<string>(pi);
    
    return parameter_str;
}

// Initialise simplex probability matrix
pair<int, double> FixedBinomialSimplexSizeGenerator::initSimplexSize(int count_plus_size, double gamma) {
	
		
	int b = sampleSimplexSize(pi, gamma, count_plus_size);
			
	return pair<int, double>(b, pi);
}

// Samples a simplex size from the simplex size probability matrix
pair<int, double> FixedBinomialSimplexSizeGenerator::generateSimplexSize(int expression_plus_size, int count_plus_size, double gamma) {
		
	int b = sampleSimplexSize(pi, gamma, count_plus_size);
			
	return pair<int, double>(b, pi);
}

BetaBinomialSimplexSizeGenerator::BetaBinomialSimplexSizeGenerator(double alpha_in, double beta_in, int num_transcripts_in, mt_rng_pt_t mt_rng_pt_in, int num_fragments_in, int slice_iterations_in, double slice_window_size_in) {
	
	alpha = alpha_in;
	beta = beta_in;
	num_transcripts = num_transcripts_in;
	is_fixed = false;
	mt_rng_pt = mt_rng_pt_in;
	num_fragments = num_fragments_in;
	slice_iterations = slice_iterations_in;
    slice_window_size = slice_window_size_in;
	
}

string BetaBinomialSimplexSizeGenerator::getParameterString(){
    
	string parameter_str;       
	parameter_str += "alpha";
    parameter_str += boost::lexical_cast<string>(alpha);
	parameter_str += "_";
	parameter_str += "beta";
    parameter_str += boost::lexical_cast<string>(beta);
    
    return parameter_str;
}

pair<int, double> BetaBinomialSimplexSizeGenerator::initSimplexSize(int count_plus_size, double gamma) {
			
	// Sample pi from beta-prior
	gamma_distribution_t gamma_dist_alpha(alpha,1);
	gamma_distribution_t gamma_dist_beta(beta,1);
	
	gamma_sampler_t sample_gamma_alpha(mt_rng_pt, gamma_dist_alpha);
	gamma_sampler_t sample_gamma_beta(mt_rng_pt, gamma_dist_beta);
	
	double sample_alpha = sample_gamma_alpha();
	double sample_beta = sample_gamma_beta();
	
	double pi = sample_alpha / (sample_alpha + sample_beta);
		
	if (pi > double_almost_one) {
		
		pi = double_almost_one;
	}
    
	if (pi < double_precision) {
		
		pi = double_precision;	
	}
	
	int b = sampleSimplexSize(pi, gamma, count_plus_size);
			
	return pair<int, double>(b, pi);
	
}

pair<int, double> BetaBinomialSimplexSizeGenerator::generateSimplexSize(int expression_plus_size, int count_plus_size, double gamma) {
	
	double pi = samplePi(expression_plus_size);
	
	if (pi > double_almost_one) {
		
		pi = double_almost_one;	
	}
    
	if (pi < double_precision) {
		
		pi = double_precision;
	}
		
	int b = sampleSimplexSize(pi, gamma, count_plus_size);	
	
	return pair<int, double>(b, pi);
}

double BetaBinomialSimplexSizeGenerator::posteriorPiLogDensity(int expression_plus_size, double pi) {
        
	return ((alpha + expression_plus_size - 1)*log(pi) + (beta + num_transcripts - expression_plus_size - 1)*log(1-pi) - log(1 - exp(num_transcripts*log(1-pi))));
    
}

double BetaBinomialSimplexSizeGenerator::samplePi(int expression_plus_size) {
		
	// Init gamma and uniform sampler
	uniform_01_sampler_t sample_uniform_01(mt_rng_pt);
	
	double pi_current = sample_uniform_01();
	double pi = pi_current;
    		
	if (pi_current < double_precision) {
		
		pi_current = double_precision;
		
	}
	
	if (pi_current > double_almost_one) {
		
		pi_current = double_almost_one;
		
	}
		
	// Output all samples for convergence assessment
	// ofstream pi_slice_out("pi_slice_out.txt", ios::app);
	
	for (int i=0; i < slice_iterations; i++) {
		
		// Sample height	
		double y = posteriorPiLogDensity(expression_plus_size, pi_current) + log(1-sample_uniform_01());
				
		// Find slice by "step-out"
		double left = pi_current - sample_uniform_01() * slice_window_size;
		double right = left + slice_window_size;
		
		int j = 1;
		int k = 1;
		
		// Truncate distribution at zero
		if (left < double_precision) {			
			left = double_precision;
			j = 0;
		}
		
		// Truncate distribution at one
		if (right > double_almost_one) {
			right = double_almost_one;
			k = 0;
		}
		
		// Expand window to the left		
		while (j == 1 && y < (posteriorPiLogDensity(expression_plus_size, left))) {
			left = left - slice_window_size;
            
			if (left < double_precision) {
				left = double_precision;
				break;
			}
		}
		
		// Expand window to the right
		while (k == 1 && y < (posteriorPiLogDensity(expression_plus_size, right))) {
			right = right + slice_window_size;
            
			if (right > double_almost_one) {				
				right = double_almost_one;
				break;
			}
		}
        
		// Sample from the window and step-in window boundaries until in slice
		pi = sample_uniform_01()*(right-left) + left;
		
		while ( y >= (posteriorPiLogDensity(expression_plus_size, pi))) {
			            
            if (pi < pi_current) {
				
				left = pi;
				pi = sample_uniform_01()*(right-left) + left;
				
			} else {
				
				right = pi;
				pi = sample_uniform_01()*(right-left) + left;			
			}                
		}
		
		pi_current = pi;
		// pi_slice_out << pi << endl;
	}
	
	// pi_slice_out.close();
	
	return pi;
}