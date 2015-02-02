
/*
gammaGenerator.cpp - This file is part of the Bayesembler (v1.2.0)


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


#include <math.h>
#include <gammaGenerator.h>
#include <boost/math/constants/constants.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>


/*
 Gamma-sampler class
 */

SymmetricGammaGenerator::SymmetricGammaGenerator(HyperPrior * hyper_prior_in, int slice_iterations_in, double slice_window_size_in, int slice_max_windows_in, boost::random::mt19937* mt_rng_pt_in) {

    hyper_prior = hyper_prior_in;
    slice_iterations = slice_iterations_in;
    slice_window_size = slice_window_size_in;
    slice_max_windows = slice_max_windows_in;
	mt_rng_pt = mt_rng_pt_in;
}


// Init gamma with sample from the hyperprior
double SymmetricGammaGenerator::initGamma() {
		
	return hyper_prior->init();
}

double SymmetricGammaGenerator::generateGamma(ExpressionValueContainer expression) {
    
	vector<double> expression_s_plus;
	
	for (int i = 0; i < expression.getPlusSize(); i++) {
		
		int idx = expression.getPlus(i);
		double value = expression.getValue(idx);
		expression_s_plus.push_back(value);
    }
	
    // Init gamma and uniform sampler
    uniform_01_sampler_t sample_uniform_01(mt_rng_pt);
	    
    double gamma_current = sample_uniform_01();

	if (gamma_current < double_underflow) {
		
		gamma_current = double_underflow;		
	}
	
	double gamma = gamma_current;
		    
	for (int i=0; i < slice_iterations; i++) {
	    
        double y = calculateLogDensity(expression_s_plus, gamma_current) + hyper_prior->calculateLogDensity(gamma_current) + log(1-sample_uniform_01());			
		
        // Find slice by "step-out"
        double left = gamma_current - sample_uniform_01() * slice_window_size;
        double right = left + slice_window_size;
	        
        int j = floor(slice_max_windows*sample_uniform_01());
        int k = slice_max_windows-1-j;
			
		// Truncate at zero
        if (left < double_underflow) {
	            
            left = double_underflow;
            j = 0;
        }

        while (j > 0 && y < (calculateLogDensity(expression_s_plus, left) + hyper_prior->calculateLogDensity(left))) {
			
			left = left - slice_window_size;
            j--;

            if (left < double_underflow) {

                left = double_underflow;
                break;
            }
		}
	
        // Expand window to the right
        while (k > 0 && y < (calculateLogDensity(expression_s_plus, right) + hyper_prior->calculateLogDensity(right))) {
			
			right = right + slice_window_size;
            k--;
        }
			
		// Sample from the window until in slice
		gamma = sample_uniform_01()*(right-left) + left;
	    
		while ( y >= (calculateLogDensity(expression_s_plus, gamma) + hyper_prior->calculateLogDensity(gamma))) {
			  
			if (gamma < gamma_current) {
				
				left = gamma;
				gamma = sample_uniform_01()*(right-left) + left;
								
			} else {
					
				right = gamma;
				gamma = sample_uniform_01()*(right-left) + left;
			}
		}

		gamma_current = gamma;	
	}
	    
    return gamma;		
	
}

double SymmetricGammaGenerator::calculateLogDensity(vector<double>& expression, double gamma) {
    
    int size = expression.size();
    
    // Init prob with normalisation constant
    double prob = boost::math::lgamma(size*gamma) - size*boost::math::lgamma(gamma);
    
    for (int i=0; i < size; i++) {
        prob += (gamma-1)*log(expression[i]);
    }
    
    return prob;
}

FixedGammaGenerator::FixedGammaGenerator(double gamma_in) {
    gamma = gamma_in;
}

double FixedGammaGenerator::initGamma() {
	
	return gamma;
}

double FixedGammaGenerator::generateGamma(ExpressionValueContainer expression) {
      
    return gamma;
}


GammaHyperPrior::GammaHyperPrior(double shape_in, double scale_in, mt_rng_pt_t mt_rng_pt_in) {

    shape = shape_in;
    scale = scale_in;
	mt_rng_pt = mt_rng_pt_in;
};

double GammaHyperPrior::init() {
	
	boost::random::gamma_distribution<> gamma_dist(shape, scale);
	boost::random::variate_generator<boost::random::mt19937*, boost::random::gamma_distribution<> > sample_gamma(mt_rng_pt, gamma_dist);
	
	double gamma = sample_gamma();
	
	return gamma;
}

double GammaHyperPrior::calculateLogDensity(double gamma) {
    	
	double prob = (shape-1)*log(gamma) - (gamma/scale);
	
    return prob;
    
};


LogNormalHyperPrior::LogNormalHyperPrior(double location_in, double scale_in, mt_rng_pt_t mt_rng_pt_in) {
    
    location = location_in;
    scale = scale_in;
	mt_rng_pt = mt_rng_pt_in;
}

double LogNormalHyperPrior::init() {
    
	boost::random::lognormal_distribution<> lognormal_dist(location, scale);
	boost::random::variate_generator<boost::random::mt19937*, boost::random::lognormal_distribution<> > sample_gamma(mt_rng_pt, lognormal_dist);
    
	double gamma = sample_gamma();
	
	return gamma;
}

double LogNormalHyperPrior::calculateLogDensity(double gamma) {
    		
	double prob = - log(gamma) - pow(log(gamma)-location,2)/(2*pow(scale,2));	
    
    return prob;
}



