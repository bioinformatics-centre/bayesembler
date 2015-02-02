
/*
gibbsSampler.cpp - This file is part of the Bayesembler (v1.2.0)


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


#include <gibbsSampler.h>
#include <valueContainer.h>
#include <sstream>
#include <fstream>
#include <vector>

GibbsSampler::GibbsSampler(SimplexSizeGenerator* simplex_size_generator_in, GammaGenerator* gamma_generator_in, ExpressionGenerator* expression_generator_in, AssignmentGenerator* assignment_generator_in, tr1::unordered_map<int,pair<int,double> > new_idx_to_old_idx_effective_length_in, int total_num_fragments_in, int num_fragments_in, int num_transcripts_in, mt_rng_pt_t mt_rng_pt_in, bool assignment_minimum_in) {

	simplex_size_generator = simplex_size_generator_in;
	gamma_generator = gamma_generator_in;
	expression_generator = expression_generator_in;
	assignment_generator = assignment_generator_in;
    new_idx_to_old_idx_effective_length = new_idx_to_old_idx_effective_length_in;
    total_num_fragments = total_num_fragments_in;
    num_fragments = num_fragments_in;
	num_transcripts = num_transcripts_in;
    mt_rng_pt = mt_rng_pt_in;
	assignment_minimum = assignment_minimum_in;
}

void GibbsSampler::setAssignmentGenerator(AssignmentGenerator * assignment_generator_in) {
    
	assignment_generator = assignment_generator_in;
}

void GibbsSampler::runSampler(GibbsOutput* gibbs_output, int gibbs_burn_in, int gibbs_samples, stringstream& log_stream) {
    
    if (num_transcripts == 1) { 

        vector<int> single_binary(1,0);
        vector<double> single_fpkm(1, ((num_fragments * 1e9)/(total_num_fragments * new_idx_to_old_idx_effective_length[0].second)));

        log_stream << "\n[" << getLocalTime() << "] " << "Single candidate - skipping gibbs sampling" << endl;
        
        for (int i=0; i < gibbs_samples; i++) {
            
            gibbs_output->addSample(single_binary, single_fpkm);
        }


            
    } else {
    
        vector<double> effective_lengths(num_transcripts, 0);
    
        for (int i=0; i < num_transcripts; i++) {
        
            effective_lengths[i] = new_idx_to_old_idx_effective_length[i].second;
        }
    
    	// Initialise sampler    
    	CountValueContainer map_counts(num_transcripts);
    	ExpressionValueContainer expression(num_transcripts);
		
    	if (assignment_minimum) {
	
    		map_counts = assignment_generator->initAssignmentMinimum(log_stream);
		
    	} else {
		
    		map_counts = assignment_generator->initAssignment(log_stream);		
    	}
        
        vector<int> count_debug = map_counts.getCounts();
        
    	double gamma = gamma_generator->initGamma();
    	pair<int, double> simplex_size = simplex_size_generator->initSimplexSize(map_counts.getPlusSize(), gamma);
        
    	// Burn-in sampler 
    	log_stream << "\n[" << getLocalTime() << "] " << "Performing " << gibbs_burn_in << " burn-in iterations" << endl;
    
        int print_frequency = 2000;
    
        for (int i = 0; i < gibbs_burn_in; i++) {
        
    		expression = expression_generator->generateExpression(simplex_size.first, gamma, map_counts);
            map_counts = assignment_generator->generateAssignment(expression);
            simplex_size = simplex_size_generator->generateSimplexSize(expression.getPlusSize(), map_counts.getPlusSize(), gamma);
    		gamma = gamma_generator->generateGamma(expression);
            
            vector<int> plus = expression.getPlus(); 

            if (i%print_frequency == 0 && i>0) {
            
                log_stream << "[" << getLocalTime() << "] " << "Completed " << i << " burn-in iterations" << endl;    
            }
    	}
    
        log_stream << "[" << getLocalTime() << "] " << "Completed " << gibbs_burn_in << " burn-in iterations" << endl;
	
    	// Run main Gibbs loop	
    	log_stream << "\n[" << getLocalTime() << "] " << "Performing " << gibbs_samples << " gibbs iterations" << endl;


    	for (int i = 0; i < gibbs_samples; i++) {
	    			
            expression = expression_generator->generateExpression(simplex_size.first, gamma, map_counts);
            map_counts = assignment_generator->generateAssignment(expression);
            simplex_size = simplex_size_generator->generateSimplexSize(expression.getPlusSize(), map_counts.getPlusSize(), gamma);
            gamma = gamma_generator->generateGamma(expression);            
            
            vector<double> fpkm = expression.getFPKM(effective_lengths, num_fragments, total_num_fragments);          
            vector<int> plus = expression.getPlus(); 
            gibbs_output->addSample(plus, fpkm);     


    		if (i%print_frequency == 0 && i>0) {
		
                log_stream << "[" << getLocalTime() << "] " << "Completed " << i << " gibbs iterations" << endl;
            
            }
    	}
    
        log_stream << "[" << getLocalTime() << "] " << "Completed " << gibbs_samples << " gibbs iterations" << endl;
        
        count_debug = map_counts.getCounts();
    
    }
}

Ensemble GibbsSampler::reestimateEnsembleExpression(Ensemble ensemble, int samples, stringstream& log_stream) {
    
    log_stream << "[" << getLocalTime() << "] " << "Re-estimating expression levels for transcripts " << ensemble.indices << endl;
    
    if (num_transcripts == 1) { 
    
        return ensemble;
                        
    } else {
                
        vector<double> effective_lengths(num_transcripts, 0);
    
        for (int i=0; i < num_transcripts; i++) {
        
            effective_lengths[i] = new_idx_to_old_idx_effective_length[i].second;
        }
    
        Eigen::VectorXd mean_cumulant(ensemble.indices.size());
        Eigen::VectorXd var_cumulant(ensemble.indices.size());
    	
        CountValueContainer map_counts(num_transcripts);
    	ExpressionValueContainer expression(num_transcripts);
        
        map_counts = assignment_generator->initEnsembleAssignment(ensemble.indices);
        double gamma = gamma_generator->initGamma();
                
        // Burn-in reestimation sampler
        for (int i=0; i < 1000; i++) {

            expression = expression_generator->generateEnsembleExpression(ensemble.indices, gamma, map_counts);
            map_counts = assignment_generator->generateAssignment(expression);
            gamma = gamma_generator->generateGamma(expression);

        }
        
        // Collect samples
        for (int i=0; i < samples; i++) {
                
            expression = expression_generator->generateEnsembleExpression(ensemble.indices, gamma, map_counts);
            map_counts = assignment_generator->generateAssignment(expression);
            gamma = gamma_generator->generateGamma(expression);


            Eigen::VectorXd fpkm_plus = expression.getFPKMplus(effective_lengths, num_fragments, total_num_fragments);
                        
            if (i==0) {
                
                mean_cumulant = fpkm_plus;
                var_cumulant = Eigen::VectorXd::Zero(ensemble.indices.size());
            
            } else {
                
                Eigen::VectorXd delta = fpkm_plus.array() - mean_cumulant.array();        
                mean_cumulant.array() += delta.array()/(i+1);
                var_cumulant.array() += delta.array() * (fpkm_plus.array() - mean_cumulant.array());
            }
        }
        
        ensemble.means = mean_cumulant;
        Eigen::VectorXd variance = var_cumulant.array()/(samples-1);
        ensemble.sd = variance.array().sqrt();
        
        return ensemble;
    }    
} 




































