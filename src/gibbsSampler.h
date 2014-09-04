
/*
gibbsSampler.h - This file is part of the Bayesembler (v1.1.1)


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


#ifndef __bayesembler__gibbsSampler_h
#define __bayesembler__gibbsSampler_h

#include <gibbsOutput.h>
#include <gammaGenerator.h>
#include <assignmentGenerator.h>
#include <expressionGenerator.h>
#include <simplexSizeGenerator.h>
#include <boost/random/mersenne_twister.hpp>
#include <string>
#include <utils.h>


using namespace std;

typedef boost::random::mt19937* mt_rng_pt_t;

// GibbsSampler class
class GibbsSampler{
		
	public:
		GibbsSampler(SimplexSizeGenerator*, GammaGenerator*, ExpressionGenerator*, AssignmentGenerator*, tr1::unordered_map<int,pair<int,double> >, int, int, int, mt_rng_pt_t, bool);
        void setAssignmentGenerator(AssignmentGenerator*);
		void runSampler(GibbsOutput*, int, int, stringstream&);
        Ensemble reestimateEnsembleExpression(Ensemble, int, stringstream&);
        
	private:
		SimplexSizeGenerator* simplex_size_generator;
		GammaGenerator* gamma_generator;
		ExpressionGenerator* expression_generator;
		AssignmentGenerator* assignment_generator;
        tr1::unordered_map<int,pair<int,double> > new_idx_to_old_idx_effective_length;
        int total_num_fragments;
        int num_fragments;
		int num_transcripts;
		mt_rng_pt_t mt_rng_pt;
		bool assignment_minimum;
};

#endif