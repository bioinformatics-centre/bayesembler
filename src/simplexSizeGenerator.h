
/*
simplexSizeGenerator.h - This file is part of the Bayesembler (v1.1.1)


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


#ifndef __bayesembler__simplexSizeGenerator_h
#define __bayesembler__simplexSizeGenerator_h

#include <boost/math/distributions/beta.hpp>
#include <string>
#include <vector>
#include <utils.h>
#include <boost/random/mersenne_twister.hpp>


using namespace std;

typedef boost::random::mt19937* mt_rng_pt_t;

class SimplexSizeGenerator {
    
    public:
		virtual string getParameterString() = 0;
		virtual pair<int, double> initSimplexSize(int, double) = 0;
        virtual pair<int, double> generateSimplexSize(int, int, double) = 0;
		bool isFixed();
		int sampleSimplexSize(double, double, int);
		virtual ~SimplexSizeGenerator(){};
	
	protected:
		int num_transcripts;
		bool is_fixed;
		mt_rng_pt_t mt_rng_pt;
		int num_fragments;
		bool is_normalised;

};


class FixedBinomialFixedGammaSimplexSizeGenerator: public SimplexSizeGenerator {
        
    public:
        FixedBinomialFixedGammaSimplexSizeGenerator(double, double, int, mt_rng_pt_t, int);
		string getParameterString();
		pair<int, double> initSimplexSize(int, double);
        pair<int, double> generateSimplexSize(int, int, double);
		int sampleSimplexSize(double, double, int);

    private:
		double pi;
        double gamma;
        vector<vector<double> > simplex_prob_matrix;
		vector<double> simplex_size_prior;
};


class FixedBinomialSimplexSizeGenerator: public SimplexSizeGenerator {
        
    public:
        FixedBinomialSimplexSizeGenerator(double, int, mt_rng_pt_t, int);
		string getParameterString();
		pair<int, double> initSimplexSize(int, double);
        pair<int, double> generateSimplexSize(int, int, double);
		
    private:
		double pi;
		vector<double> simplex_size_prior;
};

class BetaBinomialSimplexSizeGenerator: public SimplexSizeGenerator {
        
    public:
		BetaBinomialSimplexSizeGenerator(double, double, int, mt_rng_pt_t, int, int, double);
		string getParameterString();
		pair<int, double> initSimplexSize(int, double);
        pair<int, double> generateSimplexSize(int, int, double);

    private:
        double alpha;
		double beta;
		boost::math::beta_distribution<> pi_prior;
		double samplePi(int);
		double posteriorPiLogDensity(int, double);
		int slice_iterations;
    	double slice_window_size;
};

#endif 