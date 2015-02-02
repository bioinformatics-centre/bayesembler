
/*
gammaGenerator.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef __bayesembler__gammaSampler__
#define __bayesembler__gammaSampler__

#include <string>
#include <valueContainer.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <vector>

using namespace std;


// Class forward declarations

class HyperPrior;
class Prior;
class SliceUnivariateSampler;

/*
    Gamma-sampler class
 */

// Abstract base class for gamma sampler 
class GammaGenerator {
    
    public:
		virtual double initGamma() = 0;
        virtual double generateGamma(ExpressionValueContainer) = 0;
		bool isFixed();	
		
	protected:
		bool fixed;
};

class FixedGammaGenerator: public GammaGenerator{
    
    private:
        
        double gamma;
    
    public:
        FixedGammaGenerator(double);
		double initGamma();
        double generateGamma(ExpressionValueContainer);
};

class SymmetricGammaGenerator: public GammaGenerator{
        
    private:

        double calculateLogDensity(vector<double>&, double);

		typedef boost::random::uniform_01<boost::random::mt19937*> uniform_01_sampler_t;
	    HyperPrior * hyper_prior;
	    Prior * prior;
	    int slice_iterations;
	    double slice_window_size;
	    int slice_max_windows;
		boost::random::mt19937* mt_rng_pt;
    
    public:

		SymmetricGammaGenerator(HyperPrior*, int, double, int, boost::random::mt19937*);
		double initGamma();
        double generateGamma(ExpressionValueContainer expression);
};

///**
// * @brief HyperPrior class
// *
// * Abstract base class for input to @see UnivariateSampler
// */

class HyperPrior {

public:

	virtual double init() = 0;
    virtual double calculateLogDensity(double)=0;
    virtual ~HyperPrior(){};

protected:

	typedef boost::random::mt19937* mt_rng_pt_t;
	mt_rng_pt_t mt_rng_pt;
};

class GammaHyperPrior : public HyperPrior {
    
public:
    GammaHyperPrior(double, double, mt_rng_pt_t);
	double init();
    double calculateLogDensity(double);
    
private:

    double shape;
    double scale;
};


class LogNormalHyperPrior : public HyperPrior {

public:

    LogNormalHyperPrior(double, double, mt_rng_pt_t);
    double init();
    double calculateLogDensity(double);
    
private:
    
    double location;
    double scale;
};



#endif