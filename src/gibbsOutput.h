
/*
gibbsOutput.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef __bayesembler__gibbsOutput_h
#define __bayesembler__gibbsOutput_h

#include "utils.h"
#include <vector>
#include <tr1/unordered_map>
#include <Eigen/Dense>
#include <string>

using namespace std;

class GibbsOutput {


    public:
        
        GibbsOutput(int, int);
        void addSample(vector<int>, vector<double>);
        Ensemble getTopMarginals(double);
        string writeEnsembleGtf(Ensemble , vector<Candidate>, tr1::unordered_map<int,pair<int,double> >, int, double);
        string writeCandidateGtf(vector<Candidate>, tr1::unordered_map<int,pair<int,double> >, int);
                                
    private:    
        
        int num_iterations;
        int num_transcripts;
        
        Eigen::VectorXi marginal_occurence;
        Eigen::VectorXd marginal_mean_cumulant;
        Eigen::VectorXd marginal_var_cumulant;
    
};

#endif
