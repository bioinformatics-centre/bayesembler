
/*
assignmentGenerator.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef __bayesembler__assignmentGenerator_h
#define __bayesembler__assignmentGenerator_h

#include <valueContainer.h>
#include <utils.h>
#include <string>
#include <tr1/unordered_set>
#include <boost/random/mersenne_twister.hpp>
#include <Eigen/Dense>
#include <vector>

using namespace std;

typedef boost::random::mt19937* mt_rng_pt_t;
	
class AssignmentGenerator {
	
	public:
		AssignmentGenerator(mt_rng_pt_t, Eigen::MatrixXd *, vector <int> *);

        void calcMinimumReadCover(stringstream&);
        int getMinimumReadCoverSize();
		CountValueContainer initAssignment(stringstream&);
		CountValueContainer	initAssignmentMinimum(stringstream&);
        CountValueContainer initEnsembleAssignment(vector<int> indices);
		CountValueContainer generateAssignment(ExpressionValueContainer);
			
	private:
		
        Eigen::MatrixXd * collapsed_map_eigen;
        vector <int> * collapsed_count_vec;
    
		int num_transcripts;
		mt_rng_pt_t mt_rng_pt;
		tr1::unordered_set<int> minimum_set;
                
};
#endif