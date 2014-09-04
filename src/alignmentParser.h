
/*
alignmentParser.h - This file is part of the Bayesembler (v1.1.1)


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


#ifndef ALIGNMENT_PARSER_H
#define ALIGNMENT_PARSER_H

#include <string>
#include <utils.h>
#include <sequencingModel.h>
#include <Eigen/Dense>
#include <tr1/unordered_set>
#include <map>
#include <list>

using namespace std;

class AlignmentParser{
    
    public:
    
    	CollapsedMap calculateFragTranProbabilities(vector<Candidate>&, vector<FragmentAlignment*>*, SequencingModel *, bool, stringstream&, tr1::unordered_map<string, int> *);    
        void fetchFragmentLengths(map<uint,uint> *, vector<Candidate> &, vector<FragmentAlignment*>*, stringstream&);
    
    private:

        struct FragmentMatch {
            
            int first_exon_idx; 
            int last_exon_idx;
            int fragment_length;
            int three_prime_position;
        };


        pair<bool, FragmentMatch> matchFragmentToCandidate(FragmentAlignment, Candidate, string);  
};

#endif
    
    