
/*
sequencingModel.h - This file is part of the Bayesembler (v1.1.1)


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


#ifndef SEQUENCING_MODEL_H
#define SEQUENCING_MODEL_H

#include <utils.h>
#include <fragmentLengthModel.h>
#include <string>
#include <api/BamAux.h>
#include <vector>

using namespace std;

class SequencingModel {
    
    public:

        SequencingModel(FragmentLengthModel *, int, stringstream&);
        double calculateThreePrimeProb(int, int);
        double calculateLengthProb(int, int);
        double getEffectiveLength(int);

    private:
        
        double lamdba; 
        int max_transcript_length;   
        vector<double> frag_length_prob;
        vector<double> frag_length_cum_prob;
        vector<double> three_prime_cum_prob;

};

#endif