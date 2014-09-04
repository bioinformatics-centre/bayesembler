
/*
sequencingModel.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include <math.h>
#include <sequencingModel.h>


SequencingModel::SequencingModel(FragmentLengthModel * fragment_length_model_pt, int max_transcript_length, stringstream& log_stream) {

	log_stream << "\n[" << getLocalTime() << "] " << "Tabulating fragment length and 3'-position probabilities up to transcript lengths of " << max_transcript_length << endl;
    
	// 3'-position cumulated probability vector
	three_prime_cum_prob = vector< double > (max_transcript_length, 0);
	frag_length_cum_prob = vector< double >  (max_transcript_length, 0);
	frag_length_prob = vector< double > (max_transcript_length, 0);
	
	// Initialise for fragment length = 1
	// frag_length_prob[0] = pdf(frag_length_dist, 1);
	frag_length_prob[0] = fragment_length_model_pt->pmf(1);
	frag_length_cum_prob[0] = frag_length_prob[0];
	three_prime_cum_prob[0] = frag_length_prob[0];
	
	// Cumulate 
	for (int i=1; i < max_transcript_length; i++) {
	
		// frag_length_prob[i] = pdf(frag_length_dist,(i+1)); 
		frag_length_prob[i] = fragment_length_model_pt->pmf(i+1); 
		frag_length_cum_prob[i] = frag_length_prob[i] + frag_length_cum_prob[i-1];
		three_prime_cum_prob[i] = frag_length_cum_prob[i] + three_prime_cum_prob[i-1];	
	}
}

double SequencingModel::calculateThreePrimeProb(int three_prime_pos, int tran_length) {
    
	double three_prime_prob = frag_length_cum_prob[three_prime_pos-1]/three_prime_cum_prob[tran_length-1];
    
    return three_prime_prob;
}

double SequencingModel::calculateLengthProb(int fragment_length, int three_prime_pos) {
    
	double length_prob = frag_length_prob[fragment_length-1]/frag_length_cum_prob[three_prime_pos-1];
    
    return length_prob;
}

double SequencingModel::getEffectiveLength(int transcript_length) {
    
    return three_prime_cum_prob[transcript_length-1];
}
