
/*
valueContainer.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include "valueContainer.h"

vector<int> ValueContainer::getPlus() {
	
	return plus;
}

vector<int> ValueContainer::getNull() {
	
	return null;
}

int ValueContainer::getPlus(int index) {
	
	return plus[index];
}

int ValueContainer::getNull(int index) {
	
	return null[index];
}

int ValueContainer::getPlusSize() {
	
	return plus.size();
}

int ValueContainer::getNullSize() {
	
	return null.size();
}

void ValueContainer::erasePlus(int index) {
	
	plus.erase(plus.begin() + index);	
}

void ValueContainer::eraseNull(int index) {
	
	null.erase(null.begin() + index);	
}

void ValueContainer::addToPlus(int index) {
	
	bool has_index = false;
	for (int i=0; i < plus.size(); i++) {
		
		if (plus[i] == index) {
			
			has_index = true; 
			break;
		}
	}
	
	if (!has_index) {
		plus.push_back(index);
	}
	
}

void ValueContainer::addToNull(int index) {

	bool has_index = false;
	for (int i=0; i < null.size(); i++) {
		
		if (null[i] == index) {
			
			has_index = true; 
			break;
		}
	}
	
	if (!has_index) {
		null.push_back(index);
	}
}

CountValueContainer::CountValueContainer(){};

CountValueContainer::CountValueContainer(int num_transcripts_in) {
	
	counts = vector<int> (num_transcripts_in, 0);
	num_transcripts = num_transcripts_in;

}


void CountValueContainer::addToCount(int count, int index) {
	
	counts[index] += count;
}

void CountValueContainer::increment(int index) {
	
	counts[index]++;
}


int CountValueContainer::getCount(int index) {
	
	return counts[index];
}

vector<int> CountValueContainer::getCounts() {
	
	return counts;
}

void CountValueContainer::fetchNullSet() {
	
	// Sort indices
	sort(plus.begin(), plus.end());
	
    // Check first elements
    int first_plus = plus.front();
    if (first_plus != 0) {
        for (int i = 0; i < first_plus; i++) {
			
			null.push_back(i);
		}
    }

    // Check middle elements
	for (int i = 1; i < plus.size(); i++) {
		
		if ((plus[i-1] + 1) !=  plus[i]) {
			
			for (int j = (plus[i-1]+1); j < plus[i]; j++) {
				
				null.push_back(j);
			}
		}
	}
	
	// Check last elements
	int last_plus = plus.back();
	
	if ( last_plus != (num_transcripts-1)) {
		for (int  i = (last_plus+1); i < num_transcripts; i++) {
			
			null.push_back(i);
		}
	}
    
    assert((plus.size() + null.size()) == num_transcripts);
}

ExpressionValueContainer::ExpressionValueContainer(){};

ExpressionValueContainer::ExpressionValueContainer(int num_transcripts_in) {
	
	binaries = vector<int> (num_transcripts_in, 0);
	values = vector<double> (num_transcripts_in, 0);    
	num_transcripts = num_transcripts_in;
}

void ExpressionValueContainer::setBinaryOn(int index) {
	
	binaries[index] = 1;
}

int ExpressionValueContainer::getBinary(int index) {
	
	return binaries[index];
}

vector<int> ExpressionValueContainer::getBinaries() {
	
	return binaries;
}

void ExpressionValueContainer::setValue(double value, int index) {
	
	values[index] = value;
}

double ExpressionValueContainer::getValue(int index) {
	
	return values[index];
}

vector<double> ExpressionValueContainer::getValues() {
	
	return values;
}

void ExpressionValueContainer::normalise(double normalisation_constant) {
	
	assert(normalisation_constant >= double_underflow);
	
	for (int i=0; i < plus.size(); i++) {
        
		values[plus[i]] /= normalisation_constant;
	}
}

vector<double> ExpressionValueContainer::getFPKM(vector<double> effective_lengths, int num_fragments, int total_num_fragments) {
	
	assert(effective_lengths.size() == num_transcripts);
	
	vector<double> fpkm_values(num_transcripts, 0); 
	
	// Normalise expression to effective length
	for (int i=0; i < plus.size(); i++) {
		int idx = plus[i];
		fpkm_values[idx] = (values[idx] * num_fragments * 1e9) / (effective_lengths[idx] * total_num_fragments);
	}

	return fpkm_values;
}


Eigen::VectorXd ExpressionValueContainer::getFPKMplus(vector<double> effective_lengths, int num_fragments, int total_num_fragments) {
	
	assert(effective_lengths.size() == num_transcripts);
	
    Eigen::VectorXd fpkm_values(plus.size()); 
	
	sort(plus.begin(), plus.end());
    
	// Normalise expression to effective length
	for (int i=0; i < plus.size(); i++) {
        
		int idx = plus[i];
		fpkm_values(i) = (values[idx] * num_fragments * 1e9) / (effective_lengths[idx] * total_num_fragments);
    }

	return fpkm_values;
}

vector<double> ExpressionValueContainer::getValuesNormalisedToEffectiveLength(vector<double> effective_lengths) {
	
	assert(effective_lengths.size() == num_transcripts);
	
	vector<double> normalised_values(num_transcripts, 0); 
	double norm_const = 0;
	
	// Normalise expression to effective length
	for (int i=0; i < plus.size(); i++) {
		int idx = plus[i];
		double norm_value = values[idx]/effective_lengths[idx];
		normalised_values[idx] = norm_value;
		norm_const += norm_value;
	}

	// Re-normalise 
	for (int i=0; i < plus.size(); i++) {
		normalised_values[plus[i]] /= norm_const;
	}
	
	return normalised_values;
}
