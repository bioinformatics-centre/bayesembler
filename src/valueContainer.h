
/*
valueContainer.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef __bayesembler__valueContainer_h
#define __bayesembler__valueContainer_h

#include <utils.h>
#include <vector>
#include <Eigen/Dense>

using namespace std;

class ValueContainer {

	public:
		vector<int> getPlus();
		vector<int> getNull();
		int getPlus(int);
		int getNull(int);
		int getPlusSize();
		int getNullSize();
		void addToPlus(int);
		void addToNull(int);
		void erasePlus(int index);
		void eraseNull(int index);
	
	protected:
		vector<int> plus;
		vector<int> null;
		int num_transcripts;
	
};


class CountValueContainer : public ValueContainer {
	
	public:	
		CountValueContainer();
		CountValueContainer(int);
		void addToCount(int, int);
		void increment(int);
		int getCount(int);
		vector<int> getCounts();
		void fetchNullSet();
	
	private:
		vector<int> counts;
		
};

class ExpressionValueContainer : public ValueContainer {
	
	public:	
		ExpressionValueContainer();
		ExpressionValueContainer(int);
		void setBinaryOn(int);
		void setValue(double, int);	
		int getBinary(int);
		vector<int> getBinaries();
		double getValue(int);
		vector<double> getValues();
		void normalise(double);
		vector<double> getFPKM(vector<double>, int, int);
        Eigen::VectorXd getFPKMplus(vector <double>, int, int);
		vector<double> getValuesNormalisedToEffectiveLength(vector<double>);
		
	private:
		vector<int> binaries;
		vector<double> values;
		double normalisation_constant;

};
#endif