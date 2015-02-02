
/*
fragmentLengthModel.h - This file is part of the Bayesembler (v1.2.0)


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


#ifndef FRAGMENT_LENGTH_MODEL_H
#define FRAGMENT_LENGTH_MODEL_H

#include "utils.h"
#include <map>
#include <tr1/unordered_map>
#include <string>

using namespace std;

class FragmentLengthModel {

	public: 
	
		virtual double pmf(uint)=0;
		virtual ~FragmentLengthModel(){};

};

class GaussianFragmentLengthModel : public FragmentLengthModel {

	public: 

		GaussianFragmentLengthModel();
		GaussianFragmentLengthModel(double, double);
		GaussianFragmentLengthModel(map<uint,uint>);
		pair<double,double> estimateParameters(map<uint,uint>);

		template<typename DataType>
		double estimateMedian(map<DataType,uint>);
	
		double pmf(uint);
		
	private: 

		double mean; 
		double sd;
		double var; 
		double norm_const; 

};

class EmpiricalFragmentLengthModel : public FragmentLengthModel {

	public:

		EmpiricalFragmentLengthModel(string);
		double pmf(uint);

	private:

		tr1::unordered_map<uint,double> distribution;
};

#endif