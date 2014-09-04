
/*
fragmentLengthModel.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include "fragmentLengthModel.h"
#include <math.h>
#include <fstream>
#include <sstream>

#include "boost/math/constants/constants.hpp"
#include <boost/filesystem.hpp>

GaussianFragmentLengthModel::GaussianFragmentLengthModel() {}

GaussianFragmentLengthModel::GaussianFragmentLengthModel(double mean_in, double sd_in) {

	mean = mean_in;
	sd = sd_in;
	var = pow(sd, 2);

	norm_const = 1/(sd*sqrt(2*boost::math::constants::pi<double>()));

	cout << "[" << getLocalTime() << "] Using Gaussian fragment length distribution with parameters: Mean=" << mean << " and SD=" << sd << endl;

}

GaussianFragmentLengthModel::GaussianFragmentLengthModel(map<uint,uint> fragment_lengths) {

	pair<double,double> parameters = estimateParameters(fragment_lengths);

	mean = parameters.first;
	sd = 1.4826 * parameters.second;

	var = pow(sd, 2);
	norm_const = 1/(sd*sqrt(2*boost::math::constants::pi<double>()));

	cout << "[" << getLocalTime() << "] Using Gaussian fragment length distribution with parameters: Mean=" << mean << " and SD=" << sd << endl;
}


double GaussianFragmentLengthModel::pmf(uint fragment_length) {

	return norm_const * exp(-(pow((double) fragment_length - mean, 2)/(2*var)));
}

template <typename DataType>
double GaussianFragmentLengthModel::estimateMedian(map<DataType,uint> observations) {

	int observation_count = 0;

	for (typename map<DataType,uint>::iterator observation_iter = observations.begin(); observation_iter != observations.end(); observation_iter++) {

		observation_count += observation_iter->second;
	}

	assert(observation_count > 0);

	// Find median observation index for uneven observation number and first median observation for even observation number
	int median_observation_idx;

	if (observation_count % 2 == 0) {

		median_observation_idx = (int) round((double) observation_count/2.0);

	} else {

		median_observation_idx = (int) ceil((double) observation_count/2.0);
	}

	// Iterate to median observation index
	typename map<DataType,uint>::iterator observation_iter = observations.begin();
	int observation_count_running = observation_iter->second;

	while (observation_count_running < median_observation_idx) {

		observation_iter++;
		observation_count_running += observation_iter->second;
	} 

	if (observation_count % 2 == 0 and observation_count_running == median_observation_idx) {

		DataType median = observation_iter->first; 
		observation_iter++; 
		assert(observation_iter != observations.end());
		median += observation_iter->first; 
		return (double) median / 2;

	} else {

		assert(observation_count_running >= median_observation_idx);

		return (double) observation_iter->first; 
	}
}

pair<double,double> GaussianFragmentLengthModel::estimateParameters(map<uint,uint> fragment_lengths) {

	double median = estimateMedian(fragment_lengths);

	map<double,uint> abs_distances;
	int observation_count = 0;

	for (map<uint,uint>::iterator fragment_lengths_iter = fragment_lengths.begin(); fragment_lengths_iter != fragment_lengths.end(); fragment_lengths_iter++) {

		assert(fragment_lengths_iter->first > 0);
		observation_count += fragment_lengths_iter->second;

		double distance = abs((double) fragment_lengths_iter->first - median); 

		if (abs_distances.count(distance) == 0) {

			abs_distances[distance] = fragment_lengths_iter->second;

		} else {

			abs_distances[distance] += fragment_lengths_iter->second;

		}	
	}

	double median_abs_deviation = estimateMedian(abs_distances);

	if (observation_count < 10000) {

		cout << "\n[" << getLocalTime() << "] WARNING: Only " << observation_count << " observation(s) were used for estimation of the fragment length distribution!\n" << endl;
	}

	cout << "[" << getLocalTime() << "] Estimated fragment length \"median\"=" << median << " and \"median absolute deviation\"=" << median_abs_deviation << " using " << observation_count << " observations" << endl;

	return pair<double,double> (median, median_abs_deviation);
}

EmpiricalFragmentLengthModel::EmpiricalFragmentLengthModel(string distribution_filename) {

	assert(boost::filesystem::exists(distribution_filename));  

	ifstream distribution_file(distribution_filename.c_str());
    assert(distribution_file.is_open());
    
    string distribution_line;    
    
    uint fragment_length;
    double probability;
    uint line_count = 0;

    while (getline(distribution_file, distribution_line)) {

    	line_count++;
        stringstream linestream(distribution_line);

	    // Read the integers using the operator >>
	    linestream >> fragment_length >> probability;

	    assert(distribution.count(fragment_length) == 0);
	    distribution[fragment_length] = probability;
    }

    cout << "[" << getLocalTime() << "] Loaded empirical fragment length distribution containing " << line_count << " fragment length values - the remaining will have point-zero probability " << endl;
}

double EmpiricalFragmentLengthModel::pmf(uint fragment_length) {

	if (distribution.count(fragment_length)) {

		return distribution[fragment_length];
	
	} else {

		return 0;
	}
}
