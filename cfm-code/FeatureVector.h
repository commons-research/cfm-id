//*#########################################################################
// # Mass Spec Prediction and Identification of Metabolites
// #
// # FeatureVector.h
// #
// # Description: 	Code for computing features for fragmentations.
// #
// # Copyright (c) 2018
// # All rights reserved.
//
// # This file is part of the cfm-id project.
// # The contents are covered by the terms of the GNU Lesser General Public
// # License, which is included in the file license.txt, found at the root
// # of the cfm source tree.
// #########################################################################*/

#pragma once

#include <iostream>
#include <vector>

// Structure to hold a sparse computed feature vector
typedef unsigned int feature_t;

class FeatureVector {

private:
	std::vector<feature_t> _fv;
	unsigned int _fv_idx;

public:
	FeatureVector() { _fv_idx = 0; };

	FeatureVector(const FeatureVector &old) {
		_fv_idx = old._fv_idx;
		_fv     = old._fv;
	};

	void addFeature(double value);

	void addFeatureAtIdx(double value, unsigned int idx);

	void addFeatures(const std::vector<double> &values);

	void addFeatures(const std::vector<int> &values);

	unsigned int getTotalLength() const { return _fv_idx; };

	feature_t getFeature(int idx) const { return _fv[idx]; };

	std::vector<feature_t>::const_iterator getFeatureBegin() const { return _fv.begin(); };

	std::vector<feature_t>::const_iterator getFeatureEnd() const { return _fv.end(); };

	unsigned int getNumSetFeatures() const { return _fv.size(); };

	void writeDebugInfo(std::ostream &out) const;

	bool equals(const FeatureVector &other_fv) const;
};
