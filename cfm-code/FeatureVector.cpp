//*#########################################################################
// # Mass Spec Prediction and Identification of Metabolites
// #
// # FeatureVector.cpp
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

#include "FeatureVector.h"

#include <iostream>

void FeatureVector::addFeature(double value) {
	if (value != 0.0)
		_fv.push_back(_fv_idx++);
	else
		_fv_idx++;
}

void FeatureVector::addFeatureAtIdx(double value, unsigned int idx) {
	if (_fv_idx <= idx) _fv_idx = idx + 1;
	if (value != 0.0) _fv.push_back(idx);
}

void FeatureVector::addFeatures(const std::vector<double> &values) {
	for (const auto &value : values) { this->addFeature(value); }
}

void FeatureVector::addFeatures(const std::vector<int> &values) {
	for (const auto &value : values) { this->addFeature((double)value); }
}

bool FeatureVector::equals(const FeatureVector &other_fv) const {
	auto other_fv_fv = other_fv._fv;
	bool equal       = _fv.size() == other_fv_fv.size() && std::equal(_fv.begin(), _fv.end(), other_fv_fv.begin()) &&
	             _fv_idx == other_fv._fv_idx;
	return equal;
}

// print debug info
void FeatureVector::writeDebugInfo(std::ostream &out) const {
	// out<< "fv_vector_size: " << fv.size() << std::endl;
	out << "fv_idx : " << _fv_idx << " fv none zero idx: ";
	for (int i = 0; i < _fv.size(); i++) out << _fv[i] << " ";
}
