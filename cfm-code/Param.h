/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param.h
#
# Description: 	Class for parameterization of fragmentation probabilities
#				within the bayesian network fragmentation trees.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#pragma once
#include "FeatureVector.h"

#include <boost/container/vector.hpp>
#include <memory>
#include <string>

// Exception to throw when the input feature vector configuration doesn't match the parameters
class ParamFeatureMismatchException : public std::exception {

	const char *what() const noexcept override { return "Mismatch between feature vector length and num parameters"; }
};

class Param {
protected:
	std::vector<float> weights;

	unsigned int num_energy_levels;
	std::vector<std::string> feature_list;
	int expected_num_input_features;

	// Initialisation options
	virtual void randomUniformInit();
	virtual void randomNormalInit();
	void zeroInit();
	void fullZeroInit();

public:
	// Constructor to initialise parameter weight size from a feature list
	Param(const std::vector<std::string> &a_feature_list, int a_num_energy_levels);

	// Constructor to create skeleton parameter copy
	//(doesn't actually copy weights, just resizes to same)
	Param(Param &param_instance);

	// Constructor for loading parameters from file
	explicit Param(std::string &filename);

	// Append a set of parameters to the current parameters
	// Either all energy levels to the next slot (if energy < 0),
	// or just the energy level specified.
	void appendNextEnergyParams(Param &next_param, int energy = -1);

	// Copy the highest energy parameters to create another higher energy level
	void appendRepeatedPrevEnergyParams();

	// Compute the theta value for an input feature vector and energy based
	// on the current weight settings
	virtual float computeTheta(const FeatureVector &fv, int energy);

	// Set the value of a weight
	void setWeightAtIdx(float value, int index) { weights[index] = value; };

	// Save parameters to file
	virtual void saveToFile(std::string &filename);

	// Access functions
	[[nodiscard]] float getWeightAtIdx(int index) const { return weights[index]; };

	[[nodiscard]] const std::vector<float> &getWeights() const { return weights; };
	[[nodiscard]] std::vector<float> &getWeights() { return weights; }

	// this will be changed once we add dropouts for linear model, for now
	// it only return nullptr
	//  NOTE use boost vector of bool for mpi
	//  because std vector can not return bool&
	virtual std::unique_ptr<boost::container::vector<bool>> getDropoutsPtr() { return nullptr; };

	virtual std::unique_ptr<std::vector<float>> getDropoutsProbPtr() { return nullptr; };

	[[nodiscard]] unsigned int getNumWeights() const { return weights.size(); };

	[[nodiscard]] unsigned int getNumWeightsPerEnergyLevel() const { return weights.size() / num_energy_levels; };

	[[nodiscard]] unsigned int getNumEnergyLevels() const { return num_energy_levels; };

	[[nodiscard]] const std::vector<std::string> &getFeatureNames() const { return feature_list; };

	virtual void getBiasIndexes(std::vector<unsigned int> &bias_indexes) { bias_indexes.push_back(0); };

	virtual void initWeights(int init_type);

	virtual void rollDropouts() {};

	void readFromFile(const std::string &filename);

	// Function To set weights from a vector
	// Used for Unit tests
	virtual void setWeights(const std::vector<float> &values) {
		weights.clear();
		weights = values;
	}
};