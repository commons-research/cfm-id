/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FeatureCalculator.h
#
# Description: 	Code for computing features for fragmentations.
#
# Copyright (c) 2018
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#pragma once

#include "Feature.h"
#include "FeatureVector.h"
#include "Util.h"

#include <GraphMol/FragCatalog/FragCatParams.h>

#include <boost/lexical_cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Exception to throw when the input feature configuration file is invalid
class InvalidConfigException : public std::exception {

	[[nodiscard]] const char *what() const noexcept override { return "Invalid Feature Configuration File"; }
};

class FeatureCalculationException : public std::exception {
private:
	std::string _message;

public:
	explicit FeatureCalculationException(std::string message) noexcept : _message(std::move(message)) {};

	[[nodiscard]] const char *what() const noexcept override {
		std::cout << "Error computing feature vector: " << _message << std::endl;
		return _message.c_str();
	}

	~FeatureCalculationException() noexcept override = default;
};

// Class to compute a feature vector
class FeatureCalculator {
private:
	// Indexes of feature classes that are selected for use
	std::vector<int> _used_break_feature_idxs;
	std::vector<int> _used_fragement_feature_idxs;
	// List of feature classes ready to be used
	static const boost::ptr_vector<BreakFeature> &breakFeatureCogs();

	// List of fragmentation features ready to be used
	static const boost::ptr_vector<FragmentFeature> &fragmentFeatureCogs();

	// Helper function - Configure feature for use
	void configureFeature(std::string &name);

public:
	// Constructor: Initialise the calculator using a config file listing features
	explicit FeatureCalculator(std::string &config_filename);

	// Constructor: Initialise the calculator using a list of feature names
	explicit FeatureCalculator(std::vector<std::string> &feature_list);

	// Compute the expected number of total features
	unsigned int getNumFeatures();

	// Retrieve the list of feature names being used
	std::vector<std::string> getFeatureNames();

	// Retrieve a list of valid feature names (for testing)
	static const std::vector<std::string> getValidFeatureNames();

	// Compute the feature vector for the input ion and nl (with labeled Root
	// atoms)
	// - NB: responsibility of caller to delete.
	FeatureVector *computeFeatureVector(const RootedROMol *ion, const RootedROMol *nl,
	                                    const romol_ptr_t &precursor_ion);

	bool includesFeature(const std::string &fname);
};
