/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features.cpp
#
# Description: 	Code for computing features for fragmentations.
#
#				Assume that we have a config file that lists
the feature
#				vectors to compute (line separated text).
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.
# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include <DataStructs/SparseIntVect.h>
#include <memory>

#include "Feature.h"

// Helper functions for multiple features
void replaceUncommonWithX(std::string &symbol, bool use_full_symbol_set) {

	// Replace uncommon symbols with X
	std::unique_ptr<std::vector<std::string>> ok_symbols;
	if (use_full_symbol_set)
		ok_symbols = std::make_unique<std::vector<std::string>>(OKsymbols());
	else
		ok_symbols = std::make_unique<std::vector<std::string>>(OKSymbolsLess());

	for (auto &str : *ok_symbols) {
		if (symbol == str) { return; }
	}
	symbol = "X";
}

// assume last is always "X"
int getSymbolsIndex(const std::string &symbol, bool use_full_symbol_set) {
	int index = 0;
	std::unique_ptr<std::vector<std::string>> ok_symbols;
	if (use_full_symbol_set)
		ok_symbols = std::make_unique<std::vector<std::string>>(OKsymbols());
	else
		ok_symbols = std::make_unique<std::vector<std::string>>(OKSymbolsLess());

	for (auto &str : *ok_symbols) {
		if (symbol == str) { break; }
		index++;
	}
	return index;
}
