/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootPairs.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "IonRootPairs.h"

void IonRootPairs::compute(FeatureVector &fv, const std::unique_ptr<RootedROMol> &ion,
                           const std::unique_ptr<RootedROMol> &nl) const {
	int ring_break;
	nl->mol->getProp("IsRingBreak", ring_break);
	std::vector<path_t> paths;
	computeRootPaths(paths, ion, 2, false);
	addRootPairFeatures(fv, paths, ring_break);
}
