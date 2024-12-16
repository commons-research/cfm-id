/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGenerator.h
#
# Description: 	FragmentTree class for holding the results of a generated
#				fragment tree.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MILP_H__
#define __MILP_H__

#include <GraphMol/ROMol.h>
#include <vector>

class MILP {
private:
	RDKit::ROMol *_mol;
	int _fragmentidx;
	int _broken_ringidx; // Store the idx of any broken rings (or -1 if there are none).
	bool _verbose;
	int _status = 0;

	// Helper functions:
	// Allows traversal of a ring one bond at a time
	RDKit::Bond *getNextBondInRing(RDKit::Bond *bond, RDKit::Atom *atom, std::vector<int> &ring_bond_flags);

	// Checks whether an atom should be allowed a lone pair bond (not including those already using theirs to create an
	// extra single bond)
	static int getAtomLPLimit(RDKit::Atom *atom);

	void printConstraint(int num_terms, int *colno, bool ge, int val) const;

public:
	MILP(RDKit::ROMol *a_mol, int a_fragmentidx, int a_broken_ringidx, bool a_verbose)
	    : _mol(a_mol), _fragmentidx(a_fragmentidx), _broken_ringidx(a_broken_ringidx), _verbose(a_verbose) {};

	MILP(RDKit::ROMol *a_mol, int a_fragmentidx, bool a_verbose)
	    : _mol(a_mol), _fragmentidx(a_fragmentidx), _broken_ringidx(-1), _verbose(a_verbose) {};

	int runSolver(std::vector<int> &output_bmax, bool allow_lp_q, int max_free_pairs, bool allow_rearrangement);
};

#endif // __MILP_H__
