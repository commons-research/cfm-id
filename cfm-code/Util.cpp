#include "Util.h"
#include "FunctionalGroups.h"
#include <GraphMol/AtomIterators.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <INCHI-API/inchi.h>
#include <memory>
#include <vector>

double getMassTol(double abs_tol, double ppm_tol, double mass) {
	double mass_tol = (mass / 1000000.0) * ppm_tol;
	if (mass_tol < abs_tol) mass_tol = abs_tol;
	return mass_tol;
}

double getMonoIsotopicMass(const romol_ptr_t &mol) {

	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	double mass              = 0.0;
	int natoms               = mol->getNumAtoms();
	for (int i = 0; i < natoms; i++) {
		RDKit::Atom *atom  = mol->getAtomWithIdx(i);
		std::string symbol = atom->getSymbol();
		mass += pt->getMostCommonIsotopeMass(symbol);
		mass += atom->getTotalNumHs() * pt->getMostCommonIsotopeMass("H");
	}

	// Adjust the mass by one electron according to the charge
	int charge = RDKit::MolOps::getFormalCharge(*mol.get());
	if (charge == 1) mass -= MASS_ELECTRON;
	if (charge == -1) mass += MASS_ELECTRON;

	return mass;
}

// Helper function to find an atom with the given label
std::unique_ptr<RDKit::Atom> getLabeledAtom(romol_ptr_t mol, const char *label) {
	RDKit::ROMol::AtomIterator ai;
	int root = 0;
	for (ai = mol.get()->beginAtoms(); ai != mol.get()->endAtoms(); ++ai) {
		(*ai)->getProp(label, root);
		if (root) break;
	}
	if (root)
		return std::unique_ptr<RDKit::Atom>((*ai)->copy());
	else
		return nullptr;
}

int moleculeHasSingleRadical(const std::unique_ptr<RDKit::ROMol> &romol) {
	int num_radicals = 0;
	for (const auto &atom : romol->atoms()) {
		int ionic_frag_q;
		atom->getProp("IonicFragmentCharge", ionic_frag_q);
		if (ionic_frag_q != 0) continue;
		num_radicals += atom->getNumRadicalElectrons();
	}
	return (num_radicals == 1);
}

int addIonicChargeLabels(const std::unique_ptr<RDKit::ROMol> &mol) {
	std::vector<int> mapping;
	int num_frags = RDKit::MolOps::getMolFrags(*mol, mapping);
	int num_ionic = 0;
	for (auto &atom : mol->atoms()) {
		atom->setProp("IonicFragmentCharge", 0);
		if (num_frags > 1 && atom->getDegree() == 0 && atom->getFormalCharge() != 0) {
			atom->setProp("IonicFragmentCharge", atom->getFormalCharge());
			num_ionic++;
		}
	}
	return num_ionic;
}

void alterNumHs(std::unique_ptr<RDKit::Atom> &atom, int H_diff) {
	int nHs = atom->getTotalNumHs();
	atom->setNoImplicit(true);
	atom->setNumExplicitHs(nHs + H_diff);
}

romol_ptr_t createMolPtr(const char *smiles_or_inchi) {
	std::unique_ptr<RDKit::RWMol> rwmol;
	if (std::string(smiles_or_inchi).substr(0, 6) == "InChI=") {
		RDKit::ExtraInchiReturnValues rv;
		rwmol = std::unique_ptr<RDKit::RWMol>(RDKit::InchiToMol(smiles_or_inchi, rv));
	} else {
		rwmol = std::unique_ptr<RDKit::RWMol>(RDKit::SmilesToMol(smiles_or_inchi));
	}

	// we cast the RWMol to a ROMol to avoid the need to change the function signature

	auto mol = std::unique_ptr<RDKit::ROMol>(rwmol.release());
	addIonicChargeLabels(mol);
	return std::move(mol);
}

void softmax(std::vector<double> &weights, std::vector<double> &probs) {
	probs.clear();
	double sum = 0.0;
	for (auto weight : weights) {
		double tmp = std::exp(weight);
		probs.push_back(tmp);
		sum += tmp;
	}
	for (auto &prob : probs) { prob /= sum; }
}

void labelNitroGroup(const std::unique_ptr<RDKit::ROMol> &mol) {
	// NOTE this is a context specific solution for nitro group single bond oxygen
	std::unique_ptr<RDKit::FragCatParams> fparams =
	    std::unique_ptr<RDKit::FragCatParams>(new RDKit::FragCatParams(PI_BOND_FGRPS_PICKLE));
	const RDKit::MOL_SPTR_VECT &fgrps = fparams->getFuncGroups();
	for (auto &fgrp : fgrps) {
		std::string fg_name;
		fgrp->getProp("_Name", fg_name);

		for (auto &atom : mol->atoms()) {
			atom->setProp(fg_name, 0);
			atom->setProp(fg_name + "Charge", 0);
		}
		// The format for each match is (queryAtomIdx, molAtomIdx)
		std::vector<RDKit::MatchVectType> fgp_matches;
		RDKit::SubstructMatch(*mol, *fgrp, fgp_matches);
		for (auto &fgp_match : fgp_matches) {
			for (auto &match : fgp_match) {
				mol->getAtomWithIdx(match.second)->setProp(fg_name, 1);
				mol->getAtomWithIdx(match.second)
				    ->setProp(fg_name + "Charge", mol->getAtomWithIdx(match.second)->getFormalCharge());
			}
		}
	}
}

int getValence(const std::unique_ptr<RDKit::Atom> &atom) {
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	// Fetch or compute the valence of the atom in the input molecule (we disallow valence changes for now)
	int valence              = -1;
	unsigned int num_val     = pt->getValenceList(atom->getSymbol()).size();
	int def_val              = pt->getDefaultValence(atom->getSymbol());

	// special case for nitrogroup
	int on_nitro_group;
	atom->getProp("NitroGroup", on_nitro_group);
	if (atom->getSymbol() == "O" && atom->getFormalCharge() == -1 && on_nitro_group) {
		valence = 1;
		return valence;
	}
	if (num_val == 1 && def_val != -1) {
		valence = def_val; // Hack to cover many cases - which can otherwise get complicated
	} else {
		// This seems to work in most cases....
		valence = atom->getExplicitValence() + atom->getImplicitValence() + atom->getNumRadicalElectrons();
		if (4 - pt->getNouterElecs(atom->getAtomicNum()) > 0) {
			valence += atom->getFormalCharge();
		} else {
			valence -= atom->getFormalCharge();
		}
	}
	return valence;
}
