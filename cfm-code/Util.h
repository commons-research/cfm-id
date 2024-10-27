#pragma once

#include <GraphMol/Atom.h>
#include <GraphMol/ROMol.h>
#include <memory>
#include <random>

typedef std::unique_ptr<RDKit::ROMol> romol_ptr_t;
static const double MASS_ELECTRON = 0.000548579909;

struct beliefs_t {                       // IN LOG DOMAIN!
	std::vector<std::vector<double>> tn; // Transition:  Indexed by transition, then depth
	std::vector<std::vector<double>> ps; // Persistence: Indexed by fragment, then depth
};

// Helper Function to figure out the mass tolerance to use based on
// absolute and ppm tolerances for the current mass
double getMassTol(double abs_tol, double ppm_tol, double mass);

// Structure for storing a shared pointer to a molecule, as well as atom pointers
// to the root atoms (i.e. the atoms that were at either end of a broken bond)
class RootedROMol {
public:
	RootedROMol() = default; // Default Constructor
	RootedROMol(romol_ptr_t a_mol, std::unique_ptr<RDKit::Atom> a_root)
	    : mol(std::move(a_mol)), root(std::move(a_root)) {}
	romol_ptr_t mol;
	std::unique_ptr<RDKit::Atom> root = std::make_unique<RDKit::Atom>();
};
;

// Helper function to compute the monoisotopic mass of a molecule
double getMonoIsotopicMass(const romol_ptr_t &mol);

// Helper function to find an atom with the given label
std::unique_ptr<RDKit::Atom> getLabeledAtom(romol_ptr_t mol, const char *label);

// Helper function to check for radical electrons within a molecule
int moleculeHasSingleRadical(const std::unique_ptr<RDKit::ROMol> &mol);

// Helper function to identify and label ionic charges
int addIonicChargeLabels(const std::unique_ptr<RDKit::ROMol> &mol);

// Helper function to alter number of Hs on an atom (accounting for implicit Hs)
void alterNumHs(std::unique_ptr<RDKit::Atom> &atom, int H_diff);

romol_ptr_t createMolPtr(const char *smiles_or_inchi);

// Helper function label NirtoGroup
void labelNitroGroup(const std::unique_ptr<RDKit::ROMol> &mol);

// Helper function to get valence
int getValence(const std::unique_ptr<RDKit::Atom> &atom);

inline double logAdd(double log_x, double log_y) {

	// ensure log_y >= log_x, can save some expensive log/exp calls
	if (log_x > log_y) {
		double t = log_x;
		log_x    = log_y;
		log_y    = t;
	}

	double rval = log_y - log_x;

	// only replace log(1+exp(log(y)-log(x))) with log(y)-log(x)
	// if the the difference is small enough to be meaningful
	if (rval > 100.0) return log_y;
	rval = std::log(1.0 + std::exp(rval));
	rval += log_x;
	return rval;
}

void softmax(std::vector<double> &weight, std::vector<double> &prob);

// Init static members
static std::random_device util_rd;
static std::mt19937 util_rng(util_rd());

// A Big DBL number
static const double A_BIG_DBL               = 1000000000.0;
static const double WEIGHT_SELECTION_SCALER = 10000.0;
