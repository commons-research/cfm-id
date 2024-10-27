#pragma once
#include "Spectrum.h"
#include "Util.h"
#include <map>
#include <string>

class EmptyIsotopeSpectrumException : public std::exception {

	const char *what() const noexcept override {
		return "Error in isotope spectrum generation: - no peaks generated above threshold.";
	}
};

struct ipeak {
	double mass;
	double rel_area;
};

// map from element abbreviation to index in the elements table
typedef std::map<std::string, size_t> ElemMap;

// map from element index to the count of occurences in the formula
typedef std::map<size_t, long> FormMap;

typedef std::vector<ipeak> Pattern;               // index: peak_number
typedef std::vector<Pattern> SuperAtomList;       // index: bit_number
typedef std::vector<SuperAtomList> SuperAtomData; // index: element_number

// Exception to throw when something goes wrong during the isotope calculation
class IsotopeCalculationException : public std::exception {

	const char *what() const noexcept override { return "Exception occurred computing isotope peaks"; }
};

class IsotopeCalculator {

public:
	IsotopeCalculator(double a_intensity_thresh, std::string &isotope_pattern_file)
	    : verbose(false), intensity_thresh(a_intensity_thresh) {
		init_data(isotope_pattern_file.c_str());
	};

	void computeIsotopeSpectrum(Spectrum &output, const romol_ptr_t &mol, long charge);

	void setVerbose() { verbose = true; };

	double getIntensityThresh() const { return intensity_thresh; };

private:
	double intensity_thresh;
	SuperAtomData sad; // atom_idx -> isotope information
	ElemMap em;        // Atom symbol -> atom_idx

	// This function sets the formula map (atom_idx -> count) structure for the input molecule
	void setFormulaMap(FormMap &output, const romol_ptr_t &mol);

	// The remainder of the functions are copied directly from emass (with minor mods)
	void init_data(const char *filename);

	void convolute_basic(Pattern &h, const Pattern &g, const Pattern &f);

	void prune(Pattern &f, double limit);

	void calculate(Pattern &tmp, Pattern &result, FormMap &fm, double limit, long charge);

	void print_pattern(Pattern &result, int digits);

	// This function writes the output to our spectrum output
	void print_to_output(Spectrum &output, Pattern &result);

	bool verbose;
};
