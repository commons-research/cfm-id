/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Spectrum.h
#
# Description: 	Class for spectrum.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

#include <ostream>
#include <set>
#include <vector>

typedef std::pair<int, double> annotation_t; //<Fragment Id, Score>

static bool sort_annotations_by_score(const annotation_t &u, const annotation_t &v) { return u.second > v.second; }

class Peak {
private:
	double _mass      = 0.0;
	double _intensity = 0.0;
	std::vector<annotation_t> _annotations;

public:
	Peak() = default;

	Peak(double a_mass, double an_intensity) : _mass(a_mass), _intensity(an_intensity) {};

	double mass() const { return _mass; };
	double intensity() const { return _intensity; };
	[[nodiscard]] const std::vector<annotation_t> &annotations() const { return _annotations; };
	std::vector<annotation_t> &annotations() { return _annotations; };
	void new_mass(double a_mass) { _mass = a_mass; };
	void new_intensity(double an_intensity) { _intensity = an_intensity; };
};

static bool sort_peaks_by_intensity(const Peak &u, const Peak &v) { return u.intensity() > v.intensity(); }

static bool sort_peaks_by_mass(const Peak &u, const Peak &v) { return u.mass() < v.mass(); }

class Spectrum {
private:
	std::vector<Peak> _peaks;
	bool _is_normalized;
	bool _is_sorted;

	double getMaxIntensity() const;

public:
	Spectrum() : _is_normalized(false), _is_sorted(false) {};

	// copy constructor
	Spectrum(const Spectrum &spectrum) {
		for (const auto &peak : spectrum) push_back(peak);

		// overwrite the those flags
		// do not move it up, those flags are set to false by pusb_back function
		_is_normalized = spectrum.is_normalized();
		_is_sorted     = spectrum.is_sorted();
	};

	bool is_normalized() const { return _is_normalized; };
	bool is_sorted() const { return _is_sorted; };

	// Iterating over the spectrum = iterating over the peaks
	typedef std::vector<Peak>::const_iterator const_iterator;

	const_iterator begin() const { return _peaks.begin(); };

	const_iterator end() const { return _peaks.end(); };

	void push_back(const Peak &pk) {
		_peaks.push_back(pk);
		_is_normalized = false;
		_is_sorted     = false;
	};

	void clear() {
		_peaks.clear();
		_peaks.resize(0);
		_is_normalized = false;
		_is_sorted     = false;
	};

	unsigned int size() const { return _peaks.size(); };

	const Peak *getPeak(int peak_idx) const { return &_peaks[peak_idx]; };

	void addPeakAnnotation(int peak_idx, annotation_t annot) { _peaks[peak_idx].annotations().push_back(annot); };

	void updateAnnotationId(int peak_idx, int annot_idx, int new_id) {
		_peaks[peak_idx].annotations()[annot_idx].first = new_id;
	}

	const std::vector<Peak> &getPeaks() const { return _peaks; };

	bool isNormalizedAndSorted() const { return _is_normalized && _is_sorted; };

	void postProcess(double perc_thresh, int min_peaks, int max_peaks, double min_relative_intensity_prec = 0.0);

	void normalizeAndSort();

	void clean(double abs_mass_tol, double ppm_mass_tol);

	void sortAndNormalizeAnnotations();

	void outputToStream(std::ostream &out, bool do_annotate, int mz_precision, bool normalize_to_max = true) const;

	void outputToMspStream(std::ostream &out, std::string id, int ionization_mode, int energy,
	                       std::string &smiles_or_inchi, int mz_precision) const;

	void outputToMgfStream(std::ostream &out, std::string id, int ionization_mode, int energy, double mw,
	                       std::string &smiles_or_inchi, int mz_precision) const;

	int removePeaksWithNoFragment(std::vector<double> &frag_masses, double abs_tol, double ppm_tol);

	void quantisePeaksByMass(int num_dec_places);

	void convertToLogScale();

	void convertToLinearScale();

	void getDisplayedFragmentIds(std::set<int> &ids, bool normalize_to_max = true) const;
};

#endif // __SPECTRUM_H__
