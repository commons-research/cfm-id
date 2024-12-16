set(BASE_HEADERS
	Util.h
	FunctionalGroups.h
	Config.h
	FeatureVector.h
	Feature.h
	FeatureCalculator.h
	Isotope.h
	Spectrum.h
	Version.h
	Solver.h
	Param.h
	NNParam.h
	MspReader.h
	Message.h
	MILP.h
	FragmentTreeNode.h
)
set(BASE_SOURCES
	Util.cpp
	FunctionalGroups.cpp
	Config.cpp
	FeatureVector.cpp
	Feature.cpp
	FeatureCalculator.cpp
	Isotope.cpp
	Spectrum.cpp
	Solver.cpp
	Param.cpp
	NNParam.cpp
	MspReader.cpp
	Message.cpp
	MILP.cpp
	FragmentTreeNode.cpp
)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Features FEATURES_SRC_DIR)

add_library(cfm-code STATIC ${BASE_HEADERS} ${BASE_SOURCES})  # Include both headers and sources
target_link_libraries(cfm-code ${REQUIRED_LIBS})
