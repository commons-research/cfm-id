set(BASE_HEADERS
	Util.h
	FunctionalGroups.h
	Config.h
	FeatureVector.h
	Feature.h
	FeatureCalculator.h
)
set(BASE_SOURCES
	Util.cpp
	FunctionalGroups.cpp
	Config.cpp
	FeatureVector.cpp
	Feature.cpp
	FeatureCalculator.cpp
)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Features FEATURES_SRC_DIR)

add_library(cfm-code STATIC ${BASE_HEADERS} ${BASE_SOURCES})  # Include both headers and sources
target_link_libraries(cfm-code ${REQUIRED_LIBS})
