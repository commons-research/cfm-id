##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-predict/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################

set ( SRC_FILES  main.cpp )

add_executable ( cfm-predict ${SRC_FILES} )
target_link_libraries ( cfm-predict cfm-code ${REQUIRED_LIBS})

install ( TARGETS cfm-predict
          DESTINATION ${CFM_OUTPUT_DIR} )
