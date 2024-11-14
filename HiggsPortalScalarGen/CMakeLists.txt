find_ups_product( ifdhc )
find_ups_product( ifdh_art )

# ifdhc
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
# ifdhc_art  # as a service
cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB NO_DEFAULT_PATH )
include_directories( $ENV{IFDH_ART_FQ_DIR}/inc )

#include_directories(${IFDH_ART_INC})
include_directories(${IFDHC_FQ_DIR}/inc)

find_ups_product(genie)

art_make(BASENAME_ONLY
          LIBRARY_NAME uboonecode_HiggsPortalScalarGen
          LIB_LIBRARIES
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nutools_RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            art_Persistency_Common canvas
            art_Persistency_Provenance canvas
            art_Utilities canvas
            cetlib cetlib_except
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CLHEP}
            ${ROOT_BASIC_LIB_LIST}
            ${SQLITE3}
            ${ROOT_EG}
            ${ROOT_GEOM}
            ${IFDH}
            ${IFDH_SERVICE}
            ${Boost_SYSTEM_LIBRARY}
            ${DK2NU_TREE}
            ${GFWMSG}
            ${GFWPARDAT}
            ${LOG4CPP}
      MODULE_LIBRARIES
            uboonecode_HiggsPortalScalarGen
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nutools_RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            art_Persistency_Common canvas
            art_Persistency_Provenance canvas
            art_Utilities canvas
            cetlib cetlib_except
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CLHEP}
            ${ROOT_BASIC_LIB_LIST}
            ${SQLITE3}
            ${ROOT_EG}
            ${ROOT_GEOM}
            ${IFDH}
            ${IFDH_SERVICE}
            ${DK2NU_TREE}
            ${GFWMSG}
            ${GFWPARDAT}
            ${LOG4CPP}
            )

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

