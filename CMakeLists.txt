install_fhicl()
install_source()
add_subdirectory(GENIE)
add_subdirectory(HepMC)
add_subdirectory(HeavySterileNeutrino)
add_subdirectory(HiggsPortalScalarGen)
add_subdirectory(HNL)
add_subdirectory(HeavyNeutralLeptonGen)
art_make( BASENAME_ONLY
          MODULE_LIBRARIES 	  		   
		   larcorealg_Geometry
		   ubcore_Geometry
		   larcore_Geometry_Geometry_service
		   nutools_ParticleNavigation lardataobj_Simulation
		   nusimdata_SimulationBase
		   lardata_Utilities
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   nusimdata_SimulationBase
		   nutools_RandomUtils_NuRandomService_service
		   ${LARDATA_LIB}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                   ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		   art_Persistency_Common 
		   art_Persistency_Provenance 
		   art_Utilities 
		   canvas
		   ${MF_MESSAGELOGGER}
		   
		   ${FHICLCPP}
		   ${CLHEP}
		   ${CETLIB}
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}
		   ${ROOT_EG}
		)

