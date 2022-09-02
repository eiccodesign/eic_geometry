cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

find_package( DD4hep REQUIRED COMPONENTS DDCore DDG4 )
#find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )

#-----------------------------------------------------------------------------------
set(PROJECT_NAME hadron_endcap)
dd4hep_add_plugin(${PROJECT_NAME} SOURCES src/*.cpp
  )
target_link_libraries(${PROJECT_NAME}
  PUBLIC DD4hep::DDCore  DD4hep::DDRec
  #ActsCore ActsPluginDD4hep
  )

#-----------------------------------------------------------------------------------
# Moving compact, scripts, and hepmc_generation to install share directory
install(DIRECTORY compact scripts hepmc_generation ip6
    DESTINATION share/${PROJECT_NAME}
    )
# Moving main XML file to install share directory
INSTALL(FILES endcapP_insert.xml
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/
    )
# Moving setup script to install directory
INSTALL(FILES scripts/insert_setup.sh
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )