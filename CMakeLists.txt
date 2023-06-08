cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

PROJECT(hadron_endcap
  LANGUAGES CXX
  )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

# Export compile commands as json for run-clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(DD4hep_required_components DDCore DDRec DDG4)

# The items passed to target_link_libraries below need DD4hep:: prepended
set(DD4hep_required_libraries ${DD4hep_required_components})
list(TRANSFORM DD4hep_required_libraries PREPEND DD4hep::)

find_package(DD4hep REQUIRED COMPONENTS ${DD4hep_required_components})

find_package(ActsDD4hep)
if(ActsDD4hep_FOUND)
  add_compile_definitions(USE_ACTSDD4HEP)
  set(ActsDD4hep ActsDD4hep::ActsDD4hep)
else()
  find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep)
  set(ActsDD4hep ActsCore ActsPluginDD4hep)
endif()
find_package(fmt REQUIRED)


#-----------------------------------------------------------------------------------
dd4hep_configure_output()

file(GLOB sources CONFIGURE_DEPENDS src/*.cpp)

dd4hep_add_plugin(${PROJECT_NAME}
  SOURCES ${sources}
  USES ${ActsDD4hep} ROOT::Core ROOT::Gdml
  )
target_link_libraries(${PROJECT_NAME}
  PUBLIC ${DD4hep_required_libraries}
  )

#-----------------------------------------------------------------------------------
# Moving compact, scripts, and hepmc_generation to install share directory
install(DIRECTORY compact scripts hepmc_generation
    DESTINATION share/${PROJECT_NAME}
    FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
    )
# Moving main XML file to install share directory
INSTALL(FILES hadron_endcap.xml
    DESTINATION share/${PROJECT_NAME}
    )
