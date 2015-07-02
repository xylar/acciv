
PROJECT(acciv)


SET(INCLUDE_DIRS
    ./Utility
    ${HDF5_INCLUDE_DIR}
    ${FFTW_INCLUDE_DIR}
)


SET(SRC_LIB 
    ./Utility/tools/URandom.cpp
    ./Utility/core/UString.cpp
    ./Utility/platform/UStreamFile.cpp
)

SET(SRC 
    ./ScalarField2D.cpp
    ./VectorField2D.cpp
    ./QuadTree.cpp
    ./Correlator.cpp
    ./ParameterFileReader.cpp
    ./TiePointSet.cpp
    ./ScatteredVectorData2D.cpp
    ./ACCIVPass.cpp
    ./MaskedImage.cpp
    ./ScatteredScalarData2D.cpp
)

SET(INC 
    ./Utility/tools/URandom.h
    ./Utility/containers/UList.h
    ./Utility/containers/UArray.h
    ./Utility/containers/UAllocationPool.h
    ./Utility/containers/UMap.h
    ./Utility/core/Utility.h
    ./Utility/core/UConstants.h
    ./Utility/core/UTypes.h
    ./Utility/core/UString.h
    ./Utility/platform/UStreamFile.h
    ./vector2.h
    ./TiePointSet.h
    ./ScalarField2D.h
    ./PowelMinimizer.h
    ./MakeDirectory.h
    ./ScatteredVectorData2D.h
    ./VectorField2D.h
    ./LineMethodMinimizer.h
    ./Correlator.h
    ./QuadTree.h
    ./BrentMinimizer.h
    ./SquareBoundedMinimizer.h
    ./ParticleIntegrator.h
    ./ACCIVPass.h
    ./MaskedImage.h
    ./BulirschStoerIntegrator.h
    ./ScatteredScalarData2D.h
    ./ScatteredData2D.h
    ./GridData2D.h
    ./ParameterFileReader.h
)


include_directories(${INCLUDE_DIRS})

SET(LINK_LIBS
    ${HDF5_LIBRARIES}
    ${FFTW_LIBRARIES}
)

ADD_LIBRARY(acciv-utils STATIC ${SRC_LIB} ${INCLUDE_FILES})



include_directories(${INCLUDE_DIRS})
SET(SRC_ACCIV ./ACCIV.cpp ${SRC})
ADD_EXECUTABLE(acciv ${SRC_ACCIV} ${INC})
target_link_libraries(acciv "acciv-utils" ${LINK_LIBS} )


include_directories(${INCLUDE_DIRS})
SET(SRC_SMOOTHFIT ./smoothFit.cpp ${SRC})
ADD_EXECUTABLE(acciv-smoothFit ${SRC_SMOOTHFIT} ${INC})
target_link_libraries(acciv-smoothFit "acciv-utils" ${LINK_LIBS} )


include_directories(${INCLUDE_DIRS})
SET(SRC_MAKEMOVIE ./makeMovie.cpp ${SRC})
ADD_EXECUTABLE(acciv-makeMovie ${SRC_MAKEMOVIE} ${INC})
target_link_libraries(acciv-makeMovie "acciv-utils" ${LINK_LIBS} )


include_directories(${INCLUDE_DIRS})
SET(SRC_MAKEGEOMFACTOR ./makeGeometryFactors.cpp ${SRC})
ADD_EXECUTABLE(acciv-makeGeometryFactors ${SRC_MAKEGEOMFACTOR} ${INC})
target_link_libraries(acciv-makeGeometryFactors "acciv-utils" ${LINK_LIBS} )


install(TARGETS acciv acciv-smoothFit acciv-makeMovie acciv-makeGeometryFactors RUNTIME DESTINATION bin)
