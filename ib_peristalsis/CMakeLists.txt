ADD_EXECUTABLE(ib_peristalsis ib_ex.cpp)
SET_TARGET_PROPERTIES(ib_peristalsis PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ib_peristalsis" OUTPUT_NAME ib_ex)
TARGET_SOURCES(ib_peristalsis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VCTwoFluidStaggeredStokesOperator.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../FullFACPreconditioner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../RBGS.f ${CMAKE_CURRENT_SOURCE_DIR}/../INSVCTwoFluidStaggeredHierarchyIntegrator.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../FRMThnForcing.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../IBMultiphaseHierarchyIntegrator.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../IBMultiphaseEulerianForceFunction.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../IBMultiphaseCrossLinks.cpp)
TARGET_INCLUDE_DIRECTORIES(ib_peristalsis PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../")
TARGET_LINK_LIBRARIES(ib_peristalsis PRIVATE IBAMR::IBAMR2d)

