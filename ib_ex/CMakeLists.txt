ADD_EXECUTABLE(ib_ex ib_ex.cpp)
SET_TARGET_PROPERTIES(ib_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ib_ex" OUTPUT_NAME ib_ex)
TARGET_LINK_LIBRARIES(ib_ex PRIVATE MultiphaseLib)

