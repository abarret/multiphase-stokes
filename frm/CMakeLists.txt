ADD_EXECUTABLE(frm four_roll_mill.cpp)
TARGET_LINK_LIBRARIES(frm MultiphaseLib)

ADD_EXECUTABLE(frm_convergence convergence_tester.cpp)
TARGET_LINK_LIBRARIES(frm_convergence MultiphaseLib)

