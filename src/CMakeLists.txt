







add_executable(dootoy1 dootoy1.cpp)
add_executable(dootoy2 dootoy2.cpp)
add_executable(cp-test2 cp-test2.cpp)
add_executable(dootoycp dootoycp.cpp)
add_executable(dootoycp-lhcb dootoycp-lhcb.cpp)
add_executable(dootoycp_spline dootoycp_spline.cpp)
add_executable(dootoycp_different dootoycp_different.cpp)
add_executable(dootoycp_omega dootoycp_omega.cpp)
add_executable(dootoycp_float dootoycp_float.cpp)
add_executable(dootoycp_gen dootoycp_gen.cpp)
add_executable(dootoycp_dgamma dootoycp_dgamma.cpp)
add_executable(pulldist pulldist.cpp)


target_link_libraries(dootoy2 ${LIBRARIES})
target_link_libraries(dootoycp ${LIBRARIES})
target_link_libraries(dootoycp-lhcb ${LIBRARIES})
target_link_libraries(dootoycp_spline ${LIBRARIES})
target_link_libraries(dootoycp_different ${LIBRARIES})
target_link_libraries(dootoycp_omega ${LIBRARIES})
target_link_libraries(dootoycp_float ${LIBRARIES})
target_link_libraries(dootoycp_gen ${LIBRARIES})
target_link_libraries(dootoycp_dgamma ${LIBRARIES})
target_link_libraries(pulldist ${LIBRARIES})
