function(pi_target_link_common_libraries target)
  pi_target_link_libraries(drm_${target})
endfunction()

function(drm_target_link_libraries target)
  pi_target_link_common_libraries(drm_${target})
endfunction()

function(drm_add_executable exe)
  add_executable(drm_${exe} ${ARGN})
  set_target_properties(drm_${exe} PROPERTIES OUTPUT_NAME ${exe})
endfunction()

function(drm_add_test exe)
  add_test(check-drm_${exe} ${exe})
endfunction()

add_library(chudnovsky chudnovsky.cc chudnovsky.h)
