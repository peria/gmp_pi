# --------------
# Function definition
function(drm_add_library target)
  add_library(drm_${target}_lib ${ARGN})
endfunction()

function(drm_target_link_common_libraries target)
  pi_target_link_libraries(drm_${target}_lib)
endfunction()

function(drm_target_link_libraries target)
  # add other common libraries defined in drm/
  drm_target_link_common_libraries(${target})
endfunction()

function(drm_add_executable exe)
  add_executable(drm_${exe}_exe ${ARGN})
  set_target_properties(drm_${exe}_exe PROPERTIES OUTPUT_NAME drm_${exe})
endfunction()

function(drm_add_test exe)
  drm_add_executable(${exe} ${exe}.cc)
  pi_target_link_libraries(drm_${exe}_exe gtest gtest_main)
  add_test(drm_${exe}_test ${exe})
endfunction()

# ----------
# Dependencies

drm_add_library(chudnovsky chudnovsky.cc)
drm_target_link_libraries(chudnovsky)

drm_add_executable(chudnovsky chudnovsky_main.cc)
target_link_libraries(drm_chudnovsky_exe drm_chudnovsky_lib)

