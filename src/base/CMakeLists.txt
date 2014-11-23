# --------------
# Function definition
function(base_add_library target)
  add_library(base_${target}_lib ${ARGN})
endfunction()

function(base_target_link_common_libraries target)
  pi_target_link_libraries(base_${target}_lib ${ARGN})
endfunction()

function(base_target_link_libraries target)
  base_target_link_common_libraries(${target} ${ARGN})
endfunction()

function(base_add_test exe)
  base_add_executable(${exe} ${exe}.cc)
  pi_target_link_libraries(base_${exe}_exe gtest gtest_main)
  add_test(base_${exe}_test ${exe})
endfunction()

# ----------
# Dependencies
base_add_library(time time.cc)
