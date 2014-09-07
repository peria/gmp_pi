#!/bin/zsh

declare -a cmake_opts
cmake_opts+=("-GNinja")
cmake_opts+=("-DCMAKE_C_COMPILER=clang")
cmake_opts+=("-DCMAKE_CXX_COMPILER=clang++")
cmake_opts+=("-DCMAKE_BUILD_TYPE=Release")

declare -a cmake_cxx_flags
cmake_cxx_flags+=("-std=c++11")
cmake_cxx_flags+=("-O3")

cmake_opts+=("-DCMAKE_CXX_FLAGS=${cmake_cxx_flags[*]}")

out_dir=out/Release

cd "$(dirname "$0")/.."
mkdir -p out
mkdir -p ${out_dir}
cd ${out_dir}
cmake "${cmake_opts[@]}" "../../src"
ninja "$@"
