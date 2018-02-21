#!/usr/bin/env bash

pushd `dirname $0`
mkdir -p build && \
cd build && \
cmake .. && \
make
BUILD_EXIT_CODE=$?

popd

exit "$BUILD_EXIT_CODE"