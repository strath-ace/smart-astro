#!/usr/bin/env bash

pushd .

sudo apt-get -qq update \
&& sudo apt-get install --assume-yes gfortran
BEFORE_INSTALL_EXIT_CODE=$?

popd

exit "$BEFORE_INSTALL_EXIT_CODE"