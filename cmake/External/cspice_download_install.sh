#!/usr/bin/env bash

DIRECTORY=$(cd `dirname $0` && pwd)
pushd `dirname $0`

BUILD_EXIT_CODE=0

wget http://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z  \
&& ls \
&& gzip -d cspice.tar.Z \
&& tar xfv cspice.tar \
&& rm cspice.tar \
&& mkdir -p $DIRECTORY/../../build/external/install/include \
&& mkdir -p $DIRECTORY/../../build/external/install/include/cspice \
&& mkdir -p $DIRECTORY/../../build/external/install/lib \
&& mkdir -p $DIRECTORY/../../build/external/install/lib/cspice \
&& sudo mv cspice/include/* $DIRECTORY/../../build/external/install/include/cspice/ \
&& sudo mv cspice/lib/*     $DIRECTORY/../../build/external/install/lib/cspice/ \
&& rm -Rf cspice

BUILD_EXIT_CODE=$?

popd

exit "$BUILD_EXIT_CODE"