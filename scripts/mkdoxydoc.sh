#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

cd apidocs
wget https://github.com/HaoZeke/doxyYoda/releases/download/0.0.2/doxyYoda_0.0.2.tar.gz
tar -xf doxyYoda_0.0.2.tar.gz
rm -rf doxyYoda_0.0.2.tar.gz
cd ../
doxygen apidocs/Doxygen-featom.cfg
python3 -m http.server --bind 0.0.0.0 8000 -d html
