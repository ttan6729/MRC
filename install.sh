#!/bin/bash
git submodule update --init
cd src
make
cp MRC ../ && cd ../

mv PgRC pgrc
cd pgrc
mkdir build
cd build
#/projects/NGSData/RNACompression/cmake-3.17.2/bin/cmake ..
cmake ..
make PgRC
cp PgRC ../../
cd ../../

cd minicom 
sh install.sh

#cd FaStore
#make
#cp scripts/* ./
