#!/bin/bash
cd src
make
cp MRC ../ && cd ../
cd PgRC
make

cd FaStore
make
cp scripts/* ./