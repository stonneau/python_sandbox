#! /bin/bash


python -m compileall ./src/
mkdir -p $1/cwc
mkdir -p $1/cwc/tools
cp ./src/*.py $1/cwc
cp ./src/tools/*.py $1/cwc


exit 0



