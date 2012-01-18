#!/bin/bash
export ENABLE_ITAPS=1
wget http://www.scorec.rpi.edu/FMDB/source/buildFMDBSerial.env
wget http://www.scorec.rpi.edu/FMDB/source/buildFMDBSerial.sh
wget http://www.scorec.rpi.edu/FMDB/source/auxillaryBuildScripts.tar
tar xvf auxillaryBuildScripts.tar
chmod 755 buildFMDBSerial.sh
source buildFMDBSerial.env 
./buildFMDBSerial.sh

