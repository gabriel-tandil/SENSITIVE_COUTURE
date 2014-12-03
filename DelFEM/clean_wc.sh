#!/bin/sh
cd `dirname $0`
find . -name *.1 -maxdepth 5 -type f -exec rm -f {} \;
find . -name *.o -maxdepth 5 -type f -exec rm -f {} \;
find . -name Thumbs.db -type f -maxdepth 5 -exec rm -f {} \;
find . -name .DS_Store -type f -maxdepth 5 -exec rm -f {} \;
find . -name *.pbxuser -type f -maxdepth 5 -exec rm -f {} \;
find . -name *.perspectivev3 -type f -maxdepth 5 -exec rm -f {} \;
find . -name *.dsp -maxdepth 5 -type f -exec rm -f {} \;
find . -name *.dsw -maxdepth 5 -type f -exec rm -f {} \;
find . -name .svn -type d -maxdepth 5 -exec rm -rf {} \;
find . -name build -type d -maxdepth 5 -exec rm -rf {} \;
find . -name "\\moc_*.cpp" -type f -maxdepth 5 -exec rm -f {} \;