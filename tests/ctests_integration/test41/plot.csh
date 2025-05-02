#/bin/csh

gunzip *.gz
rm -f p.out

root -b -q $G4INSTALL/tests/ctests_integration/test41/Plot.C >& p.out

gzip *.out *.log
#
