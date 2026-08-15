#/bin/csh

gunzip -q *.gz

root -b -q $G4INSTALL/tests/ctests_integration/test37/Plot.C >! p.out

gzip -q *.log *.out
#
