#/bin/csh

gunzip *.gz

ln -s $G4INSTALL/tests/ctests_integration/test44/Exp_Data/*.txt ./

foreach tPart (p he4 c12)
    setenv PARTICLE $tPart
    root -b -q ${G4INSTALL}/tests/ctests_integration/test44/Plot.C
end

gzip *.log *.out *.root
