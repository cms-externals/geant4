#/bin/csh

if ( $?REFERENCE == 0 ) then
setenv REFERENCE `date '+%m_%d_%Y-%H:%M:%S'`
endif

cd $VFEM/test44

mkdir -p $REFERENCE
cd $REFERENCE
rm -f *.txt *.gz

set    work = "$G4BWORK/test44"
set    dir  = "$G4INSTALL/tests/ctests_integration/test44/"

ln -s ${dir}/Exp_Data/*.txt ./

setenv PHYSLIST  QBBC

set tPart = (p he4 c12)
foreach phys (opt0 opt3 opt4) 
    foreach part ($tPart)
	set file = "${tPart}_water_${phys}.log"
	if ( -e "$file" )  then
	    rm -f $file
	endif
	${work} ${dir}${part}_water_${phys}.in >& ${part}_water_${phys}.log
	mv test44.root ${part}_${phys}.root
    end
end

gzip *.log *.out
