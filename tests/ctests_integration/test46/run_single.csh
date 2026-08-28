#/bin/csh

rm -f $1.out

echo 'Run test46 at ' $HOST >& $1.out
echo 'PHYSLIST      ' $PHYSLIST >>& $1.out
echo 'PRIMARYBEAM   ' $PRIMARYBEAM >>& $1.out
echo 'Input macro   ' $1 '.in' >>& $1.out 

$G4BWORK/test46 $G4INSTALL/tests/ctests_integration/test46/$1.in >>& $1.out
mv test46.root $1.root

#
