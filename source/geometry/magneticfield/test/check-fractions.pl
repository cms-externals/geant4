#!/usr/bin/perl

# eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;
			# process any FOO=bar switches

$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

printf "# %10s  %12s  %7s  %7s\n", 'Move', 'Step', 'Diff', 'fraction';
$maxfraction = -1.000;

while (<>) {
    chomp;	# strip record separator
    if( m/Move displaced/) { 
        s/Move displaced//; s/\- further than step=//; s/by//; s/, *a fraction of//; 
       @Fld = split(' ', $_, -1);
       # if (/^# /) { print $_; printf "\n"; }
       if (!/^# /) {
	   $move = $Fld[(1)-1];
	   $step = $Fld[(2)-1];
	   $diff = $Fld[(3)-1];
	   $fraction = $Fld[(4)-1];
	   #  Diff is always positive - it exceeds step
	   # if ( fraction > 3.0e-6) { printf "%12.7f  %12.7f  %7.3e  %7.3e\n", move, step, diff, fraction;; }
	   if ($fraction > $maxfraction) {	#???
	       printf "%12.7f  %12.7f  %7.3e  %7.3e\n", $move, $step, $diff, $fraction;
	       $maxfraction = $fraction 
	   }
       } 
    } 
}
