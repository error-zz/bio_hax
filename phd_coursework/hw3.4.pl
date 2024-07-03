 #!/usr/bin/perl -w

use strict;
no warnings;

my %cpt;

# CPT ONE IMPORT
# x=0
$cpt{"x=0"}{"y1=1"} = 0.1;
$cpt{"x=0"}{"y1=0"} = 1 - $cpt{"x=0"}{"y1=1"};
$cpt{"x=0"}{"y2=1"} = 0.3;
$cpt{"x=0"}{"y2=0"} = 1 - $cpt{"x=0"}{"y2=1"};
$cpt{"x=0"}{"y3=1"} = 0.5;
$cpt{"x=0"}{"y3=0"} = 1 - $cpt{"x=0"}{"y3=1"};
# x=1
$cpt{"x=1"}{"y1=1"} = 0.8;
$cpt{"x=1"}{"y1=0"} = 1 - $cpt{"x=1"}{"y1=1"};
$cpt{"x=1"}{"y2=1"} = 0.6;
$cpt{"x=1"}{"y2=0"} = 1 - $cpt{"x=1"}{"y2=1"};
$cpt{"x=1"}{"y3=1"} = 0.4;
$cpt{"x=1"}{"y3=0"} = 1 - $cpt{"x=1"}{"y3=1"};

# RUN
print "Y1\tY2\tY3\tY\tP(Y|X=0)\tP(Y|X=1)\tP(Z1|Y)\tP(Z2|Y)\n";
	
	my $tmp1 = $cpt{"x=0"}{"y1=0"} * $cpt{"x=0"}{"y2=0"} * $cpt{"x=0"}{"y3=0"};
	my $tmp2 = $cpt{"x=1"}{"y1=0"} * $cpt{"x=1"}{"y2=0"} * $cpt{"x=1"}{"y3=0"};
	print "0\t0\t0\t1\t$tmp1\t$tmp2\t0.1\t0.9\n";

	my $tmp1 = $cpt{"x=0"}{"y1=1"} * $cpt{"x=0"}{"y2=0"} * $cpt{"x=0"}{"y3=0"};
	my $tmp2 = $cpt{"x=1"}{"y1=1"} * $cpt{"x=1"}{"y2=0"} * $cpt{"x=1"}{"y3=0"};
	print "1\t0\t0\t2\t$tmp1\t$tmp2\t0.2\t0.8\n";

	my $tmp1 = $cpt{"x=0"}{"y1=0"} * $cpt{"x=0"}{"y2=1"} * $cpt{"x=0"}{"y3=0"};
	my $tmp2 = $cpt{"x=1"}{"y1=0"} * $cpt{"x=1"}{"y2=1"} * $cpt{"x=1"}{"y3=0"};
	print "0\t1\t0\t3\t$tmp1\t$tmp2\t0.3\t0.7\n";

	my $tmp1 = $cpt{"x=0"}{"y1=0"} * $cpt{"x=0"}{"y2=0"} * $cpt{"x=0"}{"y3=1"};
	my $tmp2 = $cpt{"x=1"}{"y1=0"} * $cpt{"x=1"}{"y2=0"} * $cpt{"x=1"}{"y3=1"};
	print "0\t0\t1\t4\t$tmp1\t$tmp2\t0.4\t0.6\n";

	my $tmp1 = $cpt{"x=0"}{"y1=1"} * $cpt{"x=0"}{"y2=1"} * $cpt{"x=0"}{"y3=0"};
	my $tmp2 = $cpt{"x=1"}{"y1=1"} * $cpt{"x=1"}{"y2=1"} * $cpt{"x=1"}{"y3=0"};
	print "1\t1\t0\t5\t$tmp1\t$tmp2\t0.5\t0.5\n";

	my $tmp1 = $cpt{"x=0"}{"y1=1"} * $cpt{"x=0"}{"y2=0"} * $cpt{"x=0"}{"y3=1"};
	my $tmp2 = $cpt{"x=1"}{"y1=1"} * $cpt{"x=1"}{"y2=0"} * $cpt{"x=1"}{"y3=1"};
	print "1\t0\t1\t6\t$tmp1\t$tmp2\t0.6\t0.4\n";

	my $tmp1 = $cpt{"x=0"}{"y1=0"} * $cpt{"x=0"}{"y2=1"} * $cpt{"x=0"}{"y3=1"};
	my $tmp2 = $cpt{"x=1"}{"y1=0"} * $cpt{"x=1"}{"y2=1"} * $cpt{"x=1"}{"y3=1"};
	print "0\t1\t1\t7\t$tmp1\t$tmp2\t0.7\t0.3\n";

	my $tmp1 = $cpt{"x=0"}{"y1=1"} * $cpt{"x=0"}{"y2=1"} * $cpt{"x=0"}{"y3=1"};
	my $tmp2 = $cpt{"x=1"}{"y1=1"} * $cpt{"x=1"}{"y2=1"} * $cpt{"x=1"}{"y3=1"};
	print "1\t1\t1\t8\t$tmp1\t$tmp2\t0.8\t0.2\n";


