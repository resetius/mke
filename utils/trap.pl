
sub norm($$$$) {
	my ($x, $y, $z, $r) = @_;
	my $rr = sqrt($x * $x + $y * $y);
	($x, $y, $z) = ($r * $x / $rr, $r * $y / $rr, $z);
	return ($x, $y, $z);
}

my $r1  = 1.0;
my $zz2 = 1 / sqrt(2.0);
my $r2  = sqrt(1.0 - $zz2 * $zz2);

my ($x0, $y0, $z0) = (0.0, 1.0, 0.0);
my ($x1, $y1, $z1) = (1.0, -1.0 / sqrt(3.0), 0.0);
my ($x2, $y2, $z2) = (-1.0, -1.0 / sqrt(3.0), 0.0);

($x1, $y1, $z1) = norm($x1, $y1, $z1, $r1);
($x2, $y2, $z2) = norm($x2, $y2, $z2, $r1);

print "{$x0, $y0, $z0},\n";
print "{$x1, $y1, $z1},\n";
print "{$x2, $y2, $z2},\n";

my ($x3, $y3, $z3) = (0.0, 1.0, $zz2);
my ($x4, $y4, $z4) = (1.0, -1.0 / sqrt(3.0), $zz2);
my ($x5, $y5, $z5) = (-1.0, -1.0 / sqrt(3.0), $zz2);

($x3, $y3, $z3) = norm($x3, $y3, $z3, $r2);
($x4, $y4, $z4) = norm($x4, $y4, $z4, $r2);
($x5, $y5, $z5) = norm($x5, $y5, $z5, $r2);

print "{$x3, $y3, $z3},\n";
print "{$x4, $y4, $z4},\n";
print "{$x5, $y5, $z5},\n";

my ($x6, $y6, $z6) = (($x2 + $x1) / 2.0, ($y2 + $y1) / 2.0, ($z2 + $z1) / 2.0);
my ($x7, $y7, $z7) = (($x0 + $x2) / 2.0, ($y0 + $y2) / 2.0, ($z0 + $z2) / 2.0);
my ($x8, $y8, $z8) = (($x1 + $x0) / 2.0, ($y1 + $y0) / 2.0, ($z1 + $z0) / 2.0);

print "{$x6, $y6, $z6},\n";
print "{$x7, $y7, $z7},\n";
print "{$x8, $y8, $z8},\n";
