#!/usr/bin/perl

use strict;

use DBI;
use Digest::MD5  qw(md5 md5_hex md5_base64);
use POSIX;

#print "drivers:\n";
#my @drivers = DBI->available_drivers;
#for my $i (@drivers) {
#	print STDERR "$i\n";
#}
#print "---\n";

my $db_str   = "dbi:Pg:dbname=science";
my $username = "postgres";
my $password = "";

my $dbh = DBI->connect($db_str, $username, $password, {PrintError => 1, AutoCommit => 0});
if ($DBI::err != 0) {
	print $DBI::errstr . "\n";
	exit($DBI::err);
}

my $id = "CREATE SEQUENCE experiment_id_seq";
my $table1 = "CREATE TABLE experiment (
	id INT UNIQUE,
	name VARCHAR(256), -- link to experiment_table
	t TIMESTAMP
)";
my $alter = "ALTER TABLE experiment ALTER COLUMN id SET DEFAULT NEXTVAL('experiment_id_seq')";

my $table2 = "CREATE TABLE barvortex_mke (
	experiment_id INT, -- link to experiment
	t TIMESTAMP,
	mesh_builder_cmd TEXT,
	mesh_builder_build TEXT,
	mesh_builder_points INT,
	mesh_builder_triangles INT,
	mesh TEXT,
	tau FLOAT8,
	sigma FLOAT8,
	mu FLOAT8,
	k1 FLOAT8,
	k2 FLOAT8,
	theta FLOAT8,
	rp TEXT,
	coriolis TEXT,
	initial TEXT,
	build TEXT,
	other TEXT,
	cmd TEXT, -- command line
	calc_table VARCHAR(256)
)";
my $table2_index = "CREATE INDEX calc_table_idx ON barvortex_mke(calc_table)";

$dbh->do($id); $dbh->commit();
$dbh->do($table1); $dbh->commit();
$dbh->do($alter); $dbh->commit();
$dbh->do($table2); $dbh->commit();
$dbh->do($table2_index); $dbh->commit();

sub create_uniq_name($) {
	my ($ins) = @_;
	my $hash = md5_hex $ins;
	print STDERR "hash => $hash\n";
	return $hash;
}

sub create_calc_table($) {
	my ($ins) = @_;
	my $uniq_name = create_uniq_name($ins); #md5 of input data
	my $uniq_table_name = "barvortex_mke_$uniq_name";

	my $ans = $dbh->selectall_arrayref("SELECT experiment_id FROM barvortex_mke WHERE calc_table=?",
		undef,$uniq_table_name);
	$dbh->commit();
	my $upd = 0;
	if (scalar @$ans) {
		print STDERR "updating experiment!\n";
		$upd = 1;
	} else {
		print STDERR "create new experiment!\n";
	}

	$dbh->do("DROP TABLE IF EXISTS $uniq_table_name");
	my $table3 = "CREATE TABLE $uniq_table_name (
		step INT UNIQUE,
		t FLOAT8,
		wt FLOAT8,
		nr FLOAT8,
		mn FLOAT8,
		mx FLOAT8,
		v TEXT
	)";
	$dbh->do($table3);
	if (not $upd) {
		# todo > fill experiment
		my $ts = time();
		my $time = POSIX::strftime '%Y-%m-%d %H:%M:%S', localtime($ts);
		$dbh->do($ins,undef,$uniq_table_name,$time);
	}
	$dbh->commit();
	return $uniq_table_name;
}

sub create_insert_string($)
{
	my ($h) = @_;
	my $s = "INSERT INTO barvortex_mke (";
	my $fst = 1;
	while (my ($key, $value) = each(%$h)) {
		$s .= "$key,";
	}
	$s .= "calc_table,t) VALUES(";
	while (my ($key, $value) = each(%$h)) {
		$s .= "'$value',";
	}
	$s .= "?,?)";
	#print STDERR "=> $s\n";
	return $s;
}

sub insert_data
{
	my ($tabname, $step, $t, $wt, $nr, $mn, $mx, $v) = @_;
	my $s = "INSERT INTO $tabname VALUES(?,?,?,?,?,?,?)";
	$dbh->do($s, undef, $step, $t, $wt, $nr, $mn, $mx, $v);
	$dbh->commit();
}

my $read_data = 0;
my $cur = "";
my %fields = ();
$fields{'other'} = `hg id`;

if (not scalar @ARGV) {
	print "ENTER MESH ARGUMENTS!\n";
}

my $mesh_args = $ARGV[0];
my $step_interval = $ARGV[1];

open(PIPE, "./bin/sphere $mesh_args | ");
while(<PIPE>) {
	if (not $read_data) {
		if ($_ =~ m/^#([^:]+):(.*)/) {
			$fields{"mesh_builder_$1"} = $2;
		} else {
			$read_data = 1;
		}
	}

	if ($read_data) {
		$cur .= $_;
	}
};
close(PIPE);
$fields{"mesh"} = $cur;
open(FILE, ">tt.txt");
print FILE $cur;
close(FILE); 

my $cmd = "./bin/test_barvortex -t 1 --task kornev1 -f tt.txt -v 10 --time 100 2>&1 | ";
open(PIPE, $cmd);
print "run $cmd\n";

my $uniq_table_name;
my $t    = 0;
my $nr   = 0;
my $mn   = 0;
my $mx   = 0;
my $step = 1;
my $total = 0;
my $wt   = 0;

$read_data = 0;
$cur = "";

while(<PIPE>) {
	if (not $read_data) {
		if ($_ =~ m/^#.*/) {
			if ($_ =~ m/^#([^:]+):(.*)/) {
				$fields{$1}=$2;
			}
		} else {
			# data begins
			$read_data = 1;
			$uniq_table_name = create_calc_table(create_insert_string(\%fields));
		}
	}

	if ($read_data) {
		if ($_ =~ m/^# end.*\n/) {
			if (not $step_interval or $total % $step_interval == 0) {
				print STDERR "insert $step, $t, $wt, $nr, $mn, $mx \n";
				insert_data($uniq_table_name, $step, $t, $wt, $nr, $mn, $mx, $cur);
				$step += 1;
			}
			$cur = "";
			$total += 1;
		} elsif ($_ =~ m/t=([^;]+); nr=([^;]+); min=([^;]+); max=([^;]+); work=([^;]+);/) {
			$t  = $1;
			$nr = $2;
			$mn = $3;
			$mx = $4;
			$wt = $5;
		} else {
			$cur .= $_;
		}
	}
}

close(PIPE);

