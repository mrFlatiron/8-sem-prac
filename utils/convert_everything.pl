#!/usr/bin/perl

use utf8;
use strict;
use warnings;

opendir my $dir, "." or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

chdir "./cdiff-out/gnuplot/h";

system ("../../../../utils/convert_to_map.pl");
system ("../../../../utils/convert_h_to_png.pl");

chdir "../v";

system ("../../../../utils/convert_v_to_png.pl");

chdir "../../..";

chdir "./sokolov-out/gnuplot/h";

system ("../../../../utils/convert_to_map.pl");
system ("../../../../utils/convert_h_to_png.pl");

chdir "../v";

system ("../../../../utils/convert_v_to_png.pl");

chdir "../../..";
