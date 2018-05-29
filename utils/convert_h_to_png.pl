#!/usr/bin/perl

use utf8;
use strict;
use warnings;

opendir my $dir, "." or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

mkdir "images";

for my $fname (@files)
{
  if (!($fname =~ m/h_.*/)) {next;}
  system ("gnuplot -e \"set terminal png; set output \'images/$fname.png\'; set pm3d map; set palette defined (0.1 \'black\', 0.8 \'blue\', 1 \'yellow\'); splot \'$fname\' notitle;\" ");
}
