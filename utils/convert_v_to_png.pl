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
  if (!($fname =~ m/v_.*/))
  {
    next;
  }
  
  system ("gnuplot -e \"set terminal png; set output \'images/$fname.png\'; plot \'$fname\' u 1:2:3:4 with vectors head filled lt 2 notitle;\" ");
}
