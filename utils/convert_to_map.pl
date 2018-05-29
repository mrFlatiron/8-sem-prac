#!/usr/bin/perl

use utf8;
use strict;
use warnings;

opendir my $dir, "." or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

my $convert_script_name = "../../../../utils/addblanklines.awk";

for my $fname (@files)
{
  if (!($fname =~ m/h_.*/)) {next;}

  system ("awk -f \"$convert_script_name\" <$fname > temp");

  system ("mv temp $fname");
}
