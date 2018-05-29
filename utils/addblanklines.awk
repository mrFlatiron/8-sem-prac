/^[[:blank:]]*#/ {next} # ignore comments (lines starting with #)  
NF < 3 {next} # ignore lines which don’t have at least 3 columns  
$2 != prev {printf "\n"; prev=$2} # print blank line  
{print} # print the line 
