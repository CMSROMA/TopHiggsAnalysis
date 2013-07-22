#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitAllSignaltHAnalysis.pl -p <prefix> ";}

## example: 
## system("python cmst3_submit_manyfilesperjob.py MC_or_Data MC_Dataset_Name 5 TopHiggsApp 8nh tH isMC?");

print  "submitting signal ...\n";

system("python cmst3_submit_manyfilesperjob.py MC tH125q_blvu_Yt1_H126toWW 10 TopHiggsApp 8nh V05 1");
system("python cmst3_submit_manyfilesperjob.py MC tH125q_blvu_YtMinus1_H126toWW 10 TopHiggsApp 8nh V05 1");
#system("python cmst3_submit_manyfilesperjob.py MC  5 TopHiggsApp 8nh tH 1");
#system("python cmst3_submit_manyfilesperjob.py MC  5 TopHiggsApp 8nh tH 1");

print "\nDONE \n";


