#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitAllProcessestHAnalysis.pl -p <prefix> ";}

## example: 
## system("python cmst3_submit_manyfilesperjob.py MC_or_Data MC_Dataset_Name 5 TopHiggsApp 8nh tH isMC?");

print  "submitting mueg ...\n";

system("python cmst3_submit_manyfilesperjob.py Data MuEG-Run2012AB 15 TopHiggsApp 8nh tH 0");
system("python cmst3_submit_manyfilesperjob.py Data MuEG-Run2012C 15 TopHiggsApp 8nh tH 0");
system("python cmst3_submit_manyfilesperjob.py Data MuEG-Run2012D 15 TopHiggsApp 8nh tH 0");
sleep 180;

print  "submitting doubleelectron ...\n";

system("python cmst3_submit_manyfilesperjob.py Data DoubleElectron-Run2012AB 15 TopHiggsApp 8nh tH 0");
system("python cmst3_submit_manyfilesperjob.py Data DoubleElectron-Run2012C 15 TopHiggsApp 8nh tH 0");
system("python cmst3_submit_manyfilesperjob.py Data DoubleElectron-Run2012D 15 TopHiggsApp 8nh tH 0");

sleep 180;

print  "submitting doublemuon ...\n";

system("python cmst3_submit_manyfilesperjob.py Data DoubleMu-Run2012AB 15 TopHiggsApp 8nh tH 0");
system("python cmst3_submit_manyfilesperjob.py Data DoubleMu-Run2012C 15 TopHiggsApp 8nh tH 0");
system("python cmst3_submit_manyfilesperjob.py Data DoubleMu-Run2012D 15 TopHiggsApp 8nh tH 0");

print "\nDONE \n";


