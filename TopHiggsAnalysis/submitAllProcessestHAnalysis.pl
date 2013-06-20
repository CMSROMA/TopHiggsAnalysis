: #!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitAllProcessestHAnalysis.pl -p <prefix> ";}

print  "submitting top...\n";
system("python cmst3_submit_manyfilesperjob.py MC_or_Data_Folder MC_Dataset_Name 5 TopHiggsApp 8nh tH 1");

#system("python cmst3_submit_manyfilesperjob.py Summer12_V14_52X Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola 5 HiggsApp 8nh $prefix 1 2012");

#system("python cmst3_submit_manyfilesperjob.py Summer12_V14_52X T_t-channel_TuneZ2star_8TeV-powheg-tauola 5 HiggsApp 8nh $prefix 1 2012");

#system("python cmst3_submit_manyfilesperjob.py Summer12_V14_52X Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola 5 HiggsApp 8nh $prefix 1 2012");

#system("python cmst3_submit_manyfilesperjob.py Summer12_V14_52X TTJets_TuneZ2star_8TeV-madgraph-tauola 5 HiggsApp 8nh $prefix 1 2012");

#system("python cmst3_submit_manyfilesperjob.py Summer12_V14_52X TTTo2L2Nu2B_8TeV-powheg-pythia6 5 HiggsApp 8nh $prefix 1 2012");

print  "done with top.\n";

sleep 600;

print "\nDONE \n";

