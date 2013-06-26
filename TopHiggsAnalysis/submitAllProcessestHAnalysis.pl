x#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitAllProcessestHAnalysis.pl -p <prefix> ";}

## example: 
## system("python cmst3_submit_manyfilesperjob.py MC_or_Data MC_Dataset_Name 5 TopHiggsApp 8nh tH isMC?");

print  "submitting disobons ...\n";

system("python cmst3_submit_manyfilesperjob.py MC WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola 5 TopHiggsApp 8nh tH 1");

sleep 600;

print  "submitting TT+boson ...\n";

system("python cmst3_submit_manyfilesperjob.py MC TTWJets 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC TTZJets 5 TopHiggsApp 8nh tH 1");

sleep 600;

print  "submitting tribosons ...\n";

system("python cmst3_submit_manyfilesperjob.py MC WWWJets_8TeV-madgraph 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC WWZNoGstarJets_8TeV-madgraph 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC WZZNoGstarJets_8TeV-madgraph 5 TopHiggsApp 8nh tH 1");

sleep 600;

print  "submitting top ...\n";

system("python cmst3_submit_manyfilesperjob.py MC TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC T_t-channel_TuneZ2star_8TeV-powheg-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC T_s-channel_TuneZ2star_8TeV-powheg-tauola 5 TopHiggsApp 8nh tH 1");

sleep 600;

system("python cmst3_submit_manyfilesperjob.py MC Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola 5 TopHiggsApp 8nh tH 1");

print  "submitting w+jets ...\n";
system("python cmst3_submit_manyfilesperjob.py MC WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball 5 TopHiggsApp 8nh tH 1");

print  "submitting dy+jets ...\n";
system("python cmst3_submit_manyfilesperjob.py MC DYJetsToLL_M-10To50filter_8TeV-madgraph 5 TopHiggsApp 8nh tH 1");
system("python cmst3_submit_manyfilesperjob.py MC DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball 10 TopHiggsApp 8nh tH 1");

print "\nDONE \n";


