# /bin/bash addOutputfilesJobs.sh                                                                                        

MCSAMPLE="tH125q_blvu_Yt1_H126toWW WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola TTWJets TTZJets WWWJets_8TeV-madgraph WWZNoGstarJets_8TeV-madgraph WZZNoGstarJets_8TeV-madgraph TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola T_t-channel_TuneZ2star_8TeV-powheg-tauola T_s-channel_TuneZ2star_8TeV-powheg-tauola Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball DYJetsToLL_M-10To50filter_8TeV-madgraph DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"

#MCSAMPLE="TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola"
#MCSAMPLE=""

echo " ... merging mc files"

mkdir finalresults/
mkdir finalresults/MC

for SAMPLE in $MCSAMPLE; do

    hadd finalresults/MC/${SAMPLE}-datasetAll.root ../results/MC/${SAMPLE}-datasetEEE.root ../results/MC/${SAMPLE}-datasetEEM.root ../results/MC/${SAMPLE}-datasetMME.root ../results/MC/${SAMPLE}-datasetMMM.root   

done

echo " ... merging data files"

#DATASAMPLE="MuEG-Run2012 DoubleElectron-Run2012 DoubleMu-Run2012"
#DATASAMPLE="MuEG-Run2012"
DATASAMPLE=""

mkdir finalresults/Data

for SAMPLE in $DATASAMPLE; do

    hadd finalresults/Data/${SAMPLE}-datasetAll.root finalresults/Data/${SAMPLE}-datasetEEE.root finalresults/Data/${SAMPLE}-datasetEEM.root finalresults/Data/${SAMPLE}-datasetMME.root finalresults/Data/${SAMPLE}-datasetMMM.root


done

