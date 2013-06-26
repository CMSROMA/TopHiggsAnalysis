# /bin/bash addOutputfilesJobs.sh                                                                                        

#MCSAMPLE="tH125q_blvu_Yt1_H126toWW WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola TTWJets TTZJets WWWJets_8TeV-madgraph WWZNoGstarJets_8TeV-madgraph WZZNoGstarJets_8TeV-madgraph TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola T_t-channel_TuneZ2star_8TeV-powheg-tauola T_s-channel_TuneZ2star_8TeV-powheg-tauola Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball DYJetsToLL_M-10To50filter_8TeV-madgraph DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"

#MCSAMPLE="TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola"
MCSAMPLE="DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"

echo " ... merging mc files"

mkdir finalresults2/
mkdir finalresults2/MC

for SAMPLE in $MCSAMPLE; do

    hadd finalresults2/MC/${SAMPLE}-datasetAll.root ../results2/MC/${SAMPLE}-datasetEEE.root ../results2/MC/${SAMPLE}-datasetEEM.root ../results2/MC/${SAMPLE}-datasetMME.root ../results2/MC/${SAMPLE}-datasetMMM.root   

done

echo " ... merging data files"

#DATASAMPLE="MuEG-Run2012 DoubleElectron-Run2012 DoubleMu-Run2012"
#DATASAMPLE="MuEG-Run2012"
DATASAMPLE=""

mkdir finalresults2/Data

for SAMPLE in $DATASAMPLE; do

    hadd finalresults2/Data/${SAMPLE}-datasetAll.root ../results2/Data/${SAMPLE}-datasetEEE.root ../results2/Data/${SAMPLE}-datasetEEM.root ../results2/Data/${SAMPLE}-datasetMME.root ../results2/Data/${SAMPLE}-datasetMMM.root

done

hadd finalresults2/Data/DataAll.root finalresults2/Data/MuEG-Run2012-datasetAll.root finalresults2/Data/DoubleElectron-Run2012-datasetAll.root finalresults2/Data/DoubleMu-Run2012-datasetAll.root