# /bin/bash finalMerge.sh                                                                                        

MCSAMPLE="tH125q_blvu_Yt1_H126toWW WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola TTWJets TTZJets WWWJets_8TeV-madgraph WWZNoGstarJets_8TeV-madgraph WZZNoGstarJets_8TeV-madgraph TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola T_t-channel_TuneZ2star_8TeV-powheg-tauola T_s-channel_TuneZ2star_8TeV-powheg-tauola Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball DYJetsToLL_M-10To50filter_8TeV-madgraph DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"

#MCSAMPLE="TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola"
#MCSAMPLE="tH125q_blvu_Yt1_H126toWW"

echo " ... merging mc files"

VERSION="V05"
FOLDER="results_"${VERSION}

FINALFOLDER=/cmsrm/pc24_2/jorda/data/final${FOLDER}/
PREFOLDER=/cmsrm/pc24_2/jorda/data/${FOLDER}/

mkdir ${FINALFOLDER}
mkdir ${FINALFOLDER}/MC

for SAMPLE in $MCSAMPLE; do

    hadd ${FINALFOLDER}/MC/${SAMPLE}-datasetAll.root ${PREFOLDER}/MC/${SAMPLE}-datasetEEE.root ${PREFOLDER}/MC/${SAMPLE}-datasetEEM.root ${PREFOLDER}/MC/${SAMPLE}-datasetMME.root ${PREFOLDER}/MC/${SAMPLE}-datasetMMM.root ${PREFOLDER}/MC/${SAMPLE}-datasetEE.root ${PREFOLDER}/MC/${SAMPLE}-datasetMM.root

done

#echo " ... merging data files"

#mkdir ${FINALFOLDER}/Data

#hadd ${FINALFOLDER}/Data/DataAll.root ${PREFOLDER}/Data/DoubleElectron-Run2012-datasetEEE.root ${PREFOLDER}/Data/DoubleMu-Run2012-datasetMMM.root ${PREFOLDER}/Data/DoubleElectron-Run2012-datasetEEM.root ${PREFOLDER}/Data/MuEG-Run2012-datasetEEM.root ${PREFOLDER}/Data/DoubleMu-Run2012-datasetMME.root ${PREFOLDER}/Data/MuEG-Run2012-datasetMME.root ${PREFOLDER}/Data/DoubleElectron-Run2012-datasetEE.root ${PREFOLDER}/Data/DoubleMu-Run2012-datasetMM.root