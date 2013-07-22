# /bin/bash addOutputfilesJobs.sh                                                                                        

#! Channels
CHANNEL="EEE EEM MME MMM EE MM"
#CHANNEL="MM"

#! Samples
#MCSAMPLE="tH125q_blvu_Yt1_H126toWW WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola TTWJets TTZJets WWWJets_8TeV-madgraph WWZNoGstarJets_8TeV-madgraph WZZNoGstarJets_8TeV-madgraph TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola T_t-channel_TuneZ2star_8TeV-powheg-tauola T_s-channel_TuneZ2star_8TeV-powheg-tauola Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball DYJetsToLL_M-10To50filter_8TeV-madgraph DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"

#MCSAMPLE="DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"
#MCSAMPLE="WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola"
#MCSAMPLE="tH125q_blvu_Yt1_H126toWW"
MCSAMPLE="TTJetsFullyLeptMGDecaysTauola tH125q_blvu_YtMinus1_H126toWW"

echo " ... merging mc files"

VERSION="V05"
FOLDER="results_"${VERSION}

FINALFOLDER=/cmsrm/pc24_2/jorda/data/${FOLDER}/

mkdir ${FINALFOLDER}
mkdir ${FINALFOLDER}/MC

for SAMPLE in $MCSAMPLE; do

    for CH in $CHANNEL; do

	hadd ${FINALFOLDER}/MC/${SAMPLE}-dataset${CH}.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_${VERSION}/${VERSION}/MC/${SAMPLE}/${SAMPLE}_*-dataset${CH}.root

    done

    hadd ${FINALFOLDER}/MC/${SAMPLE}-Counters.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_${VERSION}/${VERSION}/MC/${SAMPLE}/${SAMPLE}_*-Counters.root
    
done

#echo " ... merging data files"

#DATASAMPLE="MuEG-Run2012 DoubleElectron-Run2012 DoubleMu-Run2012"
#DATASAMPLE="MuEG-Run2012"
#DATASAMPLE=""

#mkdir ${FINALFOLDER}/Data

#for SAMPLE in $DATASAMPLE; do

#    for CH in $CHANNEL; do

#	hadd ${FINALFOLDER}/Data/${SAMPLE}-dataset${CH}.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_${VERSION}/${VERSION}/Data/${SAMPLE}AB/${SAMPLE}AB_*-dataset${CH}.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_${VERSION}/${VERSION}/Data/${SAMPLE}C/${SAMPLE}C_*-dataset${CH}.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_${VERSION}/${VERSION}/Data/${SAMPLE}D/${SAMPLE}D_*-dataset${CH}.root

#    done
    
#done

