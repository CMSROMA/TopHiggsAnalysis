# /bin/bash addOutputfilesJobs.sh                                                                                        

#MCSAMPLE="tH125q_blvu_Yt1_H126toWW WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola TTWJets TTZJets WWWJets_8TeV-madgraph WWZNoGstarJets_8TeV-madgraph WZZNoGstarJets_8TeV-madgraph TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola T_t-channel_TuneZ2star_8TeV-powheg-tauola T_s-channel_TuneZ2star_8TeV-powheg-tauola Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball DYJetsToLL_M-10To50filter_8TeV-madgraph DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"

MCSAMPLE="DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball"
#MCSAMPLE=""

echo " ... merging mc files"

mkdir results2/
mkdir results2/MC

for SAMPLE in $MCSAMPLE; do


    hadd results2/MC/${SAMPLE}-datasetEEE.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/MC/${SAMPLE}/${SAMPLE}_*-datasetEEE.root
    hadd results2/MC/${SAMPLE}-datasetEEM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/MC/${SAMPLE}/${SAMPLE}_*-datasetEEM.root
    hadd results2/MC/${SAMPLE}-datasetMME.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/MC/${SAMPLE}/${SAMPLE}_*-datasetMME.root
    hadd results2/MC/${SAMPLE}-datasetMMM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/MC/${SAMPLE}/${SAMPLE}_*-datasetMMM.root

    hadd results2/MC/${SAMPLE}-Counters.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/MC/${SAMPLE}/${SAMPLE}_*-Counters.root

    #hadd results/MC/${SAMPLE}-datasetAll.root results/MC/${SAMPLE}-datasetEEE.root results/MC/${SAMPLE}-datasetEEM.root results/MC/${SAMPLE}-datasetMME.root results/MC/${SAMPLE}-datasetMMM.root   

done

echo " ... merging data files"

#DATASAMPLE="MuEG-Run2012 DoubleElectron-Run2012 DoubleMu-Run2012"
#DATASAMPLE="MuEG-Run2012"
DATASAMPLE=""

mkdir results2/Data

for SAMPLE in $DATASAMPLE; do

    hadd results2/Data/${SAMPLE}-datasetEEE.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}AB/${SAMPLE}AB_*-datasetEEE.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}C/${SAMPLE}C_*-datasetEEE.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}D/${SAMPLE}D_*-datasetEEE.root 

    hadd results2/Data/${SAMPLE}-datasetEEM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}AB/${SAMPLE}AB_*-datasetEEM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}C/${SAMPLE}C_*-datasetEEM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}D/${SAMPLE}D_*-datasetEEM.root 

    hadd results2/Data/${SAMPLE}-datasetMME.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}AB/${SAMPLE}AB_*-datasetMME.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}C/${SAMPLE}C_*-datasetMME.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}D/${SAMPLE}D_*-datasetMME.root 

    hadd results2/Data/${SAMPLE}-datasetMMM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}AB/${SAMPLE}AB_*-datasetMMM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}C/${SAMPLE}C_*-datasetMMM.root /cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_prueba2/tH/Data/${SAMPLE}D/${SAMPLE}D_*-datasetMMM.root 

    #add results/Data/${SAMPLE}-datasetAll.root results/Data/${SAMPLE}-datasetEEE.root results/Data/${SAMPLE}-datasetEEM.root results/Data/${SAMPLE}-datasetMME.root results/Data/${SAMPLE}-datasetMMM.root


done

