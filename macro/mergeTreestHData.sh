#! /bin/sh 

mkdir -p results_data/merged

echo "Now merging EEE datasets..."
hadd results_data/merged/DoubleElectron-Run2012-datasetEEE.root results_data/testV1/DoubleElectron/*datasetEEE.root

echo "Now merging MMM datasets..."
hadd results_data/merged/DoubleMu-Run2012-datasetMMM.root results_data/testV1/DoubleMu/*datasetMMM.root

echo "Now merging EEM datasets..."
hadd results_data/merged/MuEG-Run2012-datasetEEM.root results_data/testV1/MuEG/*datasetEEM.root
hadd results_data/merged/DoubleElectron-Run2012-datasetEEM.root results_data/testV1/DoubleElectron/*datasetEEM.root

echo "Now merging MME datasets..."
hadd results_data/merged/MuEG-Run2012-datasetMME.root results_data/testV1/MuEG/*datasetMME.root
hadd results_data/merged/DoubleMu-Run2012-datasetMME.root results_data/testV1/DoubleMu/*datasetMME.root

echo "Now merging EE datasets..."
hadd results_data/merged/DoubleElectron-Run2012-datasetEE.root results_data/testV1/DoubleElectron/*datasetEE.root

echo "Now merging MM datasets..."
hadd results_data/merged/DoubleMu-Run2012-datasetMM.root results_data/testV1/DoubleMu/*datasetMM.root
