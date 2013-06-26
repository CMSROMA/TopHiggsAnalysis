#! /bin/sh


mkdir -p ../results2/MC
lumiEEE=$1
lumiMMM=$2
lumiEEM=$3
lumiMME=$3
echo "Adding weights for eee datasets for " $lumiEEE " pb-1..."
make LatinosAnalyzer
#./LatinosAnalyzer ../results2/MC/tH125q_blvu_Yt1_H126toWW-datasetEEE.root 9.8044e-06*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetEEE.root 3.6945e-06*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-datasetEEE.root 5.34644e-07*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetEEE.root 5.10368e-07*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/TTWJets-datasetEEE.root 1.1834e-06*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/TTZJets-datasetEEE.root 8.27948e-07*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/WWWJets_8TeV-madgraph-datasetEEE.root 3.72708e-07*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/WWZNoGstarJets_8TeV-madgraph-datasetEEE.root 2.84835e-07*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/WZZNoGstarJets_8TeV-madgraph-datasetEEE.root 8.73382e-08*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-datasetEEE.root 4.85102e-05*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root 2.34003e-05*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root 1.60121e-05*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root 3.84675e-06*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root 4.40992e-05*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root 1.72949e-05*$lumiEEE 1 0 
./LatinosAnalyzer ../results2/MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root 7.14418e-06*$lumiEEE 1065353216 0 
./LatinosAnalyzer ../results2/MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-datasetEEE.root 0.000921043*$lumiEEE 0 0 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-datasetEEE.root 0.00215125*$lumiEEE 0 0 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-datasetEEE.root 0*$lumiEEE 81 0 
echo "Adding weights for mmm datasets for " $lumiMMM " pb-1..."
./LatinosAnalyzer ../results2/MC/tH125q_blvu_Yt1_H126toWW-datasetMMM.root 9.8044e-06*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetMMM.root 3.6945e-06*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-datasetMMM.root 5.34644e-07*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetMMM.root 5.10368e-07*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/TTWJets-datasetMMM.root 1.1834e-06*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/TTZJets-datasetMMM.root 8.27948e-07*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/WWWJets_8TeV-madgraph-datasetMMM.root 3.72708e-07*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/WWZNoGstarJets_8TeV-madgraph-datasetMMM.root 2.84835e-07*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/WZZNoGstarJets_8TeV-madgraph-datasetMMM.root 8.73382e-08*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-datasetMMM.root 4.85102e-05*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetMMM.root 2.34003e-05*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetMMM.root 1.60121e-05*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetMMM.root 3.84675e-06*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetMMM.root 4.40992e-05*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetMMM.root 1.72949e-05*$lumiMMM 1 1 
./LatinosAnalyzer ../results2/MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetMMM.root 7.14418e-06*$lumiMMM 1065353216 1 
./LatinosAnalyzer ../results2/MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-datasetMMM.root 0.000921043*$lumiMMM 0 1 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-datasetMMM.root 0.00215125*$lumiMMM 0 1 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-datasetMMM.root 0*$lumiMMM 81 1 
echo "Adding weights for eem datasets for " $lumiEEM " pb-1..."
./LatinosAnalyzer ../results2/MC/tH125q_blvu_Yt1_H126toWW-datasetEEM.root 9.8044e-06*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetEEM.root 3.6945e-06*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-datasetEEM.root 5.34644e-07*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetEEM.root 5.10368e-07*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/TTWJets-datasetEEM.root 1.1834e-06*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/TTZJets-datasetEEM.root 8.27948e-07*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/WWWJets_8TeV-madgraph-datasetEEM.root 3.72708e-07*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/WWZNoGstarJets_8TeV-madgraph-datasetEEM.root 2.84835e-07*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/WZZNoGstarJets_8TeV-madgraph-datasetEEM.root 8.73382e-08*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-datasetEEM.root 4.85102e-05*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetEEM.root 2.34003e-05*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEM.root 1.60121e-05*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEM.root 3.84675e-06*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetEEM.root 4.40992e-05*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEM.root 1.72949e-05*$lumiEEM 1 2 
./LatinosAnalyzer ../results2/MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEM.root 7.14418e-06*$lumiEEM 1065353216 2 
./LatinosAnalyzer ../results2/MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-datasetEEM.root 0.000921043*$lumiEEM 0 2 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-datasetEEM.root 0.00215125*$lumiEEM 0 2 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-datasetEEM.root 0*$lumiEEM 81 2 
echo "Adding weights for mme datasets for " $lumiMME " pb-1..."
./LatinosAnalyzer ../results2/MC/tH125q_blvu_Yt1_H126toWW-datasetMME.root 9.8044e-06*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetMME.root 3.6945e-06*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-datasetMME.root 5.34644e-07*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetMME.root 5.10368e-07*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/TTWJets-datasetMME.root 1.1834e-06*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/TTZJets-datasetMME.root 8.27948e-07*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/WWWJets_8TeV-madgraph-datasetMME.root 3.72708e-07*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/WWZNoGstarJets_8TeV-madgraph-datasetMME.root 2.84835e-07*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/WZZNoGstarJets_8TeV-madgraph-datasetMME.root 8.73382e-08*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-datasetMME.root 4.85102e-05*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetMME.root 2.34003e-05*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetMME.root 1.60121e-05*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetMME.root 3.84675e-06*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetMME.root 4.40992e-05*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetMME.root 1.72949e-05*$lumiMME 1 3 
./LatinosAnalyzer ../results2/MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetMME.root 7.14418e-06*$lumiMME 1065353216 3 
./LatinosAnalyzer ../results2/MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-datasetMME.root 0.000921043*$lumiMME 0 3 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-datasetMME.root 0.00215125*$lumiMME 0 3 
./LatinosAnalyzer ../results2/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-datasetMME.root 0*$lumiMME 81 3 
echo "done weighting."
