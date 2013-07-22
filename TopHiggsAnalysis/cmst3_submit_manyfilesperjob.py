
#! /usr/bin/env python
import os
import sys
import time
import commands
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py process dataset njobs applicationName queue prefix isMC
#######################################
if len(sys.argv) != 8:
    print "usage cmst3_submit_manyfilesperjob.py process datasetname nfileperjob applicationName queue prefix isMC"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
output = dataset
nfileperjob = int(sys.argv[3])
application = sys.argv[4]
queue = sys.argv[5] # choose among cmt3 8nm 1nh 8nh 1nd 1nw 
prefix = sys.argv[6]
isMC = int(sys.argv[7])

if isMC != 0:
    inputlist = "cmst3_53X/"+process+"/"+dataset+".list"
else:
    inputlist = "cmst3_53X/"+process+"/"+dataset+".list"

# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################

diskoutputdir = "/cmsrm/pc24_2/jorda/data/TopHiggs5.3.X_V05/"
diskoutputmain = diskoutputdir+"/"+prefix+"/"+process+"/"+output

# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir -p "+prefix+"/"+process+"/"+output)
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/log/")
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/input/")
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/src/")
#outputroot = outputmain+"/root/"

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm24 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
os.system("cp -r config "+prefix)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+prefix+"/"+process+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(nfileperjob,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)
            

    inputfile.close()

    # prepare the script to run
    outputname = prefix+"/"+process+"/"+output+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('cp -r '+pwd+"/"+prefix+'/config $WORKDIR\n')
    outputfile.write('cp -r '+pwd+"/"+'/data $WORKDIR\n')
    outputfile.write('cp -r '+pwd+"/"+'/elebdtweights $WORKDIR\n')
    outputfile.write('cp -r '+pwd+"/"+'/leptonMVA $WORKDIR\n')

    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc462\n') 
    outputfile.write('cd /afs/cern.ch/work/j/jorda/prueba3/CMSSW_5_3_6/src/TopHiggsAnalysis\n')

    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+" "+str(isMC)+" "+dataset+"\n")

    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm24:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+prefix+"/"+process+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+prefix+"/"+process+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    ijob = ijob+1

    # and now cope with the max number of accesses (3000, keep 2500 for the safe side)
    totfiles = nfileperjob*ijob
    if(totfiles % 1000 == 0):
        time.sleep(500);
        print "sleeping 500 s during a dataset...";
    #time.sleep(1)
    continue
