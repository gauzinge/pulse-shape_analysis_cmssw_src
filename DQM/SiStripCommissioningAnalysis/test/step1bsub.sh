#!/bin/sh
# trap -p "echo got a signal" SIGABRT SIGBUS SIGILL SIGINT SIGKILL SIGQUIT SIGSEGV SIGSTOP 

echo preparing environment
WORKDIR=`pwd`
export CONFDB=cms_trk_r/1A3C5E7G:FIN@CMS_OMDS_TUNNEL
cd $5
eval `scramv1 runtime -sh`
cd $WORKDIR
export SCRATCH=`pwd`
cp $5/OfflineDbClient_new.template OfflineClientScan.py
ex OfflineClientScan.py +":%s/RUNNUMBER/$1" +wq
ex OfflineClientScan.py +":%s/DBPART/$2" +wq
echo checking out data
for i in `ls /afs/cern.ch/work/d/dstrom/CalibScan/${1} | grep ISHA${3}_VFS${4}`; 
do 
    cp /afs/cern.ch/work/d/dstrom/CalibScan/${1}/$i .; 
done
echo running in $SCRATCH
cmsRun OfflineClientScan.py || echo 
echo saving client file as /afs/cern.ch/work/d/dstrom/CalibScan/${1}/SiStripCommissioningClient_00${1}_ISHA${3}_VFS${4}.root
cp SiStripCommissioningClient_00${1}.root /afs/cern.ch/work/d/dstrom/CalibScan/${1}/SiStripCommissioningClient_00${1}_ISHA${3}_VFS${4}.root
echo saving debug.log file as /afs/cern.ch/work/d/dstrom/CalibScan/${1}/00${1}_ISHA${3}_VFS${4}_debug.log
cp *debug.log  /afs/cern.ch/work/d/dstrom/CalibScan/${1}/00${1}_ISHA${3}_VFS${4}_debug.log
