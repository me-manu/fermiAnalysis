/afs/slac.stanford.edu/g/glast/applications/xrootd/PROD/bin/xrdls -l /glast/Data/Flight/Reprocess/P305/ft1/ > tt
/afs/slac.stanford.edu/g/glast/applications/xrootd/PROD/bin/xrdls /glast/Data/Flight/Level1/LPA/prod/5.7/ft1/ > tt2
sed '/gll/!d; s/^.*\(gll.*\).*$/\1/; s#^#/glast/Data/Flight/Reprocess/P305/ft1/#' tt > xrootd-ft1-new.txt
sed '/gll/!d; s/^.*\(gll.*\).*$/\1/; s#^#/glast/Data/Flight/Level1/LPA/prod/5.7/ft1/#' tt2 >> xrootd-ft1-new.txt
rm tt
rm tt2
