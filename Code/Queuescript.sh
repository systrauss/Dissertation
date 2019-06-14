#!/bin/tcsh

#$ -M netID@nd.edu
#$ -m a
#$ -q long
#$ -N Nickname
#$ -o logs/
#$ -t runStart-runEnd
#$ -j y

module load /afs/crc.nd.edu/user/n/nsl/nuclear/startup/nsl
module load root/6.02
root-config --version
cd CodeDirectory
echo "start task ${SGE_TASK_ID}"
date
./main $SGE_TASK_ID FileOut CutFile TimingFile
echo "ended task ${SGE_TASK_ID}"
date

exit $?