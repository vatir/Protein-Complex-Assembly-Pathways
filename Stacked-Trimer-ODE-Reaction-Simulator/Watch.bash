#!/bin/bash

while :
do
G=/bin/grep
Q="/usr/local/bin/qstat -u vatir -t" 
RTJ=`$Q|$G vatir|$G -c vatir`
RJ=`$Q|$G vatir|$G -c " R"`
QJ=`$Q|$G vatir|$G -c " Q"`
CJ=`$Q|$G vatir|$G -c " C"`
R1=`$Q|$G "\\[1\\]"|$G -c R`
Q1=`$Q|$G "\\[1\\]"|$G -c Q`
R2=`$Q|$G "\\[2\\]"|$G -c R`
Q2=`$Q|$G "\\[2\\]"|$G -c Q`
R3=`$Q|$G "\\[3\\]"|$G -c R`
Q3=`$Q|$G "\\[3\\]"|$G -c Q`
R4=`$Q|$G "\\[4\\]"|$G -c R`
Q4=`$Q|$G "\\[4\\]"|$G -c Q`
R5=`$Q|$G "\\[5\\]"|$G -c R`
Q5=`$Q|$G "\\[5\\]"|$G -c Q`
R6=`$Q|$G "\\[6\\]"|$G -c R`
Q6=`$Q|$G "\\[6\\]"|$G -c Q`

/usr/bin/clear

echo "Total Running Jobs:  "  $RJ
echo "Total Queued Jobs:   "  $QJ
echo "Total Completed Jobs:"  $CJ
echo "Total Jobs:          "  $RTJ
echo ""
echo "Running [1] Jobs:    "  $R1
echo "Queued  [1] Jobs:    "  $Q1
echo "Running [2] Jobs:    "  $R2
echo "Queued  [2] Jobs:    "  $Q2
echo "Running [3] Jobs:    "  $R3
echo "Queued  [3] Jobs:    "  $Q3
echo "Running [4] Jobs:    "  $R4
echo "Queued  [4] Jobs:    "  $Q4
echo "Running [5] Jobs:    "  $R5
echo "Queued  [5] Jobs:    "  $Q5
echo "Running [6] Jobs:    "  $R6
echo "Queued  [6] Jobs:    "  $Q6

sleep 2
done
