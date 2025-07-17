#!/bin/sh
module load python/2.7.11
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
python3 ~/GS-Analysis/Codes/GS_detection.py $1 106 118 14 1 $2 $3 $4 $5 $6 $7

