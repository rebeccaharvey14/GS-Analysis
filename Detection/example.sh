#!/bin/sh

time_head='2019-09-27 00:00'
time_tail='2019-09-27 18:00'

python3 ~/GS-Analysis/Codes/GS_detection.py $1 11 15 6 1 $2 $3 $4

#$1 = year
#$2 = dt
#$4 = time_head
#$5 = time_tail