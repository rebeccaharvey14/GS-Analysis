#!/bin/bash

probe_str='_MMS1'
time_head='2015-09-10 12:34'
time_tail='2015-09-10 19:15'
namestr='_20150910_20150911'

python3 /home/rharvey/Documents/GS-Analysis/Codes/GS_detection.py $1 $2 $3 $4 11 15 6 1 $5 $6

# $1 = time_head
# $2 = time_tail
# $3 = namestr
# $4 = probe_str