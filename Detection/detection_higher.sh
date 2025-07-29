#!/bin/bash

probe_str='_MMS1'

########### MMS ###########
time_head='2019-09-27 00:00'
time_tail='2019-09-27 18:00'

time_head='2021-01-27 04:00'
time_tail='2021-01-27 16:00'

time_head='2018-09-02 18:00'
time_tail='2018-09-03 12:00'

time_head='2008-05-26 14:00'
time_tail='2008-05-27 04:00'

time_head='2009-06-08 02:00'
time_tail='2009-06-08 22:00'

time_head='2020-11-28 13:00'
time_tail='2020-11-28 17:00'
namestr='_20201128_20201128'

################### DETECTION ####################
./343.sh $time_head $time_tail $namestr $probe_str
./385.sh $time_head $time_tail $namestr $probe_str
./429.sh $time_head $time_tail $namestr $probe_str
./488.sh $time_head $time_tail $namestr $probe_str
./549.sh $time_head $time_tail $namestr $probe_str
./608.sh $time_head $time_tail $namestr $probe_str
##################################################