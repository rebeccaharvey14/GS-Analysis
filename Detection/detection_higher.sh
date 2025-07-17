#!/bin/bash

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

################### DETECTION ###################
./343.sh $year $dt $time_head $time_tail $namestr
./385.sh $year $dt $time_head $time_tail $namestr
./429.sh $year $dt $time_head $time_tail $namestr
./488.sh $year $dt $time_head $time_tail $namestr
./549.sh $year $dt $time_head $time_tail $namestr
./608.sh $year $dt $time_head $time_tail $namestr
#################################################