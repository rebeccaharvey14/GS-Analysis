#!/bin/bash

########### THM ###########
time_head='2008-05-13 00:00'
time_tail='2008-05-13 12:00'

time_head='2008-05-26 10:00'
time_tail='2008-05-27 17:00'

time_head='2008-05-22 12:00'
time_tail='2008-05-22 22:00'

time_head='2008-05-16 20:00'
time_tail='2008-05-17 10:00'

time_head='2008-05-28 18:00'
time_tail='2008-05-29 06:00'

time_head='2008-06-01 08:00'
time_tail='2008-06-02 04:00'

time_head='2009-06-04 02:00'
time_tail='2009-06-04 23:00'

time_head='2009-06-08 02:00'
time_tail='2009-06-08 22:00'

time_head='2009-06-10 02:00'
time_tail='2009-06-10 16:00'

time_head='2009-06-16 00:00'
time_tail='2009-06-17 00:00'

time_head='2009-06-12 00:00'
time_tail='2009-06-13 00:00'

time_head='2009-06-14 00:00'
time_tail='2009-06-14 18:00'

time_head='2009-06-19 22:00'
time_tail='2009-06-21 00:00'

time_head='2009-05-31 02:00'
time_tail='2009-05-31 21:00'

time_head='2018-09-13 17:00'
time_tail='2018-09-13 22:00'

time_head='2018-09-14 14:00'
time_tail='2018-09-15 04:00'

time_head='2022-01-05 00:00'
time_tail='2022-01-05 06:00'

time_head='2022-01-06 00:00'
time_tail='2022-01-06 06:00'

time_head='2022-01-07 06:00'
time_tail='2022-01-07 12:00'

time_head='2016-12-19 00:00'
time_tail='2016-12-19 17:00'

########### MMS ###########
time_head='2017-09-30 18:00'
time_tail='2017-10-02 08:00'

time_head='2017-10-03 12:00'
time_tail='2017-10-05 06:00'

time_head='2017-10-06 06:00'
time_tail='2017-10-08 00:00'

time_head='2017-10-09 04:00'
time_tail='2017-10-10 22:00'

time_head='2017-11-05 03:30'
time_tail='2017-11-05 12:00'

time_head='2017-11-08 01:30'
time_tail='2017-11-08 06:00'

time_head='2015-09-08 10:43'
time_tail='2015-09-08 18:00'

time_head='2015-09-10 12:34'
time_tail='2015-09-10 19:15'
namestr='_20150910_20150911'

################## DETECTION ###################
./11.sh $year $dt $time_head $time_tail $namestr
./14.sh $year $dt $time_head $time_tail $namestr
./20.sh $year $dt $time_head $time_tail $namestr
./32.sh $year $dt $time_head $time_tail $namestr
./42.sh $year $dt $time_head $time_tail $namestr
./54.sh $year $dt $time_head $time_tail $namestr
./63.sh $year $dt $time_head $time_tail $namestr
./75.sh $year $dt $time_head $time_tail $namestr
./85.sh $year $dt $time_head $time_tail $namestr
./97.sh $year $dt $time_head $time_tail $namestr
./106.sh $year $dt $time_head $time_tail $namestr
./118.sh $year $dt $time_head $time_tail $namestr
./127.sh $year $dt $time_head $time_tail $namestr
./171.sh $year $dt $time_head $time_tail $namestr
./213.sh $year $dt $time_head $time_tail $namestr
./257.sh $year $dt $time_head $time_tail $namestr
./299.sh $year $dt $time_head $time_tail $namestr
#################################################

################### EXAMPLE #####################
# $1 = year
# $2 = dt
# $3 = time_head
# $4 = time_tail

# time_head='2019-09-27 00:00'
# time_tail='2019-09-27 18:00'

# python3 ~/GS-Analysis/Codes/GS_detection.py $1 11 15 6 1 $2 $3 $4
# python3 ~/GS_v3.1.1_MMS.py year($1) 14 23 10 1 dt($2) time_head($3) timetail($4)
