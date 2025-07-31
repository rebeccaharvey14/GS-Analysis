#!/bin/bash

########### THM ###########
probe_str='_THMC'

# time_head='2008-05-13 00:00'
# time_tail='2008-05-13 12:00'

# time_head='2008-05-26 10:00'
# time_tail='2008-05-27 17:00'

# time_head='2008-05-22 12:00'
# time_tail='2008-05-22 22:00'

# time_head='2008-05-16 20:00'
# time_tail='2008-05-17 10:00'

# time_head='2008-05-28 18:00'
# time_tail='2008-05-29 06:00'

# time_head='2008-06-01 08:00'
# time_tail='2008-06-02 04:00'

# time_head='2009-06-04 02:00'
# time_tail='2009-06-04 23:00'

# time_head='2009-06-08 02:00'
# time_tail='2009-06-08 22:00'

# time_head='2009-06-10 02:00'
# time_tail='2009-06-10 16:00'

# time_head='2009-06-16 00:00'
# time_tail='2009-06-17 00:00'

# time_head='2009-06-12 00:00'
# time_tail='2009-06-13 00:00'

# time_head='2009-06-14 00:00'
# time_tail='2009-06-14 18:00'

# time_head='2009-06-19 22:00'
# time_tail='2009-06-21 00:00'

# time_head='2009-05-31 02:00'
# time_tail='2009-05-31 21:00'

# time_head='2018-09-13 17:00'
# time_tail='2018-09-13 22:00'

# time_head='2018-09-14 14:00'
# time_tail='2018-09-15 04:00'

# time_head='2022-01-05 00:00'
# time_tail='2022-01-05 06:00'

# time_head='2022-01-06 00:00'
# time_tail='2022-01-06 06:00'

# time_head='2022-01-07 06:00'
# time_tail='2022-01-07 12:00'

# time_head='2016-12-19 00:00'
# time_tail='2016-12-19 17:00'



# ########### MMS ###########
probe_str='_MMS1'

# time_head='2017-09-30 18:00'
# time_tail='2017-10-02 08:00'

# time_head='2017-10-03 12:00'
# time_tail='2017-10-05 06:00'

# time_head='2017-10-06 06:00'
# time_tail='2017-10-08 00:00'

# time_head='2017-10-09 04:00'
# time_tail='2017-10-10 22:00'

time_head='2017-11-05 03:30'
time_tail='2017-11-05 12:00'
namestr='_20171105_20171106'

# time_head='2017-11-08 01:30'
# time_tail='2017-11-08 06:00'

# time_head='2015-09-08 10:43'
# time_tail='2015-09-08 18:00'

# time_head='2015-09-10 12:34'
# time_tail='2015-09-10 19:15'
# namestr='_20150910_20150911'

################### DETECTION ###################
./11.sh $time_head $time_tail $namestr $probe_str
./14.sh $time_head $time_tail $namestr $probe_str
./20.sh $time_head $time_tail $namestr $probe_str
./32.sh $time_head $time_tail $namestr $probe_str
./42.sh $time_head $time_tail $namestr $probe_str
./54.sh $time_head $time_tail $namestr $probe_str
./63.sh $time_head $time_tail $namestr $probe_str
./75.sh $time_head $time_tail $namestr $probe_str
./85.sh $time_head $time_tail $namestr $probe_str
./97.sh $time_head $time_tail $namestr $probe_str
./106.sh $time_head $time_tail $namestr $probe_str
./118.sh $time_head $time_tail $namestr $probe_str
./127.sh $time_head $time_tail $namestr $probe_str
./171.sh $time_head $time_tail $namestr $probe_str
./213.sh $time_head $time_tail $namestr $probe_str
./257.sh $time_head $time_tail $namestr $probe_str
./299.sh $time_head $time_tail $namestr $probe_str
##################################################

################### EXAMPLE ######################
# $1 = time_head
# $2 = time_tail
# $3 = namestr
# $4 = probe_str

# time_head='2019-09-27 00:00'
# time_tail='2019-09-27 18:00'

# python3 /home/rharvey/Documents/GS-Analysis/Codes/GS_detection.py $1 $2 $3 $4 11 15 6 1 $5 $6
# python3 /home/rharvey/Documents/GS-Analysis/Codes/GS_detection.py time_head($1) timetail($2) 14 23 10 1 namestr($3) probe_str($4)
