#!/bin/bash
# Download And Preprocess Data 

################################ MMS #######################################
############################################################################
time_head='2017-11-05 03:30'
time_tail='2017-11-05 12:00'
namestr='_20171105_20171106'
probe='1'
year='2015'

python Codes/data_MMS.py $time_head $time_tail $namestr $probe $year


################################ THM #######################################
############################################################################
time_head='2009-06-19 22:00'
time_tail='2009-06-21 00:00'
namestr='_20090619_20090621'
probe='c'
year='2009'

#python Codes/data_THM.py $time_head $time_tail $namestr $probe $year