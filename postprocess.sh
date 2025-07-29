#!/bin/bash

#$1 = time_head
#$2 = time_tail
#$3 = namestr
#$4 = probe_str

time_head='2015-09-10 12:34'
time_tail='2015-09-10 19:15'
namestr='_20150910_20150911'
probe_str='_MMS1'

python3 Codes/1_combine_raw_results.py $time_head $time_tail $namestr $probe_str
python3 Codes/2_combineDuplicatedEvent.py $namestr $probe_str
python3 Codes/3_getMoreFluxRopeInfo.py $namestr $probe_str
python3 Codes/4_walen_test.py $namestr $probe_str
python3 Codes/5_generate_csv.py $namestr $probe_str
# python3 Codes/5_generate_html.py $namestr $probe_str
# python3 Codes/6_plotFluxRopeCandidate.py $time_head
# python3 Codes/7_event_page_html.py
# python3 Codes/output_jz.py $time_head
# python3 Codes/curve_fit_Jz_B_5nT.py $time_head

# Save detection output files
# zip_dir=/home/rharvey/Documents/GS-Analysis/output/${year}/${name_str:1}$probe_str
# zip_file=/home/rharvey/Documents/GS-Analysis/output/${year}/${name_str:1}$probe_str.zip
# echo Saving detection output files to $zip_file
# zip -r $zip_file $zip_dir