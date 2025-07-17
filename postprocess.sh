#!/bin/bash

#$1 = year
#$2 = dt
#$3 = namestr
#$4 = probe_str

python3 1_combine_raw_result_to_single_file.py $1 $2 $3 $4
python3 2_combineDuplicatedEvent.py $1 $2 $3 $4
python3 3_getMoreFluxRopeInfo.py $1 $2 $3 $4
python3 4_walen_test.py $1
python3 5_generate_csv_html.py $1 $4 $4
# python3 6_plotFluxRopeCandidate.py $1
# python3 7_event_page_html.py
# python3 Output_jz.py $1
# python3 curve_fit_Jz_B_5nT.py $1

# Save detection output files
zip_dir=/home/rharvey/Documents/GS-Analysis/output/${year}/${name_str:1}$probe_str
zip_file=/home/rharvey/Documents/GS-Analysis/output/${year}/${name_str:1}$probe_str.zip
echo Saving detection output files to $zip_file
zip -r $zip_file $zip_dir