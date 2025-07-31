# GS-Analysis
## Flux rope detection from Grad-Shafranov reconstruction ##

### Utilizes adaptive inputs for time range *[start, end]*, time cadence *dt*, and probe name ###

**download_preprocess.sh**
*downloads data using [pySPEDAS](https://github.com/spedas/pyspedas) and preprocesses data into a DataFrame*
- data_MMS.py -- fast survey plasma & magnetic field data
- data_THM.py -- fast survey (when available) plasma & fast survey magnetic field data; slow survey plasma velocity data interpolated + upsampled when fast survey data is unavailable

**detection.sh**
- GS_detection.py			 	 -- detection step called by search windows; utilizes extended-GS equation
- Detection/*.sh 				 -- original detection search windows (from Y. Chen)
- Detection/detection_higher*.sh -- higher duration search windows (from R. Harvey)

**postprocess.sh**
*combines and post-process event candidate lists for one time range*
- 1_combine_raw_results.py 	  -- combines all raw detection files
- 2_combineDuplicatedEvent.py -- eliminates duplicate events, enacts residual criteria, cleans event candidate list
- 3_getMoreFluxRopeInfo.py 	  -- gets more detailed info on identified events
- 4_walen_test.py 			  -- adds more criteria: wait time, shock label, turnTime HCS, days to Exhaust
- 5_generate_csv.py 		  -- produces csv of event list and associated parameters; runs concurrently with 5_generate_html.py
- 5_generate_html.py 		  -- puts everything in web html format; runs concurrently with 5_generate_csv.py
- 6_plotFluxRopeCandidate.py  -- plots time series and hodograms of flux rope events

<!-- - 7_event_page_html.py ---->
<!-- - output_jz.py -- outputs histogram of axial current density *j~z* -->
<!--- - curve_fit_Jz_B_5nT.py ---->

**fluxrope.py**
*R. Harvey's version of fluxrope.py modified for magnetosheath data*
- miscellaneous functions used in detection and postprocessing