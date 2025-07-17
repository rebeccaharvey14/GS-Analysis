# GS-Analysis
Flux rope detection from Grad-Shafranov reconstruction

## contains: ##

**download_preprocess.sh**
*downloads data using pySPEDAS and preprocesses data from FGM, FPI, ESA instrument*
- p2_MMS.py -- fast survey data for plasma & magnetic field data
- data_THM_v2.py -- fast survey data for plasma & magnetic field data
- data_THM_v3.py -- fast survey data for magnetic field data, some plasma data; slow survey data for velocity(+interpolation)
- NOTE: use data_THM_v3.py when velocity data has significant gaps

**detection.sh**
- GS_detection.py -- detection step called by search windows
- Detection/*.sh -- original detection search windows (from Y. Chen)
- Detection/detection_higher*.sh -- higher duration search windows (my addition)

**postprocess.sh**
- 1_combine_raw_result_to_single_file.py -- combines all detection result files
- 2_combineDuplicatedEvent.py -- removes slotting step and uses residuals step as first criteria for event list
- 3_getMoreFluxRopeInfo.py -- gets more detailed info on events
- 4_walen_test.py -- adds more criteria: wait time, shock label, turnTime HCS, days to Exhaust
- 5_WEB_generate_html_form.py -- gets everything in web format; produces csv of event list
- 6_WEB_plotFluxRopeCandidate.py -- plots time series and hodograms of flux rope events
- Output_jz.py -- outputs histogram of axial current density jz

**FluxRopeDetectionPackage**
- fluxrope -- my version of fluxrope.py modified for magnetosheath data ("updated" refers to adaptive datetimeStart(&End), dt inputs)
