import os, sys
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import numpy as np
import re 
import subprocess
import geopy.distance
import glob
import pickle

# file = glob.glob('/mnt/users/scratch/leo/scratch/era5/odbs/2/*10393*')
file = glob.glob('/mnt/users/scratch/leo/scratch/era5/odbs/1/era5.conv._*10393*')
print(file)
file = file[0]
# comm = "odb sql 'select distinct statid'  -i FILE --no_alignment ".replace('FILE', file)   
# comm = "odb sql 'select lat,lon'  -i FILE --no_alignment ".replace('FILE', file)   
# comm = "odb sql 'select *'  -i FILE --no_alignment ".replace('FILE', file) 
# comm = "odb sql 'select MIN(date)'  -i FILE --no_alignment ".replace('FILE', file)  # --> '19570901' '19500101'
# comm = "odb sql 'select MAX(date)'  -i FILE --no_alignment ".replace('FILE', file)  # --> '19750630' '19781231'

# comm = "odb sql 'select date, time, varno, obsvalue, vertco_reference_1, vertco_reference_2 WHERE varno=2 and date = 20201210 and vertco_reference_1 = 100000'  -i FILE --no_alignment ".replace('FILE', file)  

comm = "odb sql 'select date, lat, lon, time, varno, obsvalue, vertco_reference_1, vertco_reference_2, report_event1 WHERE (varno=1 and date = 20201210) or (varno=2 and date = 20201210) or (varno=3 and date = 20201210) or (varno=4 and date = 20201210)'  -i FILE --no_alignment ".replace('FILE', file)  # --> '19750630' '19781231'

proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
b = proc.stdout.read().decode('utf-8').split('\n')
df = pd.DataFrame([sub.split("\t") for sub in b])
pickle.dump( b, open( "./lindenberg_hires.p", "wb" ) )

print(b)
# statids = [ eval(c) for c in np.unique( b[1:-1] ) ] 
# print(statids)

'''
Header 8374. Begin offset: 5573955157, end offset: 5574701845, number of rows in block: 10000, byteOrder: same
0. name: type, type: INTEGER, codec: constant, value=263.000000
1. name: expver, type: STRING, codec: constant_string, value='    0001'
2. name: class, type: INTEGER, codec: constant, value=23.000000
3. name: stream, type: INTEGER, codec: constant, value=1025.000000
4. name: andate, type: INTEGER, codec: constant, value=20150527.000000
5. name: antime, type: INTEGER, codec: int32, range=<0.000000,120000.000000>
6. name: reportype, type: INTEGER, codec: int8, range=<16022.000000,16045.000000>
7. name: numtsl@desc, type: INTEGER, codec: constant, value=25.000000
8. name: timeslot@timeslot_index, type: INTEGER, codec: int8, range=<5.000000,17.000000>
9. name: seqno@hdr, type: INTEGER, codec: int32, range=<1675771.000000,7067620.000000>
10. name: source@hdr, type: STRING, codec: constant_string, value='BUFRDATA'
11. name: bufrtype@hdr, type: INTEGER, codec: constant, value=2.000000
12. name: subtype@hdr, type: INTEGER, codec: int8, range=<101.000000,109.000000>
13. name: groupid@hdr, type: INTEGER, codec: constant, value=17.000000
14. name: obstype@hdr, type: INTEGER, codec: constant, value=5.000000
15. name: codetype@hdr, type: INTEGER, codec: int8, range=<35.000000,109.000000>
16. name: sensor@hdr, type: INTEGER, codec: constant, value=0.000000
17. name: statid@hdr, type: STRING, codec: constant_string, value='   10393'
18. name: date@hdr, type: INTEGER, codec: constant, value=20150527.000000
19. name: time@hdr, type: INTEGER, codec: int32, range=<45217.000000,170000.000000>
20. name: report_status@hdr, type: BITFIELD [active:1;passive:1;rejected:1;blacklisted:1;use_emiskf_only:1] , codec: int8, range=<1.000000,4.000000>
21. name: report_event1@hdr, type: BITFIELD [no_data:1;all_rejected:1;bad_practice:1;rdb_rejected:1;redundant:1;stalt_missing:1;qc_failed:1;overcast_ir:1;thinned:1;latlon_corrected:1;stalt_corrected:1] , codec: int8, range=<0.000000,16.000000>
22. name: report_rdbflag@hdr, type: BITFIELD [lat_humon:1;lat_qcsub:1;lat_override:1;lat_flag:2;lat_hqc_flag:1;lon_humon:1;lon_qcsub:1;lon_override:1;lon_flag:2;lon_hqc_flag:1;date_humon:1;date_qcsub:1;date_override:1;date_flag:2;date_hqc_flag:1;time_humon:1;time_qcsub:1;time_override:1;time_flag:2;time_hqc_flag:1;stalt_humon:1;stalt_qcsub:1;stalt_override:1;stalt_flag:2;stalt_hqc_flag:1] , codec: constant, value=0.000000
23. name: lat@hdr, type: REAL, codec: short_real2, range=<52.208000,52.220001>
24. name: lon@hdr, type: REAL, codec: short_real2, range=<14.119800,14.120000>
25. name: lsm@modsurf, type: REAL, codec: short_real2, range=<0.970952,0.971717>
26. name: orography@modsurf, type: REAL, codec: short_real2, range=<56.743793,56.980240>
27. name: windspeed10m@modsurf, type: REAL, codec: short_real2, range=<3.360114,4.602132>
28. name: tsfc@modsurf, type: REAL, codec: short_real2, range=<283.355835,288.638428>
29. name: albedo@modsurf, type: REAL, codec: short_real2, range=<0.152702,0.153323>
30. name: seaice@modsurf, type: REAL, codec: constant, value=0.000000
31. name: snow_depth@modsurf, type: REAL, codec: constant, value=0.000000
32. name: sonde_type@conv, type: INTEGER, codec: constant, value=80.000000
33. name: entryno@body, type: INTEGER, codec: int16, range=<1.000000,21267.000000>
34. name: obsvalue@body, type: REAL, codec: short_real2, range=<-19.779713,329670.156250>, missingValue=-2147483647.000000
35. name: varno@body, type: INTEGER, codec: int8, range=<1.000000,112.000000>
36. name: vertco_type@body, type: INTEGER, codec: constant, value=1.000000
37. name: vertco_reference_1@body, type: REAL, codec: short_real2, range=<720.000000,100610.000000>
38. name: vertco_reference_2@body, type: REAL, codec: real_constant_or_missing, value=NULL, missingValue=-2147483647.000000
39. name: ppcode@conv_body, type: INTEGER, codec: constant, value=0.000000
40. name: datum_anflag@body, type: BITFIELD [final:4;fg:4;depar:4;varqc:4;blacklist:4;ups:1;uvt:1;uhu:1;ut2:1;uh2:1;uv1:1;urr:1;usn:1;usst:1] , codec: int32, range=<0.000000,196608.000000>
41. name: datum_status@body, type: BITFIELD [active:1;passive:1;rejected:1;blacklisted:1;use_emiskf_only:1] , codec: int8, range=<1.000000,12.000000>
42. name: datum_event1@body, type: BITFIELD [vertco_missing:1;obsvalue_missing:1;fg_missing:1;rdb_rejected:1;assim_cld_flag:1;bad_practice:1;vertpos_outrange:1;fg2big:1;depar2big:1;obs_error2big:1;datum_redundant:1;level_redundant:1;not_analysis_varno:1;duplicate:1;levels2many:1;level_selection:1;vertco_consistency:1;vertco_type_changed:1;combined_flagging:1;report_rejected:1;varqc_performed:1;obserror_increased:1;contam_cld_flag:1;contam_rain_flag:1;contam_aerosol_flag:1;bad_emissivity:1;model_cld_flag:1] , codec: int32, range=<0.000000,1048576.000000>
43. name: datum_rdbflag@body, type: BITFIELD [press_humon:1;press_qcsub:1;press_override:1;press_flag:2;press_hqc_flag:1;press_judged_prev_an:2;press_used_prev_an:1;_press_unused_6:6;varno_humon:1;varno_qcsub:1;varno_override:1;varno_flag:2;varno_hqc_flag:1;varno_judged_prev_an:2;varno_used_prev_an:1] , codec: constant, value=0.000000
44. name: biascorr@body, type: REAL, codec: short_real2, range=<-0.030000,0.732950>
45. name: biascorr_fg@body, type: REAL, codec: short_real2, range=<-0.030000,0.732950>
46. name: varbc_ix@body, type: INTEGER, codec: int16, range=<0.000000,5991.000000>
47. name: qc_pge@body, type: REAL, codec: short_real2, range=<0.000000,0.151791>, missingValue=-2147483647.000000
48. name: an_depar@body, type: REAL, codec: short_real2, range=<-140.544449,507.628326>, missingValue=-2147483647.000000
49. name: an_sens_obs@body, type: REAL, codec: constant, value=0.000000
50. name: fg_depar@body, type: REAL, codec: short_real2, range=<-47.051525,537.203491>, missingValue=-2147483647.000000
51. name: an_depar@surfbody_feedback, type: REAL, codec: real_constant_or_missing, value=NULL, missingValue=-2147483647.000000
52. name: fg_depar@surfbody_feedback, type: REAL, codec: short_real2, range=<-3.319958,0.245259>, missingValue=-2147483647.000000
53. name: snow_depth@surfbody_feedback, type: REAL, codec: real_constant_or_missing, value=NULL, missingValue=-2147483647.000000
54. name: snow_density@surfbody_feedback, type: REAL, codec: real_constant_or_missing, value=NULL, missingValue=-2147483647.000000
55. name: datum_status@surfbody_feedback, type: BITFIELD [active:1;passive:1;rejected:1;blacklisted:1;use_emiskf_only:1] , codec: constant, value=4.000000
56. name: datum_sfc_event@surfbody_feedback, type: BITFIELD [statid:1;lsmask:1;stalt_missing:1;obsvalue_missing:1;fg_missing:1;fg2big:1;not_analysis_varno:1;redundant:1;report_rejected:1] , codec: int8, range=<64.000000,128.000000>
57. name: lsm@surfbody_feedback, type: REAL, codec: real_constant_or_missing, value=NULL, missingValue=-2147483647.000000
58. name: stalt@hdr, type: REAL, codec: short_real2, range=<112.000000,115.000000>
59. name: obs_error@errstat, type: REAL, codec: short_real2, range=<0.000016,392.265991>, missingValue=-2147483647.000000
60. name: final_obs_error@errstat, type: REAL, codec: short_real2, range=<0.000016,392.265991>, missingValue=-2147483647.000000
61. name: fg_error@errstat, type: REAL, codec: short_real2, range=<0.000000,77.832458>, missingValue=-2147483647.000000
62. name: eda_spread@errstat, type: REAL, codec: short_real2, range=<0.000000,3.302798>, missingValue=-2147483647.000000
'''