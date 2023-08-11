# assemblematrix

#### run it

python run_assemblematrix.py

make sure you have all your INPUT right at the beginning of the file **run_assemblematrix.py**

### INPUT (with example input)
* vertices file name with grid nodes x, y, z  
  **vertex_file = 'vertices.h4'**  
* facets file name  
  **facet_file = 'facets.h4'**  
* 1 - invert for dlnVp (y/n), estimated parameter sigma (0.01=1% vp anomalies)  
  **dlnvp = [1, 0.02]**  
* 2 - dlnVs  
  **dlnvs = [0, 0.02]**  
* 3 - dlnQs  
  **dlnQs = [0, 0.3]**  
* 4 - hypoctr corr (in km) first value is on/off, then lon/lat
  **hypocen_corr = [1, 20.0]**  
* 5 - station corr t-P (sec)  
  **sta_corr_tp = [0, 1.0]**  
* 6 - station corr dlnA-P (dimless)  
  **sta_corr_ap = [0, 0.15]**  
* 7 - origin time correction (sec)  
  *1*: one per event 
  *2*: one per cluster  
  **orig_corr = [2, 2.5]**  
* 8 - event correction dlnA-P  
  **ev_corr_ap = [0, 0.5]**  
* 9 - station corr t-S  
  **sta_corr_ts = [0, 1.0]**  
* 10 - station corr dlnA-S  
  **sta_corr_as = [0, 1.0]**  
* 11 - hypoctr corr (in km) first value is on/off, then **depth**  
  **hypocen_corr_depth = [1, 5.0]**

>> ATTENTION: these 11 numbered corrections are stored in **assemblematrix.f** as **jssw** (index) and
**dpar** (correction values). 
Point 11) was added later by MT (26.06.2017) to differentiate between lon/lat correction 
and depth seperately! 
This line: **asparse(kount+3) = -dpar(11) pz/sv ! shift Up** changes the input for depth correction
in assemblematrix.f.<<

* corrections: ellipticity, crust, station elevation, Q-dispersion  
  **corr_switches = [1, 1, 1, 0]**  
* "ident" for output files (mat.ident, out.ident, etc)  
  **out_ident = 'P_T17'**  

* master directory  
  **add_master_mat = '/Users/maria/PhD/Codes/__TEST_DATA/OUTPUT/P_32_85_STFnc_IASP91_2s_03.06.2017_the17/rdrm_cc0.8_crust2_sigma1'**  
  
  **add_list_mat = False**  

* assembled_dir where to copy the files  
  assembled_dir = os.path.join(add_master_mat, 'assemble_dir')

* run assemblematris (what happens if false?)  
  **run_assemblematrix = True**  
* produces 3 figures: nonzero; RHS; dtheor  
  **plot_dtheor = True**  
* plot two maps: stations and number of usage + dT(av)
  **plot_assemble_stations = True**  
* plot columndensity
  **plot_assemble_columndensity = True**  
  
  
## OUTPUT

**mat.< ident >**: assembled matrix file (binary)  
**aux.< ident >**: auxiliary file (ascii)  
  * krow0,kount,rhskrow0,s,ievt,kklust0,stationcode,networkcode,kstat0,iband,krtyp,klust1,igroup,chiw(igroup),dtheor  
  
  * krow0: number of rows
  * kount: nonzero?
  * rhskrow0 (rhs - right hand side): the data are weighted by group weights  
  * s: errors/uncertainty we assigned in pyffproc
  * ievt: event id assigened in add_events.py (raydata_raymaatrix)
  * kklust0: event cluster?
  * stationcode: station name
  * networkcode: netweork
  * kstat0: collecting number of unique stations
  * iband: usually 8 bands
  * krtyp: ray direction
  * klust1: event cluster ?
  * igroup: group
  * chiw (igroup): group weights
  * dtheor: theoretical arrival time
  

>> Note on cluster numbering:  
In the data/matrix files, every event has one or more data clusters; 
such data are labeled with variable **kluster**. Because different
data files may contain the same event, we avoid overlap by redefining
**kluster** from file nfile here as **localkluster**=100*nfile+kluster.
The localklusters are stored in klustlist, which has the same
ordering as the list of event numbers in evlist.
Every source-station path then belongs to an event/cluster combination
which is numbered by **kklust**. The first kklust for each ievt is
also stored, to use when we do not discriminate between clusters of the
same event: **klust1**. << 
  

**ata.< k >.< ident >**: diagonal of  A'A matrix (k=1 for Vp, 2 for Vs, 3 for Q)  

**assemblematrix.out. < ident >**: main diagnostics file, ascii  
**data. < ident >**: ascii, two columns, ndata rows: absolute value (ddif), normalized value (rhs)  
**dum. < ident >**: scratch file for matrix output, written before mat. < ident > (temporary binary output)  
**dumaux. < ident >**: scratch file for auxiliary output, written before aux. < ident > (also temporary ascii output for debugging)  
**inmapc.< k >**: input for k'th kernel density map  
**columndensity.< k >. < ident >**: input file for map.f, scaled ata.1.P (nparms,4) entries: node x,y,z,sum(column(A))/max(A), i.e. relative sensitivity at x,y,z  
**runmapc**: shell script to plot A column sum (spatial kernel density)  
**tex. < ident >**: tex file contains some station/src statistics in tables  


