# Crawls GCMT catalog for a event specified with
# its mag and utc event time
class event_finder:
    def __init__(self, gcmtpath):
        """
        At initialization, open the GCMT ascii
        catalogue and convert to a suitable list
        
        """

        gcmt_open=open(gcmtpath)

        # Read everything and split up!
        self.gcmt_read=gcmt_open.readlines()        
        refer=self.gcmt_read[::5]       
        cmts1=self.gcmt_read[1::5]
        cmts2=self.gcmt_read[2::5]
        cmts3=self.gcmt_read[3::5]
        self.event_infos=[]

        nerrgcmt=0
        for i in xrange(len(refer)):
            
            self.evinfo_temp=[]
            # Get reference parameters and cmt sol.
            ev_orig=refer[i];
            ev_cmt1=cmts1[i].split()
            ev_cmt2=cmts2[i].split()
            ev_cmt3=cmts3[i].split()

            # Cut out reference lat/lon/time
            day=ev_orig[5:15];  hrs=ev_orig[15:26];
            dep=ev_orig[42:47]; mag=ev_orig[48:51]

            # Convert origin time to UTC format
            utc_orig=UTCDateTime(day.replace('/','-')+'T'+hrs)

            try:
                utc_gcmt=utc_orig+float(ev_cmt2[1])

                # Relevant for search in catalog
                self.evinfo_temp.append(utc_orig)
                self.evinfo_temp.append(float(mag))
                self.evinfo_temp.append(float(dep))

                # Improved parameters relevant for Axisem
                self.evinfo_temp.append(utc_gcmt)
                self.evinfo_temp.append(ev_cmt1[0]) # event_ID
                self.evinfo_temp.append(ev_cmt2[3]) # latitude                
                self.evinfo_temp.append(ev_cmt2[5]) # longitude
                self.evinfo_temp.append(ev_cmt2[7]) # depth                        

                # CMT solution, relevant for Axisem
                ex=float(10**(int(ev_cmt3[0]))) # Exponent
                self.evinfo_temp.append(float(ev_cmt3[1])*ex)  # Mrr
                self.evinfo_temp.append(float(ev_cmt3[3])*ex)  # Mtt
                self.evinfo_temp.append(float(ev_cmt3[5])*ex)  # Mpp
                self.evinfo_temp.append(float(ev_cmt3[7])*ex)  # Mrt
                self.evinfo_temp.append(float(ev_cmt3[9])*ex)  # Mrp
                self.evinfo_temp.append(float(ev_cmt3[11])*ex) # Mtp

                # Add to huge list
                self.event_infos.append(self.evinfo_temp)
            except:
                nerrgcmt=nerrgcmt+1
                

        print("")
        print("----------------------------------------")
        print("ERRORS IN READING GCMT CATALOG: ",nerrgcmt)
        print("----------------------------------------")
        print("")
        
            
    def search(self, ev_utc, ev_lat, ev_lon, ev_mag):

        """
        This method does the actual searching in the catalog
        spits out a list of fount CMT events
        """        
        self.evref_utc=ev_utc
        self.evref_lat=ev_lat
        self.evref_lon=ev_lon
        self.evref_mag=ev_mag

        self.n=0
        found_events=[]
        for i in xrange(len(self.event_infos)):
            if abs(self.evref_utc-self.event_infos[i][0]) < 60.0:
                found_events.append(self.event_infos[i])               
        found_events = rm_dupl(found_events)
        if len(found_events)>1:
            new_events=[]
            for ev in found_events:
                if (abs(found_events[6]-self.evref_lon)<2) \
                    (abs(found_events[5]-self.evref_lat)<2):                   
                    new_events.append(ev)
            found_events=new_events
        elif len(found_events)==0:
            print("ERROR, GCMT EVENT NOT FOUND!",self.evref_utc)
        if len(found_events)>1:
            print('ERROR, TOO MANY GCMT EVENTS!')
            sys.exit()
                                    
        return found_events  
