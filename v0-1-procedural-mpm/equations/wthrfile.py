from enum import Enum
from datetime import datetime
from equations import globals
import os 
# DAILY_METEO
class DAILY_METEO:
    def __init__(self):
        self.avsun = 0
        self.tmax = 0
        self.tmin = 0
        self.rh = 0
        self.wind = 0
        self.rain = 0
        self.starting_wth_step = 0


class read_wthrfile:
    
    def get_wth_info (i):
        current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_wthfile = os.path.join(current_dir, globals.PARAM.fname)
        with open(input_wthfile, "r") as f1:
            data = [line for line in f1.readlines() if line.strip()]
        
        if i==0:
            for stp,row in enumerate(data[1:]):
                wth_step=stp+1 ## +1 added becuase the first line of the weather file is the header
                date_temp=row.split(',')[0]+'-'+row.split(',')[1].zfill(2)+'-'+row.split(',')[2].zfill(2)
                if date_temp==globals.PARAM.start_date:
                    datetime_object = datetime.strptime(globals.PARAM.start_date, '%Y-%m-%d').date()
                    globals.Globals.gDoy=datetime_object.timetuple().tm_yday 
                    DAILY_METEO.starting_wth_step=wth_step
                    globals.Globals.simulation_current_date=date_temp
                    break
                else:
                    print(f"**** skipping simulations for date = {date_temp} ***")   
        else:
            wth_step=DAILY_METEO.starting_wth_step+i
            date_temp=data[wth_step].split(',')[0]+'-'+data[wth_step].split(',')[1].zfill(2)+'-'+data[wth_step].split(',')[2].zfill(2)
            date_temp = datetime.strptime(date_temp, '%Y-%m-%d').date()
            globals.Globals.gDoy=date_temp.timetuple().tm_yday
            globals.Globals.simulation_current_date=date_temp


        DAILY_METEO.avsun=float(data[wth_step].split(',')[3])
        DAILY_METEO.tmax=float(data[wth_step].split(',')[4])
        DAILY_METEO.tmin=float(data[wth_step].split(',')[5])
        DAILY_METEO.rh=float(data[wth_step].split(',')[6])
        DAILY_METEO.wind=float(data[wth_step].split(',')[7])
        DAILY_METEO.rain=float(data[wth_step].split(',')[8])







# wthrfile
class wthrfile:
    def _initialize_instance_fields(self):
        self._daymet_ = DAILY_METEO()
        self._wthrfname_ = ""
    #    self._now_ = clck()

  

    class DateDir(Enum):
        NEXT = 0
        PREV = 1

    def __init__(self):
        self._initialize_instance_fields()

    def DailyMeteo(self):
        return DAILY_METEO(self._daymet_)

 
    def Filename(self):
        return self._wthrfname_
    def Filename(self, fname):
        self._wthrfname_ = fname


