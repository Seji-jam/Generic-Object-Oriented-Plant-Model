# header file: globals.h
#import wthfile.py
import json
import math

# dry matter weight (all in g m-2)
class WGT:

    def __init__(self):
        self.stem = 0
        self.gnleaves = 0
        self.ddleaves = 0
        self.roots = 0
        self.storage = 0


# assimilates (all in g CH2O m-2 day-1)
class ASSIM:

    def __init__(self):
        self.gphot = 0
        self.maint = 0
        self.growth = 0


# the model parameters
class PARAM:

    def __init__(self):
        self.start_date=""
        self.duration=0
        self.delta = 0
#        self.steps = 0
#        self.dvsstop = 0
        self.year = 0
        self.month = 0
        self.day = 0
        self.fname = ""
        self.WeatherFilepath=""
        self.output_filename=""
        self.lat = 0
        self.lon = 0
        self.angb0 = 0
        self.angb1 = 0
        self.refhgt = 0
        self.kwind = 0
        self.keddy = 0
        self.rstA1 = 0
        self.rstA2 = 0
        self.leafwidth = 0
        self.poredist = 0
        self.porosity = 0
        self.vwcwp = 0
        self.vwcsat = 0
        self.vwc01 = 0
        self.vwc02 = 0
        self.len01 = 0
        self.len02 = 0
        self.maxlen = 0
        self.alpha = 0
        self.ksat = 0
        self.dvs = 0
        self.basetemp01 = 0
        self.basetemp02 = 0
        self.lai = 0
        self.hgt = 0
        self.hgtmax = 0
        self.hgtb0 = 0
        self.hgtb1 = 0
        self.rootgrate = 0
    #    self.dm = WGT()
        self.tabdvr01 = ""
        self.tabdvr02 = ""
        self.tabst = ""
        self.tablv = ""
        self.tabrt = ""
        self.tabso = ""
        self.tabsla = ""
        self.tablespath= ""

    # simulation and site properties:

    # flux properties:

    # soil properties:

    # plant properties:

    # tabulated data:

class Globals:
    # global variables:
    gPI = 3.1415926536
    gParam = PARAM() # model parameters
#    gMet = wthfile.DAILY_METEO() # daily weather properties
    gAssim = ASSIM() # daily assimilates produced and used
    gDoy = 0 # current day of year (e.g., Jan 1 = 1)
    gTs = 0.0 # temperature (heat) sum
    simulation_current_date=""

    # Sets the model parameters from data file
    def read_file(input_file):
        #print("****reading the input data****")
        f = open(input_file)
        input_data = json.load(f)

        # simulation and site properties:
        PARAM.start_date=input_data["simulation_start_date"]
        PARAM.duration=input_data["simulation_duration"] 

        PARAM.delta=input_data["Interval"] 
#        PARAM.steps=input_data["TimeStepsStop"] 
#        PARAM.dvsstop=input_data["DevStageStop"] 
#        PARAM.year=input_data["Year"] 
#        PARAM.month=input_data["Month"] 
#        PARAM.day=input_data["Day"]
        PARAM.output_filename=input_data["output_filename"]

        PARAM.fname=input_data["WeatherFilename"] 
        PARAM.lat=input_data["Latitude"] 
        PARAM.lon=input_data["Longitude"] 
        PARAM.angb0=input_data["AngstromIntercept"] 
        PARAM.angb1=input_data["AngstromSlope"] 

        # flux properties:
        PARAM.refhgt=input_data["ReferenceHeight"] 
        PARAM.kwind=input_data["WindAttnCoefficient"] 
        PARAM.keddy=input_data["EddyAttnCoefficient"] 
        PARAM.rstA1=input_data["StomataResA1"] 
        PARAM.rstA2=input_data["StomataResA2"] 
        PARAM.leafwidth=input_data["MeanLeafWidth"] 

        # soil properties:
        PARAM.poredist=input_data["PoreSizeDist"] 
        PARAM.porosity=input_data["Porosity"] 
        PARAM.vwcwp=input_data["WiltingPoint"] 
        PARAM.vwcsat=input_data["SaturationPoint"] 
        PARAM.vwc01=input_data["WaterAmount01"] 
        PARAM.vwc02=input_data["WaterAmount02"] 
        PARAM.len01=input_data["Depth01"] 
        PARAM.len02=input_data["Depth02"] 
        PARAM.maxlen=input_data["MaxRootDepth"] 
        PARAM.alpha=input_data["HydraulicSlope"] 
        PARAM.ksat=input_data["SaturatedHydraulic"] 

        # plant properties:
        PARAM.dvs=input_data["DevStage"] 
        PARAM.basetemp01=input_data["BaseTemp01"] 
        PARAM.basetemp02=input_data["BaseTemp02"] 
        PARAM.lai=input_data["LAI"] 
        PARAM.hgt=input_data["Height"] 
        PARAM.hgtmax=input_data["HeightMax"] 
        PARAM.hgtb0=input_data["HeightIntercept"] 
        PARAM.hgtb1=input_data["HeightSlope"] 
        PARAM.rootgrate=input_data["RootsGrowRate"] 
        WGT.stem=input_data["StemWgt"] 
        WGT.gnleaves=input_data["GreenLeavesWgt"] 
        WGT.ddleaves=input_data["DeadLeavesWgt"] 
        WGT.roots=input_data["RootsWgt"] 
        WGT.storage=input_data["StorageWgt"]

        # file names for tabulated data lookups:
        PARAM.tabdvr01=input_data["tabdvr01"]
        PARAM.tabdvr02=input_data["tabdvr02"]
        PARAM.tabst=input_data["StemTable"]
        PARAM.tablv=input_data["LeavesTable"]
        PARAM.tabrt=input_data["RootsTable"]
        PARAM.tabso=input_data["StorageTable"]
        PARAM.tabsla=input_data["SLATable"]

