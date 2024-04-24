import math
from equations import globals
from equations import watbal
from equations import wthrfile  
from equations import assim  
import os
# header file: resp.h


# Updates a given value (val) with its rate of change (rate)
def Intgrl(val, rate):
    val = val + rate * globals.PARAM.delta # Euler's integration
    return (val)

def lookup_tables (table_name,dTemp):
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(current_dir, table_name), "r") as f1:
        data = [line for line in f1.readlines() if line.strip()]
    tt=[]
    lk_dict = {}
    for stp,row in enumerate(data[1:]):
        temporary_var=float(row.split(',')[0])
        temporary_rate=float(row.split(',')[1])
        lk_dict[temporary_var]=temporary_rate
        tt.append(temporary_var)   
    if (dTemp>tt[-1])|(dTemp<tt[0]):
        raise ValueError('The lookup value is not in the range provided in the input tables')  
    x1=min(tt, key=lambda x:abs(x-dTemp))
    x1_index=tt.index(x1)
    if dTemp-x1==0:
        dDvr=lk_dict[x1]
    elif (dTemp-x1>0):
        x2=tt[x1_index+1]
        dDvr=((lk_dict[x2]-lk_dict[x1])/(x2-x1))*(dTemp-x1)+lk_dict[x1]

    elif (dTemp-x1<0):
        x2=tt[x1_index-1]
        dDvr=((lk_dict[x2]-lk_dict[x1])/(x2-x1))*(dTemp-x1)+lk_dict[x1]
    return (dDvr)




# Increment to next development growth stage based on current
#   mean daily air temperature
def NextDVS():
    dTemp = (wthrfile.DAILY_METEO.tmin + wthrfile.DAILY_METEO.tmax) / 2.0 # mean
    dDvr = 0.0

    if globals.PARAM.dvs <= 1.0:
        dDvr=lookup_tables (globals.PARAM.tabdvr01,dTemp) # pre-anthesis
    else:
        dDvr=lookup_tables (globals.PARAM.tabdvr02,dTemp) # post-anthesis

#    temp_ref_dvs = RefObject(globals.PARAM.dvs);
    globals.PARAM.dvs=Intgrl(globals.PARAM.dvs, dDvr)

# Leaf area index (LAI) growth (m2 m-2).
def LAIGrowth():
    dSla = 0.0
    dSla=lookup_tables (globals.PARAM.tabsla,globals.PARAM.dvs)

    globals.PARAM.lai = globals.WGT.gnleaves * dSla

# Increment rooting depth. limit is reduction to growth (unitless).
def RootDepthGrowth(limit):
    # increment only if current depth has not reached maximum rooting
    #   depth, soil water content is less than wilting point, and
    #   only during pre-anthesis period
    if (globals.PARAM.len02 < globals.PARAM.maxlen) and (globals.PARAM.vwc02 > globals.PARAM.vwcwp) and (globals.PARAM.dvs < 1.0):
    #    temp_ref_len02 = RefObject(globals.PARAM.len02);
        globals.PARAM.len02=Intgrl(globals.PARAM.len02, limit * globals.PARAM.rootgrate)
    #    globals.PARAM.len02 = temp_ref_len02

# Increment plant height growth.
#   limit is reduction to growth (unitless).
def PlantHeightGrowth(limit):
    # use different base temperature if for post-anthesis
    dBaseTemp = globals.PARAM.basetemp01
    if globals.PARAM.dvs > 1.0:
        dBaseTemp = globals.PARAM.basetemp02

    dTemp = (wthrfile.DAILY_METEO.tmin + wthrfile.DAILY_METEO.tmax) / 2.0 # mean
    dTs = dTemp - dBaseTemp # temperature sum for the day
    if dTs < 0.0:
        dTs = 0.0 # no temp. sum if less than base temp.

#    temp_ref_gTs = RefObject(gTs);
    globals.Globals.gTs=Intgrl(globals.Globals.gTs, dTs) # increment total temperature sum
#    gTs = temp_ref_gTs
    dA = dTs * globals.PARAM.hgtb1 * globals.PARAM.hgtb0 * globals.PARAM.hgtmax * math.exp(-globals.PARAM.hgtb1 * globals.Globals.gTs)
    dB = 1.0 + globals.PARAM.hgtb0 * math.exp(-globals.PARAM.hgtb1 * globals.Globals.gTs)
    dRate = limit * dA / (dB * dB) # daily height growth
#    temp_ref_hgt = RefObject(globals.PARAM.hgt);
    globals.PARAM.hgt=Intgrl(globals.PARAM.hgt, dRate)
#    globals.PARAM.hgt = temp_ref_hgt

# Maintenance respiration: determines the assimilates used for
#   maintenance (g CH20 m-2 day-1).
def RespMaint():
    # respiration coefficients (g CH2O g-1 DM day-1):
    MAINTST = 0.015
    MAINTLV = 0.030
    MAINTRT = 0.015
    MAINTSO = 0.010
    dMaint = MAINTST * globals.WGT.stem + MAINTLV * globals.WGT.gnleaves + MAINTRT * globals.WGT.roots + MAINTSO * globals.WGT.storage
    dTemp = (wthrfile.DAILY_METEO.tmax + wthrfile.DAILY_METEO.tmin) / 2.0 # mean
    dMaint = assim.Q10(dMaint, 2.0, dTemp) # temperature correction
    # need to correct for decrease in metabolic activity with age:
    dTotal = globals.WGT.gnleaves + globals.WGT.ddleaves
    dMaint = dMaint * globals.WGT.gnleaves / dTotal
    if globals.ASSIM.gphot < dMaint:
        dMaint = globals.ASSIM.gphot # not enough assimilates for maintenance

    globals.ASSIM.maint = dMaint # store it

# Determines daily death rate for leaves (day-1)
def LeafDeathRate():
    dDDAge = 0.0 # leaf death due to age
    if globals.PARAM.dvs > 1.0: # death only after anthesis
        dN = 2.0 - globals.PARAM.dvs
        if dN < 0.1:
            dN = 0.1

        dDvr = 0.0
        dTemp = (wthrfile.DAILY_METEO.tmin + wthrfile.DAILY_METEO.tmax) / 2.0
        dDvr=lookup_tables (globals.PARAM.tabdvr02,dTemp)
        dDDAge = dDvr / dN

    dDDShade = 0.0 # leaf death due to self-shading
    LAICR = 4.0 # critical LAI for self-shading death
    if globals.PARAM.lai > LAICR:
        DDLEAF = 0.03 # max. dead leaf coefficient
        dDDShade = DDLEAF * (globals.PARAM.lai - LAICR) / LAICR
        if dDDShade > DDLEAF:
            dDDShade = DDLEAF

    dRdr = dDDAge
    if dDDAge < dDDShade: # select the maximum death rate
        dRdr = dDDShade

    return dRdr

# Growth respiration: determines the assimilates used for
#   growth (g CH20 m-2 day-1).
def RespGrowth():
    globals.ASSIM.growth = globals.ASSIM.gphot - globals.ASSIM.maint # leftover for growth

    # read in tabulated data for DM partitioning:
    dPartst = 0.0
    dPartlv = 0.0
    dPartrt = 0.0
    dPartso = 0.0
    dPartst=lookup_tables (globals.PARAM.tabst,globals.PARAM.dvs)
    dPartlv=lookup_tables (globals.PARAM.tablv,globals.PARAM.dvs)
    dPartrt=lookup_tables (globals.PARAM.tabrt,globals.PARAM.dvs)
    dPartso=lookup_tables (globals.PARAM.tabso,globals.PARAM.dvs)

    # growth respiration coefficients (g CH20 g-1 DM):
    GROST = 1.513
    GROLV = 1.463
    GRORT = 1.444
    GROSO = 1.415
    dGroTot = GROST * dPartst + GROLV * dPartlv + GRORT * dPartrt + GROSO * dPartso
    if dGroTot > 0.0:
        dGroTot = globals.ASSIM.growth / dGroTot # in g DM m-2 day-1

    dGst = dPartst * dGroTot # stem
    dGlv = dPartlv * dGroTot # leaves
    dGrt = dPartrt * dGroTot # roots
    dGso = dPartso * dGroTot # storage organs
    dGdlv = globals.WGT.gnleaves * LeafDeathRate() # leaf death

    # update all DM weights:
#   temp_ref_gnleaves = RefObject(globals.WGT.gnleaves);
    globals.WGT.gnleaves=Intgrl(globals.WGT.gnleaves, dGlv - dGdlv)
#    globals.WGT.gnleaves = temp_ref_gnleaves
#    temp_ref_ddleaves = RefObject(globals.WGT.ddleaves);
    globals.WGT.ddleaves=Intgrl(globals.WGT.ddleaves, dGdlv)
#    globals.WGT.ddleaves = temp_ref_ddleaves
#    temp_ref_stem = RefObject(globals.WGT.stem);
    globals.WGT.stem=Intgrl(globals.WGT.stem, dGst)
#    globals.WGT.stem = temp_ref_stem
#    temp_ref_roots = RefObject(globals.WGT.roots);
    globals.WGT.roots=Intgrl(globals.WGT.roots, dGrt)
#    globals.WGT.roots = temp_ref_roots
#    temp_ref_storage = RefObject(globals.WGT.storage);
    globals.WGT.storage=Intgrl(globals.WGT.storage, dGso)
#    globals.WGT.storage = temp_ref_storage

# Determines the assimilates produced and used for current
#   development stage. Returns the previous parameter (PARAM) values.
#   potT  - potential transpiration (mm day-1).
def Grow(potT):
    dLimit = watbal.GrowthReduction(potT) # water stress

    # assimilates produced:
    globals.ASSIM.gphot = 0.0
    dGphot = assim.DayCanopyAssim()
    dGphot = dGphot * dLimit # reduced growth due to water stress
    if dGphot > 0.0:
        globals.ASSIM.gphot = dGphot * 30.0 / 1000000 # in g CH2O m-2 day-1

    # assimilates used:
    RespMaint()
    RespGrowth()

    LAIGrowth() # update LAI
    PlantHeightGrowth(dLimit) # update plant height
    RootDepthGrowth(dLimit) # update rooting depth

