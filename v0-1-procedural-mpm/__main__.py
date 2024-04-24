
from equations import globals
from equations import et
from equations import wthrfile
from equations import watbal
from equations import resp
import json
import os



def main():
    print("\n \n********************************************* STARTING SIMULATIONS **********************************************************")

    input_file_name='input.json'
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, input_file_name)
    
    print("****** reading the input files including the parameters and tables  ****")
    globals.Globals.read_file(input_file)
    print("****** starting the simulations ****")

    for i in range(0,globals.PARAM.duration):



        wthrfile.read_wthrfile.get_wth_info(i) # read weather values
    #    print(wthrfile.DAILY_METEO.__dict__)
        print(f"reading weather variables of {globals.Globals.simulation_current_date}")


        output_file = os.path.join(current_dir, f'{globals.PARAM.output_filename}.csv')
        header = "date,developement_stage,LAI,green_leaves_weight,dead_leaves_weight,stem_weight,\
    root_weight,storage_weight,plant_height,top_layer_soil_moisture,rootzone_soil_moisture,\
    rootzone_depth,growth_assimilation,maintenance_respiration,\
    growth_respiration,soil_latent_heat_flux,plant_latent_heat_flux,\
    total_latent_heat_flux,soil_sensible_heat_flux,plant_sensible_heat_flux,\
    total_sensible_heat_flux,soil_PET,plant_PET,total_PET,\
    soil_AET,plant_AET,total_AET\n"
        if i==0:

            with open(output_file, "w") as f:
                    f.write(header)
                    f.close

        outputs_p1=f"{globals.Globals.simulation_current_date},{globals.PARAM.dvs},{globals.PARAM.lai},{globals.WGT.gnleaves},{globals.WGT.ddleaves},\
    {globals.WGT.stem},{globals.WGT.roots},{globals.WGT.storage},{globals.PARAM.hgt},{globals.PARAM.vwc01},\
    {globals.PARAM.vwc02},{globals.PARAM.len02}"


        dET = 0.0
        dETs = 0.0
        dETc = 0.0
        dH = 0.0
        dHs = 0.0
        dHc = 0.0
        print("***calculating heat fluxes***")
        dET, dETs, dETc,dH,dHs,dHc=et.DayHeatFluxes(dET, dETs, dETc,dH,dHs,dHc)

        # daily soil water content (in mm day-1):
        dPotE = dETs * 1000.0 / (2454000.0 * 998.0)
        dPotT = dETc * 1000.0 / (2454000.0 * 998.0)

        print("***calculating soil water balance***")
        watbal.DailyWaterContent(dPotE, dPotT)


        # obtain actual evapotranspiration:
        dEa01 = 0.0
        dEa02 = 0.0
        dEa01, dEa02=watbal.ActualE(dPotE, dEa01, dEa02)
        dTa01 = 0.0
        dTa02 = 0.0
        dTa01, dTa02=watbal.ActualT(0.5, dPotT, dTa01, dTa02)
        
        print("***calculating respiration and assimilations***")
        resp.Grow(dPotT) # plant growth (next stage)

        outputs_p2=f",{globals.ASSIM.gphot},{globals.ASSIM.maint},{globals.ASSIM.growth},{dETs/1000000},{dETc/1000000},\
    {(dETs + dETc)/1000000},{dHs/1000000},{dHc/1000000},{(dHs + dHc)/1000000},{dPotE},{dPotT},{dPotE + dPotT},{dEa01 + dEa02},\
    {dTa01 + dTa02},{dEa01 + dEa02 + dTa01 + dTa02}"


        with open(output_file, "a") as f:
            f.write(outputs_p1+outputs_p2+'\n')

        print(f"***simulations are done for {globals.Globals.simulation_current_date} and the output file is updated*** \n")

        i += 1
        resp.NextDVS()
    



    print("********************************************* END OF SIMULATIONS ********************************************************** \n \n")

if __name__ == "__main__":
    main()     



