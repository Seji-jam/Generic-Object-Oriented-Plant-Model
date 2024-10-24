import os
os.chdir(r'D:\Diane_lab\simple_model\generic_crop_model_project\Git_OOD\Generic-Object-Oriented-Plant-Model\v3-0-OOP-operational')

import numpy as np
import pandas as pd
import Root
import Soil
import Leaf
import Canopy
import AtmoWeather
from Input_Setup import Import_Data
from Output_Manager import OutputManager




User_Inputs_file = r'User_Inputs - Cotton.xlsx'
Default_Inputs_file = r'Default_Inputs - Cotton.xlsx'

Weather_File='Maricopa_WW.xlsx'
# soil_layer_inputs='soil_layer_inputs2.xlsx'


# Load crop data from both files
Inputs = Import_Data(User_Inputs_file)
Default_Inputs = Import_Data(Default_Inputs_file)
Inputs.append_data(Default_Inputs)


# soil_data = pd.read_excel(soil_layer_inputs)



# # Extract relevant columns and convert to list of dictionaries
# layers = soil_data[['Layer_Thickness', 'Soil_Moisture', 'Saturated_Hydraulic_Conductivity',
#              'Field_Capacity_', 'Saturated_Soil_Moisture', 'Permanent_Wilting_Point',
#              'Pore_Size_Distribution', 'Air_Entry']].to_dict(orient='records')

# # Create the final dictionary
# Soil_Layer_Property = {'layers': layers}
# planting_depth=Soil_Layer_Property['layers'][0]['Layer_Thickness']
# num_intervals=24
# number_of_layers=len(Soil_Layer_Property['layers'])


psd=0.18
Soil_Layer_Property={'dailydatafile': 'data.txt', 'num_intervals': 10, 'rootdepth': 0.05, 'planting_depth': 0.05, 
 'has_watertable': False, 'number_of_layers': 5, 'Evaporative_Depth':0.25,'alpha':0.18,
 'layers': [{'thick': 0.05,'vwc': 0.3,'ThetaS':0.4,'ThetaFC':0.25,'ThetaPWP':0.10,'KSAT':0.34,'Pore_Size_Distribution':psd, 'PORS':0.43,'texture': {'clay': 12, 'sand': 70, 'om': 0}},
            {'thick': 0.1, 'vwc': 0.3,'ThetaS':0.4,'ThetaFC':0.25,'ThetaPWP':0.10,'KSAT':0.34,'Pore_Size_Distribution':psd, 'PORS':0.43, 'texture': {'clay': 12, 'sand':70, 'om': 0}},
            {'thick': 0.2, 'vwc': 0.3,'ThetaS':0.4,'ThetaFC':0.25,'ThetaPWP':0.10,'KSAT':0.34,'Pore_Size_Distribution':psd, 'PORS':0.43, 'texture': {'clay': 12, 'sand':70, 'om': 0}},
            {'thick': 0.2, 'vwc': 0.3,'ThetaS':0.4,'ThetaFC':0.25,'ThetaPWP':0.10,'KSAT':0.34,'Pore_Size_Distribution':psd, 'PORS':0.43, 'texture': {'clay': 12, 'sand':70, 'om': 0}},
            {'thick': 0.3, 'vwc': 0.3,'ThetaS':0.4,'ThetaFC':0.25,'ThetaPWP':0.10,'KSAT':0.34,'Pore_Size_Distribution':psd, 'PORS':0.43, 'texture': {'clay': 12, 'sand':70, 'om': 0}}]}



# # Initialize the dictionary
# soil_parameters = {}

# # Iterate through the columns and add them to the dictionary
# for column in soil_data.columns:
#     soil_parameters[column] = soil_data[column].tolist()[0]

# if soil_parameters['number_of_soil_layers']* soil_parameters['layer_depth'] < Inputs.max_root_depth:
#     raise RuntimeError("Maximum soil depth cannot be less than maximum root depth")



weather_meteo = AtmoWeather.WeatherMeteo(Weather_File, Inputs.Latitude, Inputs.Simulation_Start_DOY)
weather_data = weather_meteo.process_data()  # This contains weather data for all timesteps
Canopy_object=Canopy.Canopy(Inputs.BTemp_Phen, Inputs.OTemp_Phen, Inputs.CTemp_Phen, Inputs.TempCurve_Res, 
                            Inputs.CropType_Photoperiod, Inputs.StartPhotoperiod_Phase, Inputs.EndPhotoperiod_Phase, 
                            Inputs.Photoperiod_Sensitivity, Inputs.MinThermal_Day_Veg, Inputs.MinThermal_Day_Rep,
             Inputs.CarbonAlloc_Shoot, Inputs.NitrogenAlloc_Shoot, Inputs.Plant_Density, Inputs.Seed_Weight, 
             Inputs.Crop_TypeDet, Inputs.EndSeedNum_DetPeriod, Inputs.Fraction_SeedNitrogen_Remobilizable, Inputs.SLA_Const, Inputs.Min_Specific_Leaf_N, Inputs.MinRootN_Conc,
             Inputs.CarbonCost_NFix, Inputs.MaxN_Uptake, Inputs.MinStemN_Conc,
             Inputs.IniLeafN_Conc, Inputs.MaxPlant_Height,
             Inputs.Legume,
             Inputs.MaxStemGrowth_DS, Inputs.MaxSeedGrowth_DS, Inputs.StemDW_Height, Inputs.Model_TimeStep)

Root_object=Root.Root(Soil_Layer_Property,Inputs.Critical_root_weight_density, Inputs.max_root_depth,  Inputs.Model_TimeStep)
# Soil_object=Soil.Soil( Inputs.Residual_Soil_Moisture, Inputs.Saturated_Soil_Moisture, Inputs.Field_Capacity, 
#                       Inputs.Initial_Soil_Moisture,
#              Inputs.Soil_Depth, Inputs.Top_Layer_Depth, Inputs.Top_Layer_Depth,
#              Inputs.clay_percentage,Inputs.sand_percentage,Inputs.Drainage_Factor,
#               Inputs.initial_soil_temp,  Inputs.soil_resistance_to_evaporation,
#              Inputs.fraction_soil_greater_than_2mm,Inputs.soil_bulk_density,Inputs.organic_N_percentage,
#              Inputs.fraction_N_for_mineralization,
#              Inputs.nitrate_concentration_ppm,Inputs.ammonium_concentration_ppm,
#              Inputs.Fertilizer_applications_count,Inputs.Fertilizer_applications_amount,
#              Inputs.Fertilizer_applications_Days_after_planting,Inputs.Fraction_volatilization,
#              Inputs.Soil_Dynamic_Temperature_Factor)

Soil_object=Soil.Soil(Soil_Layer_Property)


Leaf_object=Leaf.Leaf( Inputs.SLA_Const, Inputs.Min_Specific_Leaf_N, Inputs.Leaf_Blade_Angle, Inputs.Leaf_Width,  Inputs.C3C4_Pathway, Inputs.Ambient_CO2,
                      Inputs.Activation_Energy_JMAX, Inputs.VCMAX_LeafN_Slope, Inputs.JMAX_LeafN_Slope, Inputs.Photosynthetic_Light_Response_Factor )


output_manager = OutputManager()
output_manager.write_header(output_manager.format_header(Soil_object))


variable_check=[]
# Iterate over each timestep's weather data
for day_data in weather_data[:]:
    if Canopy_object.Development_Stage > 2:
        break
    # Extract the necessary weather data for the current timestep
    Year=int(day_data['Year'])
    Solar_Constant=day_data['Solar_Constant']
    Sin_Solar_Declination = day_data['Sin_Solar_Declination']
    Cos_Solar_Declination = day_data['Cos_Solar_Declination']
    Day_Length = day_data['Day_Length']
    Daily_Sin_Beam_Exposure = day_data['Daily_Sin_Beam_Exposure']
    Solar_Radiation = day_data['Solar_Radiation']
    tmax = day_data['Max_Temp']
    tmin = day_data['Min_Temp']
    Vapour_Pressure = day_data['Vapour_Pressure']
    Wind_Speed = day_data['Wind_Speed']
    rain = day_data['Rain']
    doy=int(day_data['Doy'])
    Days_after_planting=doy-Inputs.planting_doy+1
    if doy < Inputs.planting_doy:
        continue
    # if Root_object.Root_Depth/100 >= .75:
    #     print(doy)
    #     break
    print('Development_Stage: ',Canopy_object.Development_Stage, doy)

    # if doy > 165:
    #     break
    # The Leaf Area Update could be part of canopy as well
    Leaf_object.Update_Leaf_Area(Canopy_object.Total_Leaf_Nitrogen,Canopy_object.Nitrogen_Leaf,Canopy_object.CarbonFrac_Veg, Canopy_object.Carbon_determined_LAI)

    Leaf_object.Update_Specific_Leaf_N(Canopy_object.Nitrogen_Leaf)

    # =============================================================================
    # Instantiating Sunlit and shaded leaves for calculating Potential PHOTOSYNTHESIS AND TRANSPIRATION
    # =============================================================================
    
    Leaf_sunlit_object = Leaf.Leaf_sunlit(Leaf_object)
    Leaf_shaded_object = Leaf.Leaf_Shaded(Leaf_object)



    # Call the Leaf method with the current timestep's weather data --> for sunlit parts
    Leaf_sunlit_object.Calculate_Leaf_temp(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)
   
    Leaf_sunlit_object.Calculate_Potential_Photosynthesis(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure, Wind_Speed,Leaf_object.C3C4_Pathway)
    
    Leaf_sunlit_object.Calculate_Potential_Transpiration(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure,Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)
    
    
    # Call the Leaf method with the current timestep's weather data --> for shaded parts
    Leaf_shaded_object.Calculate_Leaf_temp(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)

    Leaf_shaded_object.Calculate_Potential_Photosynthesis(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure, Wind_Speed,Leaf_object.C3C4_Pathway)
 
    Leaf_shaded_object.Calculate_Potential_Transpiration(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure,Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)
    
    
    # =============================================================================
    # Updating Potential Canopy  PHOTOSYNTHESIS AND TRANSPIRATION
    # =============================================================================
    Canopy_object.Update_Canopy_Temp("P",Leaf_sunlit_object.Hourly_Air_Sunlit_Leaf_Temp_diff,Leaf_shaded_object.Hourly_Air_Shaded_Leaf_Temp_diff,
                            Leaf_sunlit_object.Hourly_Sunlit_Leaf_Temp, Leaf_shaded_object.Hourly_Shaded_Leaf_Temp)
    
    Canopy_object.Update_Canopy_Photosyn("P",Leaf_sunlit_object.Hourly_Photosynthesis_Sunlit,Leaf_shaded_object.Hourly_Photosynthesis_Shaded,Day_Length)

    
    Canopy_object.Update_Canopy_Transpiration("P",Leaf_sunlit_object.Hourly_Transpiration_Sunlit,Leaf_shaded_object.Hourly_Transpiration_Shaded,Day_Length)
    # print(Canopy_object.Potential_Canopy_Transpiration,Leaf_object.Hourly_Transpiration_Sunlit,Leaf_object.Hourly_Transpiration_Shaded)
    


    # =============================================================================
    # Calculating Potential Evaporation from soil
    # =============================================================================

    # Soil_object.Calculate_Soil_Water_Content(Root_object.Root_Depth)
    # Soil_object.Calculate_Soil_Water_Content(rain/10,Root_object.Root_Depth)
    # Soil_object.Current_Soil_Water_Status(rain,Root_object.Root_Depth/100)
    
    Soil_object.Calculate_Soil_Evaporation(rain,Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure,Wind_Speed,
                                                     Soil_object.soil_resistance_to_evaporation, Root_object.Root_Depth/100,
                                                     Leaf_object.Leaf_Blade_Angle,.5,Leaf_object.Leaf_area_output['Wind_Ext_Coeff'],
                            Canopy_object.Potential_Canopy_Transpiration)

# Leaf_object.Leaf_area_output['Leaf_Area_Index']

    Soil_object.Calculate_Water_Balance(rain,Canopy_object.Potential_Canopy_Transpiration,Soil_object.potential_evaporation, Root_object.Root_Depth/100)

    # =============================================================================
    # Updating  Photosynthesis, Tranpisraion, and Evaporation if Water Stress ocuurs
    # =============================================================================


    Leaf_sunlit_object.Update_LeafTemp_Photosynthesis_if_WaterStress(Soil_object.Actual_Daily_Transpiration,Canopy_object.Potential_Canopy_Transpiration,Canopy_object.Potential_Canopy_Photosynthesis,Soil_object.Actual_Daily_Evaporation,
                                                                     Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,
                           Soil_object.average_root_zone_water_content,Soil_object.water_supply_for_evaporation, Root_object.Root_Depth/100, 
                           Leaf_sunlit_object.Hourly_Sunlit_Leaf_Temp, Leaf_shaded_object.Hourly_Shaded_Leaf_Temp,
                           Leaf_object.C3C4_Pathway)
   
    
    Leaf_shaded_object.Update_LeafTemp_Photosynthesis_if_WaterStress(Soil_object.Actual_Daily_Transpiration,Canopy_object.Potential_Canopy_Transpiration,Canopy_object.Potential_Canopy_Photosynthesis,Soil_object.Actual_Daily_Evaporation,
                                                                     Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,
                           Soil_object.average_root_zone_water_content,Soil_object.water_supply_for_evaporation, Root_object.Root_Depth/100, 
                           Leaf_sunlit_object.Hourly_Sunlit_Leaf_Temp, Leaf_shaded_object.Hourly_Shaded_Leaf_Temp,
                           Leaf_object.C3C4_Pathway)    

    Canopy_object.Update_Canopy_Transpiration("A",Leaf_sunlit_object.Hourly_Transpiration_Sunlit,Leaf_shaded_object.Hourly_Transpiration_Shaded,Day_Length)
    if abs(Canopy_object.Actual_Canopy_Transpiration - Soil_object.Actual_Daily_Transpiration) > .1 :
        print(Canopy_object.Actual_Canopy_Transpiration ,Soil_object.Actual_Daily_Transpiration)
        raise ValueError()
    # Canopy_object.Water_Stress_Status_Check(Soil_object.water_supply_for_Transpiration,Soil_object.water_supply_for_evaporation,Soil_object.Actual_Daily_Evaporation,
    #                               Soil_object.average_root_zone_water_content,
    #                               Leaf_object.Hourly_Transpiration_Sunlit,Leaf_object.Hourly_Transpiration_Shaded,
    #                               Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, 
    #                               Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, tmax, tmin, Vapour_Pressure,
    #                               Wind_Speed, Canopy_object.Plant_Height,
    #                               Root_object.Root_Depth/100,
    #                               Leaf_object.Hourly_Sunlit_Leaf_Temp, Leaf_object.Hourly_Shaded_Leaf_Temp,
    #                               Soil_object.Hourly_Actual_Soil_Evap,Leaf_object.C3C4_Pathway,
    #                               Leaf_sunlit_object,
    #                               Leaf_shaded_object,
    #                               Leaf_object)
    
    
    # Canopy_object.Update_Data(Soil_object.Actual_Daily_Transpiration,
    #                           Leaf_sunlit_object.Daily_Actual_Photosynthesis_Sunlit,Leaf_shaded_object.Daily_Actual_Photosynthesis_Shaded,
    #                           Leaf_sunlit_object.Daily_Actual_Photosynthesis_Sunlit_DELTA,Leaf_shaded_object.Daily_Actual_Photosynthesis_Shaded_DELTA,
    #                           Leaf_sunlit_object.Daily_Actual_Leaf_Temp_Sunlit,Leaf_shaded_object.Daily_Actual_Leaf_Temp_Shaded)

    Canopy_object.Update_Canopy_Temp("A",Leaf_sunlit_object.Hourly_Actual_Air_Sunlit_Leaf_Temp_Diff,Leaf_shaded_object.Hourly_Actual_Air_Shaded_Leaf_Temp_Diff,
                            Leaf_sunlit_object.Hourly_Actual_Sunlit_Leaf_Temp, Leaf_shaded_object.Hourly_Actual_Shaded_Leaf_Temp)
   
    Canopy_object.Update_Canopy_Photosyn("A",Leaf_sunlit_object.Hourly_Actual_Photosynthesis_Sunlit,Leaf_shaded_object.Hourly_Actual_Photosynthesis_Shaded,Day_Length)
    Canopy_object.Update_Canopy_Photosyn_DELTA(Leaf_sunlit_object.Hourly_Actual_Photosynthesis_Sunlit_DELTA,Leaf_shaded_object.Hourly_Actual_Photosynthesis_Shaded_DELTA,Day_Length)

    


    # Soil_object.Update_Evaporation_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,Soil_object.average_root_zone_moisture,
    #                                               Root_object.Root_Depth,                                                  
    #                                                    Leaf_object.Hourly_Transpiration_Sunlit, Leaf_object.Hourly_Transpiration_Shaded,
    #                                                    Soil_object.Hourly_Actual_Soil_Evap,)
    
    
    
    Soil_object.Calculate_Soil_Temperature(tmin,tmax)


    
    # =============================================================================
    # Thermal units and initializing Canopy Carbon and Nitrogen Accumulation for Biomass formation
    # =============================================================================
    Canopy_object.calculate_thermal_units(Canopy_object.Development_Stage, tmax, tmin, Day_Length, Canopy_object.BTemp_Phen, Canopy_object.OTemp_Phen, Canopy_object.CTemp_Phen, Canopy_object.TempCurve_Res)
    
    
    Canopy_object.calculate_developement_rate(Canopy_object.Development_Stage, Canopy_object.CropType_Photoperiod, Day_Length, Canopy_object.StartPhotoperiod_Phase, Canopy_object.EndPhotoperiod_Phase,
                                              Canopy_object.Photoperiod_Sensitivity, Canopy_object.MinThermal_Day_Veg, Canopy_object.MinThermal_Day_Rep, Canopy_object.Daily_Thermal_Unit)
       
    Canopy_object.Initialize_Biomass_Formation()
    
    Canopy_object.Initialize_Nitrogen_Accumulation()
    
    Canopy_object.Initialize_Carbon_Accumulation()
    
    
    
    # =============================================================================
    # Canopy processes for Respiration, Nitrogen demand and partitioning
    # =============================================================================
    
    Canopy_object.Calculate_Respiration()
    
    Canopy_object.Calculate_Nitrogen_Fixation()
    
    Canopy_object.Calculate_Photo_Assimilates()
    
    Canopy_object.Calculate_Crop_Nitrogen_Demand(Leaf_object.specific_Leaf_n_output['Specific_Leaf_N'])
    
    ###########################################################
    Canopy_object.Calculate_Nitrogen_Partitioning(Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top'])
    
    Canopy_object.Calculate_Seed_Properties(Inputs.Crop_TypeDet, Canopy_object.EndSeedNum_DetPeriod,
                                            Canopy_object.Remobilizable_N_seedgrowth_before_seedfilling,
                                            Canopy_object.Remobilizable_N_seedgrowth_after_seedfilling, Inputs.Fraction_SeedNitrogen_Remobilizable, 
                                            Inputs.Standard_SeedNitrogen_Conc, Canopy_object.Seed_Weight)


    Canopy_object.Calculate_Senescence(Canopy_object.Carbon_determined_LAI, Leaf_object.Leaf_area_output['N_determined_LAI'])
    
    
    # =============================================================================
    # Canopy processes for Carbon flow and partitioning
    # =============================================================================
    
    

    Canopy_object.Calculate_Carbon_Flow()

    Canopy_object.Calculate_Carbon_Flow_Seed_Filling()

    Canopy_object.Calculate_Carbon_Flow_Stem_Growth()
    
    Root_object.Calculate_Root_Senescence(Canopy_object.Root_Carbon,Canopy_object.CarbonFrac_Veg, Canopy_object.MinRootN_Conc, Canopy_object.Root_CarbonReserve, Canopy_object.Nitrogen_Root)

    Canopy_object.Calculate_Carbon_Partitioning(Leaf_object.Leaf_area_output['N_determined_LAI'],Canopy_object.Carbon_determined_LAI, Root_object.nitrogen_determined_Root_Carbon)
    
    Canopy_object.Calculate_Carbon_Production_Rate(Root_object.Root_carbon_loss_rate_senescence)
    
    
    Canopy_object.Update_Carbon_Reserve_Pools()
    
    
    
    Canopy_object.Calculate_Biomass_Change_Rate()
    
    Canopy_object.Calculate_Carbon_Rates()
    
    
    # =============================================================================
    # Soil processes for caculating the state of carbon and Nitrogen for N uptake
    # =============================================================================
    Soil_object.Organic_Carbon_Composition(Root_object.Root_Depth) ## Root_depth must be in cm becuase of an imperical equation for mineralized N
    Soil_object.Organic_Nitrogen_Composition(rain,Root_object.Root_Depth/100) ## Root_depth must be in m becuase it is used with total soil profile depth which is in m

    Soil_object.Calculate_Nitrogen_Uptake(Canopy_object.Nitrogen_Demand,Canopy_object.NitrogenFixation_Reserve_Pool_ChangeRate,
                                          Root_object.Root_Depth/100) ## Root_depth must be in m becuase it is used with total soil profile depth which is in m


    # =============================================================================
    # Canopy processes for updating Nitrogen dynamic based on the nitrogen absorption from soil
    # =============================================================================
 
    Canopy_object.Calculate_Nitrogen_Dynamics(Soil_object.nitrogen_uptake)
    
    Canopy_object.Calculate_Seed_Number_Rate(Soil_object.nitrogen_uptake,Root_object.Root_carbon_loss_rate_senescence)
    
    Canopy_object.Calculate_Nitrogen_Accumulation_Rate(Soil_object.nitrogen_uptake,Root_object.Root_nitrogen_loss_rate_senescence,Inputs.Standard_SeedNitrogen_Conc)
    
    Canopy_object.Calculate_Leaf_Area_ChangeRate(Canopy_object.Carbon_determined_LAI,Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Bottom'],Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff'])

    Root_object.Calculate_Rooting_Depth(Canopy_object.RootWeight_Rate, Canopy_object.LiveRoot_Dry_Weight, Canopy_object.DeadRoot_Dry_Weight)

    Canopy_object.Calculate_Respiration_Rates(Soil_object.nitrogen_uptake, Root_object.Root_carbon_loss_rate_senescence)
    
    Canopy_object.Calculate_Carbon_Nitrogen_Returns(Root_object.Root_carbon_loss_rate_senescence, Root_object.Root_nitrogen_loss_rate_senescence,Soil_object.average_soil_temperature)
        

    # =============================================================================
    # updating Soil water balance, carbon, and Nitrogen
    # =============================================================================
    # Canopy_object.Check_Carbon_Nitrogen_Balance(Soil_object.tnupt)

    # Soil_object.Calculate_Water_Balance(Days_after_planting,rain,Canopy_object.Actual_Canopy_Transpiration,Root_object.Root_Depth,Root_object.root_depth_growth_rate)
    # Soil_object.Calculate_Water_Balance(rain,Canopy_object.Actual_Canopy_Transpiration, Root_object.Root_Depth/100)
    # print(Canopy_object.Actual_Canopy_Transpiration,'TR')
    Soil_object.Calculate_Soil_Organic_Nitrogen_ChangeRate(Days_after_planting,Root_object.root_depth_growth_rate/100,Root_object.Root_Depth/100) #Root_depth must be in m because it is used with total depth which is in m...growth rate should be compatible and that's why it should be in m/day
    
    
    Soil_object.Calculate_Soil_Organic_Carbon_ChangeRate(Canopy_object.LitterCarbon_Total,Canopy_object.LitterNitrogen_Total)
    
    

    # =============================================================================
    # updating state variables 
    # =============================================================================
    Canopy_object.Update_State_Variables(Root_object.Root_nitrogen_loss_rate_senescence)
    
    Root_object.Update_State_Variables()
    Soil_object.Update_State_Variables()

    
    # =============================================================================
    # writing outputs
    # =============================================================================
    formatted_data = output_manager.format_data(day_data, Canopy_object, Leaf_object, Root_object, Soil_object)
    output_manager.append_data(formatted_data)
    
    df_check=pd.DataFrame({'Development_Stage':[Canopy_object.Development_Stage],'Fraction_Nitrogen_to_Shoot':Canopy_object.Fraction_Nitrogen_to_Shoot})
    variable_check.append(df_check)

# =============================================================================
# 
# =============================================================================
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
results=pd.read_csv('MODEL_OUTPUTS.csv')
results['Date'] = pd.to_datetime(results['Year'].astype(int).astype(str) + results['doy'].astype(int).astype(str), format='%Y%j')


measured_height_approximations = pd.DataFrame({'Date': ['2023-08-03','2023-08-31'], 'Plant_Height': [.6,.7]})
measured_height_approximations['Date'] = pd.to_datetime(measured_height_approximations['Date'])

results.columns


variable_check=pd.concat(variable_check,axis=0,ignore_index=True)


fig, ax = plt.subplots(figsize=(12, 8))
sns.lineplot(x=variable_check['Development_Stage'],y=variable_check['Fraction_Nitrogen_to_Shoot'],color='brown',ax=ax)
fig, ax = plt.subplots(figsize=(12, 8))
sns.lineplot(x=results['Development_Stage'],y=results['LAI'],color='brown',ax=ax)


# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Fraction_Stem_Carbon'],color='brown',ax=ax)
# sns.lineplot(x=results['Development_Stage'],y=results['Fraction_Seed_Carbon'],color='k',ax=ax)
# sns.lineplot(x=results['Development_Stage'],y=results['Fraction_Leaf_Carbon'],color='green',ax=ax)
#         # self.Leaf_Carbon_ChangeRate = 12.0 / 44.0 * self.Photo_Assimilate *
#         # self.Fraction_Carbon_to_Shoot * self.Fraction_Leaf_Carbon *
#         # self.Growth_Efficiency_Veg - self.Leaf_Carbon_Loss_ChangeRate

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Carbon_Flow_to_Stem'],color='brown',ax=ax)

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Relative_Shoot_Activity'],color='brown',ax=ax)


# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Carbon_Shoot'],color='brown',ax=ax)

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Root_Carbon'],color='brown',ax=ax)

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['nitrogen_determined_Root_Carbon'],color='brown',ax=ax)

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Nitrogen_Root'],color='brown',ax=ax)

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Root_CarbonReserve'],color='brown',ax=ax)

# # Specific_Leaf_N_Top,\
# # Nitrogen_Demand_Activity_Driven,\

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Nitrogen_Demand_Activity_Driven'],color='brown',ax=ax)

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Nitrogen_uptake'],color='brown',ax=ax)


# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Development_Stage'],y=results['Total_Nitrogen_uptake'],color='brown',ax=ax)



# sns.lineplot(x=results['Date'],y=results['Plant_Height'])
# sns.scatterplot(x=measured_height_approximations['Date'],y=measured_height_approximations['Plant_Height'])

# fig, ax = plt.subplots(figsize=(12, 8))
# ax2 = ax.twinx()
# sns.scatterplot(x=results['Date'],y=results['Rain'],ax=ax2)
# sns.lineplot(x=results['Date'],y=results['soil_moisture_1'],ax=ax)


# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['Root_Depth'])

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['soil_moisture_2'])
# sns.lineplot(x=results['Date'],y=results['soil_moisture_3'],color='red')
# ax2 = ax.twinx()
# sns.scatterplot(x=results['Date'],y=results['Rain'],ax=ax2)


# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['soil_moisture_5'])

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['Nitrogen_Leaf'])

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['Potential_Canopy_Transpiration'],color='red')
# sns.lineplot(x=results['Date'],y=results['Actual_Canopy_Transpiration'],color='blue')

# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['Potential_Canopy_Photosynthesis'],color='red')
# sns.lineplot(x=results['Date'],y=results['Actual_Canopy_Photosynthesis'],color='blue')


# fig, ax = plt.subplots(figsize=(12, 8))
# sns.lineplot(x=results['Date'],y=results['Actual_Daily_Evaporation'],color='blue')

# results.columns


# import os
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import re


file_path_csv = r'D:\Trees\trees272_corrected\LAI_measured_data.csv'
measured_data = pd.read_csv(file_path_csv)
replacements = {
    'UGA230': 'UGA230',
    'Pronto': 'Pronto',
    'Tipo Chaco': 'TipoChaco',
    'Virescent nankeen': 'Virescentnankeen',
    'Coker 310': 'Coker310',
    'DeltaPine 16': 'DeltaPine16'
}
measured_data['Entry'] = measured_data['Entry'].replace(replacements)
measured_data['Treatment'] = measured_data['Treatment'].str.lower()

measured_data['Date and Time'] = pd.to_datetime(measured_data['Date and Time'])
measured_data['Date'] = measured_data['Date and Time'].dt.date
measured_data = measured_data.loc[measured_data['Position'] != 'Top']

measured_data_avg_lai = measured_data.groupby([ 'Treatment', 'Date']).agg(
    LAI_mean=('Leaf Area Index [LAI]', 'mean'),
    LAI_std=('Leaf Area Index [LAI]', 'std')
).reset_index()

measured_data_avg_lai=measured_data_avg_lai.loc[measured_data_avg_lai['Treatment']=='ww']

file_path_csv = r'D:\Diane_lab\simple_model\generic_crop_model_project\Git_OOD\Generic-Object-Oriented-Plant-Model\v3-0-OOP-operational\MODEL_OUTPUTS.csv'
simulated_data = pd.read_csv(file_path_csv, index_col=False)

simulated_data['Date'] = pd.to_datetime(simulated_data['Year'].astype(int).astype(str) + simulated_data['doy'].astype(int).astype(str), format='%Y%j')
simulated_data.columns

fig, ax = plt.subplots(figsize=(12, 8))
sns.lineplot(x='Date', y='LAI', data=simulated_data, linewidth=2,
              linestyle='-', color='k', label='Simulated LAI', ax=ax)

eb1 = ax.errorbar(
    x=measured_data_avg_lai['Date'], 
    y=measured_data_avg_lai['LAI_mean'], 
    yerr=measured_data_avg_lai['LAI_std'], 
    fmt='D', 
    markersize=10,
    markerfacecolor='white', 
    markeredgecolor='red', 
    ecolor='red',
    label='Measured LAI', 
    capsize=5,
    capthick=2,
    elinewidth=2,
)
eb1[-1][0].set_linestyle('--')

ax.set_ylim(0.0, 6.8)
ax.set_xlabel('')
ax.set_ylabel('LAI (Leaf Area Index)', fontsize=20)
ax.tick_params(axis='x', rotation=45, labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.grid(True, linestyle='--', linewidth=0.7)


ax.legend(fontsize=20, loc='upper left')













