
import numpy as np
import pandas as pd
import Root
import Soil
import Leaf
import Canopy
import AtmoWeather
from Input_Setup import Import_Data


User_Inputs_file = r'User_Inputs.xlsx'
Default_Inputs_file = r'Default_Inputs.xlsx'
Weather_File='weatherfile.xlsx'

# Load crop data from both files
Inputs = Import_Data(User_Inputs_file)
Default_Inputs = Import_Data(Default_Inputs_file)
Inputs.append_data(Default_Inputs)

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

Root_object=Root.Root(Inputs.Critical_root_weight_density, Inputs.max_root_depth,  Inputs.Model_TimeStep, Inputs.Soil_Evaporative_Depth)
Soil_object=Soil.Soil( Inputs.Residual_Soil_Moisture, Inputs.Saturated_Soil_Moisture, Inputs.Field_Capacity, 
                      Inputs.Initial_Soil_Moisture,
             Inputs.Soil_Depth, Inputs.Top_Layer_Depth, Inputs.Soil_Evaporative_Depth,
             Inputs.clay_percentage,Inputs.sand_percentage,Inputs.Drainage_Factor,
              Inputs.initial_soil_temp,  Inputs.soil_resistance_to_evaporation,
             Inputs.fraction_soil_greater_than_2mm,Inputs.soil_bulk_density,Inputs.organic_N_percentage,
             Inputs.fraction_N_for_mineralization,
             Inputs.nitrate_concentration_ppm,Inputs.ammonium_concentration_ppm,
             Inputs.Fertilizer_applications_count,Inputs.Fertilizer_applications_amount,
             Inputs.Fertilizer_applications_DAP,Inputs.Fraction_volatilization,
           )
Leaf_object=Leaf.Leaf( Inputs.SLA_Const, Inputs.Min_Specific_Leaf_N, Inputs.Leaf_Blade_Angle, Inputs.Leaf_Width,  Inputs.C3C4_Pathway, Inputs.Ambient_CO2,
                      Inputs.Activation_Energy_JMAX, Inputs.VCMAX_LeafN_Slope, Inputs.JMAX_LeafN_Slope, Inputs.Photosynthetic_Light_Response_Factor )


# Iterate over each timestep's weather data
for day_data in weather_data:
    # Extract the necessary weather data for the current timestep
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
    doy=day_data['Doy']
    Days_after_planting=doy-Inputs.planting_doy+1
    


    # The Leaf Area Update could be part of canopy as well
    Leaf_object.Update_Leaf_Area(Canopy_object.Total_Leaf_Nitrogen,Canopy_object.CarbonFrac_Veg, Canopy_object.Carbon_determined_LAI)
    Leaf_object.Update_Specific_Leaf_N(Canopy_object.Total_Leaf_Nitrogen)

    # =============================================================================
    # Instantiating Sunlit and shaded leaves for calculating Potential PHOTOSYNTHESIS AND TRANSPIRATION
    # =============================================================================
    
    Leaf_sunlit_object = Leaf.Leaf_sunlit(Leaf_object)
    Leaf_shaded_object = Leaf.Leaf_Shaded(Leaf_object)



    # Call the Leaf method with the current timestep's weather data --> for sunlit parts
    Leaf_sunlit_object.Calculate_Leaf_temp(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)
   
    Leaf_sunlit_object.Calculate_Potential_Photosynthesis(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Leaf_object.C3C4_Pathway)
    Leaf_sunlit_object.Calculate_Potential_Transpiration(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)
    
    
    # Call the Leaf method with the current timestep's weather data --> for shaded parts
    Leaf_shaded_object.Calculate_Leaf_temp(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)

    Leaf_shaded_object.Calculate_Potential_Photosynthesis(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Leaf_object.C3C4_Pathway)
 
    Leaf_shaded_object.Calculate_Potential_Transpiration(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,Canopy_object.Plant_Height,Leaf_object.C3C4_Pathway)
    
    
    # =============================================================================
    # Updating Potential Canopy  PHOTOSYNTHESIS AND TRANSPIRATION
    # =============================================================================
    Canopy_object.Update_Canopy_Temp("P",Leaf_object.Hourly_Air_Sunlit_Leaf_Temp_diff,Leaf_object.Hourly_Air_Shaded_Leaf_Temp_diff,
                            Leaf_object.Hourly_Sunlit_Leaf_Temp, Leaf_object.Hourly_Shaded_Leaf_Temp)
    
    Canopy_object.Update_Canopy_Photosyn("P",Leaf_object.Hourly_Photosynthesis_Sunlit,Leaf_object.Hourly_Photosynthesis_Shaded,Day_Length)
   
    
    Canopy_object.Update_Canopy_Transpiration("P",Leaf_object.Hourly_Transpiration_Sunlit,Leaf_object.Hourly_Transpiration_Shaded,Day_Length)
    
    


    # =============================================================================
    # Calculating Potential Evaporation from soil
    # =============================================================================

    Soil_object.Calculate_Soil_Water_Content()
    
    Soil_object.Calculate_Soil_Potential_Evaporation(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,
                                                     Soil_object.soil_resistance_to_evaporation,Leaf_object.Leaf_Blade_Angle,
                                                     Leaf_object.Leaf_area_output['Total_LAI'],Leaf_object.Leaf_area_output['Wind_Ext_Coeff'],
                            Leaf_object.Hourly_Transpiration_Shaded,Leaf_object.Hourly_Transpiration_Sunlit)




    # =============================================================================
    # Updating  Photosynthesis, Tranpisraion, and Evaporation if Water Stress ocuurs
    # =============================================================================


    Leaf_sunlit_object.Update_LeafTemp_Photosynthesis_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,
                           Soil_object.Current_Soil_Water_Content_Top_Layer,Soil_object.Soil_Evaporative_Depth,Soil_object.Root_Depth, 
                           Leaf_object.Hourly_Sunlit_Leaf_Temp, Leaf_object.Hourly_Shaded_Leaf_Temp,
                           Soil_object.Hourly_Soil_Evap,Leaf_object.C3C4_Pathway)
    Leaf_shaded_object.Update_LeafTemp_Photosynthesis_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,
                           Soil_object.Current_Soil_Water_Content_Top_Layer,Soil_object.Soil_Evaporative_Depth,Soil_object.Root_Depth, 
                           Leaf_object.Hourly_Sunlit_Leaf_Temp, Leaf_object.Hourly_Shaded_Leaf_Temp,
                           Soil_object.Hourly_Soil_Evap,Leaf_object.C3C4_Pathway)
    

    Canopy_object.Update_Canopy_Temp("A",Leaf_object.Hourly_Actual_Air_Sunlit_Leaf_Temp_Diff,Leaf_object.Hourly_Actual_Air_Shaded_Leaf_Temp_Diff,
                            Leaf_object.Hourly_Actual_Sunlit_Leaf_Temp, Leaf_object.Hourly_Actual_Shaded_Leaf_Temp)
    Canopy_object.Update_Canopy_Photosyn("A",Leaf_object.Hourly_Actual_Photosynthesis_Sunlit,Leaf_object.Hourly_Actual_Photosynthesis_Shaded,Day_Length)
    Canopy_object.Update_Canopy_Transpiration("A",Leaf_object.Hourly_Actual_Transpiration_Sunlit,Leaf_object.Hourly_Actual_Transpiration_Shaded,Day_Length)
    
    


    Soil_object.Update_Evaporation_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,Soil_object.Current_Soil_Water_Content_Top_Layer,
                                                  Soil_object.Soil_Evaporative_Depth,Soil_object.Root_Depth,                                                  
                                                       Leaf_object.Hourly_Transpiration_Sunlit, Leaf_object.Hourly_Transpiration_Shaded,
                                                       Soil_object.Hourly_Soil_Evap,)
    
    
    
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
    
    Canopy_object.Calculate_Nitrogen_Partitioning(Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top_Increment'])
    
    Canopy_object.Calculate_Seed_Properties(Inputs.Crop_TypeDet, Canopy_object.EndSeedNum_DetPeriod,
                                            Canopy_object.Remobilizable_N_seedgrowth_before_seedfilling,
                                            Canopy_object.Remobilizable_N_seedgrowth_after_seedfilling, Inputs.Fraction_SeedNitrogen_Remobilizable, 
                                            Inputs.Standard_SeedNitrogen_Conc, Canopy_object.Seed_Weight)


    Canopy_object.Calculate_Senescence(Canopy_object.Carbon_determined_LAI, Leaf_object.Leaf_area_output['N_determined_LAI'], Canopy_object.Model_TimeStep)
    
    
    # =============================================================================
    # Canopy processes for Carbon flow and partitioning
    # =============================================================================
    
    Canopy_object.Calculate_Carbon_Flow()

    Canopy_object.Calculate_Carbon_Flow_Seed_Filling()

    Canopy_object.Calculate_Carbon_Flow_Stem_Growth()

    Root_object.Calculate_Root_Senescence(Canopy_object.Root_Carbon,Canopy_object.CarbonFrac_Veg, Canopy_object.MinRootN_Conc, Canopy_object.Root_CarbonReserve, Canopy_object.Nitrogen_Root)


    Canopy_object.Calculate_Carbon_Partitioning(Leaf_object.Leaf_area_output['N_determined_LAI'],Canopy_object.Carbon_determined_LAI, Root_object.nitrogen_determined_Root_Carbon)
    
    
    Canopy_object.Update_Carbon_Reserve_Pools()
    
    
    Canopy_object.Calculate_Carbon_Production_Rate(Root_object.Root_carbon_loss_rate_senescence)
    
    Canopy_object.Calculate_Biomass_Change_Rate()
    
    Canopy_object.Calculate_Carbon_Rates()
    
    
    # =============================================================================
    # Soil processes for caculating the state of carbon and Nitrogen for N uptake
    # =============================================================================
    Soil_object.Soil_Water_Components(rain, Root_object.root_depth_current, Canopy_object.Actual_Canopy_Transpiration)

    Soil_object.Calculate_Soil_N_Dynamics(Days_after_planting)
    Soil_object.Calculate_Nitrogen_Uptake(Canopy_object.Nitrogen_Demand)
    
    
    # =============================================================================
    # Canopy processes for updating Nitrogen dynamic based on the nitrogen absorption from soil
    # =============================================================================
 
    Canopy_object.Calculate_Nitrogen_Dynamics(Soil_object.Nitrogen_uptake)
    
    Canopy_object.Calculate_Seed_Number_Rate(Soil_object.Nitrogen_uptake,Root_object.Root_carbon_loss_rate_senescence)
    
    Canopy_object.Calculate_Nitrogen_Accumulation_Rate(Soil_object.Nitrogen_uptake,Root_object.Root_nitrogen_loss_rate_senescence,Inputs.Standard_SeedNitrogen_Conc)
    
    ########################
    Canopy_object.Calculate_Leaf_Area_ChangeRate(Canopy_object.Carbon_determined_LAI,Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Bottom_Exponential_with_Depth'],Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff'])

    Root_object.Calculate_Rooting_Depth(Canopy_object.RootWeight_Rate, Canopy_object.LiveRoot_Dry_Weight, Canopy_object.DeadRoot_Dry_Weight)

    Canopy_object.Calculate_Respiration_Rates(Soil_object.Nitrogen_uptake, Root_object.Root_carbon_loss_rate_senescence)
    
    Canopy_object.Calculate_Carbon_Nitrogen_Returns(Root_object.Root_carbon_loss_rate_senescence, Root_object.Root_nitrogen_loss_rate_senescence,Soil_object.average_soil_temperature)
        

    # =============================================================================
    # updating Soil water balance, carbon, and Nitrogen
    # =============================================================================
    
    
    
    # =============================================================================
    # updating state variables 
    # =============================================================================
    Canopy_object.Update_State_Variables(Root_object.Root_nitrogen_loss_rate_senescence)
    
    Root_object.Update_State_Variables()
    Soil_object.Update_State_Variables()

        
    