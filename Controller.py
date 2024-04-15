
import numpy as np
import pandas as pd
import Root
import Soil
import Leaf
import Canopy
import AtmoWeather
from Input_Setup import ImportData


User_Inputs = r'User_Inputs.xlsx'
Default_Inputs = r'Default_Inputs.xlsx'
weatherfilename='weatherfile.xlsx'

# Load crop data from both files
User_Inputs = ImportData(User_Inputs)
Default_Inputs = ImportData(Default_Inputs)
Inputs=User_Inputs.append(Default_Inputs)

weather_meteo = AtmoWeather.WeatherMeteo(weatherfilename, Inputs.latitude, Inputs.Simulation_Start_DOY)
weather_data = weather_meteo.process_data()  # This contains weather data for all timesteps
Canopy_object=Canopy.Canopy(Inputs.BTemp_Phen, Inputs.OTemp_Phen, Inputs.CTemp_Phen, Inputs.TempCurve_Res, 
                            Inputs.CropType_Photoperiod, Inputs.StartPhotoperiod_Phase, Inputs.EndPhotoperiod_Phase, 
                            Inputs.Photoperiod_Sensitivity, Inputs.MinThermal_Day_Veg, Inputs.MinThermal_Day_Rep,
             Inputs.CarbonAlloc_Shoot, Inputs.NitrogenAlloc_Shoot, Inputs.Plant_Density, Inputs.Seed_Weight, 
             Inputs.Crop_TypeDet, Inputs.EndSeedNum_DetPeriod, Inputs.SeedN_RemobFract, Inputs.SLA_Const, Inputs.MinSLN_Photosyn, Inputs.MinRootN_Conc,
             Inputs.CarbonCost_NFix, Inputs.MaxN_Uptake, Inputs.MinStemN_Conc,
             Inputs.IniLeafN_Conc, Inputs.MaxPlant_Height,
             Inputs.legume,
             Inputs.MaxStemGrowth_DS, Inputs.MaxSeedGrowth_DS, Inputs.StemDW_Height, Inputs.Model_TimeStep)

Root_object=Root.Root(Inputs.root_to_shoot_ratio, Inputs.max_root_depth,  Inputs.Model_TimeStep, Inputs.Soil_Depth_1)
Soil_object=Soil.Soil( Inputs.Residual_Water_Content, Inputs.Saturated_Water_Content, Inputs.temperature_change_constant,  Inputs.field_capacity_water_content, 
             Inputs.Soil_Depth_1,Inputs.clay_percentage,Inputs.sand_percentage,
              Inputs.initial_soil_temp,  Inputs.daily_water_input, Inputs.soil_resistance_to_evaporation,
             Inputs.fraction_soil_greater_than_2mm,Inputs.soil_bulk_density,Inputs.organic_N_percentage,Inputs.fraction_N_for_mineralization,
             Inputs.nitrate_concentration_ppm,Inputs.ammonium_concentration_ppm,Inputs.water,
             Inputs.Fertilizer_applications_count,Inputs.Fertilizer_applications_amount,Inputs.Fertilizer_applications_DAP,Inputs.Fraction_volatilization,
             Inputs.Days_after_planting)
Leaf_object=Leaf.Leaf( Inputs.Spec_Leaf_Area, Inputs.LAI_ini, Inputs.Leaf_Blade_Angle, Inputs.Leaf_Width,  Inputs.Min_Spec_Leaf_N, Inputs.Pathway_C3C4, Inputs.Ambient_CO2,
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
    tmax = day_data['tmax']
    tmin = day_data['tmin']
    Vapour_Pressure = day_data['Vapour_Pressure']
    Wind_Speed = day_data['Wind_Speed']
    rain = day_data['Rain']
    doy=day_data['Doy']
    Days_after_planting=doy-Inputs.planting_doy+1
    


    # The Leaf Area Update could be part of canopy as well
    Leaf_object.Update_Leaf_Area()
    Leaf_object.Update_Specific_Leaf_N()

    # =============================================================================
    # Instantiating Sunlit and shaded leaves for calculating Potential PHOTOSYNTHESIS AND TRANSPIRATION
    # =============================================================================
    
    Leaf_sunlit_object = Leaf.Leaf_sunlit(Leaf_object)
    Leaf_shaded_object = Leaf.Leaf_Shaded(Leaf_object)



    # Call the Leaf method with the current timestep's weather data --> for sunlit parts
    Leaf_sunlit_object.Calculate_leaf_temp(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height)
    Leaf_sunlit_object.Calculate_Potential_Photosynthesis(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed)
    Leaf_sunlit_object.Calculate_Potential_transpiration(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,Canopy_object.Plant_Height)
    
    
    # Call the Leaf method with the current timestep's weather data --> for shaded parts
    Leaf_shaded_object.Calculate_leaf_temp(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height)
    Leaf_shaded_object.Calculate_Potential_Photosynthesis(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed)
    Leaf_shaded_object.Calculate_Potential_transpiration(Solar_Constant,Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,Canopy_object.Plant_Height)
    
    
    # =============================================================================
    # Updating Potential Canopy  PHOTOSYNTHESIS AND TRANSPIRATION
    # =============================================================================
    Canopy_object.Update_Canopy_Temp("P",Leaf_object.hourly_Air_SU_leaf_T_diff,Leaf_object.hourly_Air_SH_leaf_T_diff,
                            Leaf_object.hourly_SU_leaf_T, Leaf_object.hourly_SH_leaf_T)
    Canopy_object.Update_Canopy_Photosyn("P",Leaf_object.hourly_photosyn_SU,Leaf_object.hourly_photosyn_SH,Day_Length)
    Canopy_object.Update_Canopy_Transpiration("P",Leaf_object.hourly_transpiration_SU,Leaf_object.hourly_transpiration_SH,Day_Length)
    
    


    # =============================================================================
    # Calculating Potential Evaporation from soil
    # =============================================================================

    Soil_object.Calculate_Soil_Water_Content()
    
    Soil_object.Calculate_Soil_Potential_Evaporation(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,
                                                     Soil_object.soil_resistance_to_evaporation,Leaf_object.Leaf_Blade_Angle,
                                                     Leaf_object.leaf_area_output['TLAI'],Leaf_object.leaf_area_output['KW'],
                            Leaf_object.hourly_transpiration_SH,Leaf_object.hourly_transpiration_SU)




    # =============================================================================
    # Updating  Photosynthesis, Tranpisraion, and Evaporation if Water Stress ocuurs
    # =============================================================================


    Leaf_sunlit_object.Update_LeafTemp_Photosynthesis_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,
                           Soil_object.dwsup,Soil_object.Soil_Depth_1,Soil_object.Root_Depth, 
                           Leaf_object.hourly_SU_leaf_T, Leaf_object.hourly_SH_leaf_T,
                           Soil_object.hourly_Soil_Evap)
    Leaf_shaded_object.Update_LeafTemp_Photosynthesis_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure, Wind_Speed,Canopy_object.Plant_Height,
                           Soil_object.dwsup,Soil_object.Soil_Depth_1,Soil_object.Root_Depth, 
                           Leaf_object.hourly_SU_leaf_T, Leaf_object.hourly_SH_leaf_T,
                           Soil_object.hourly_Soil_Evap)
    

    Canopy_object.Update_Canopy_Temp("A",Leaf_object.hourly_Actual_Air_SU_leaf_T_diff,Leaf_object.hourly_Actual_Air_SH_leaf_T_diff,
                            Leaf_object.hourly_Actual_SU_leaf_T, Leaf_object.hourly_Actual_SH_leaf_T)
    Canopy_object.Update_Canopy_Photosyn("A",Leaf_object.Hourly_Actual_Photosynthesis_SU,Leaf_object.Hourly_Actual_Photosynthesis_SH,Day_Length)
    Canopy_object.Update_Canopy_Transpiration("A",Leaf_object.Hourly_Actual_Transpiration_SU,Leaf_object.Hourly_Actual_Transpiration_SH,Day_Length)
    
    


    Soil_object.Update_Evaporation_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,Soil_object.dwsup,Soil_object.Soil_Depth_1,Soil_object.Root_Depth,                                                  
                                                       Leaf_object.hourly_transpiration_SU, Leaf_object.hourly_transpiration_SH,
                                                       Soil_object.hourly_Soil_Evap,)
    
    Soil_object.Calculate_Soil_Temperature(tmin,tmax)


    
    # =============================================================================
    # Thermal units and initializing Canopy Carbon and Nitrogen Accumulation for Biomass formation
    # =============================================================================
    Canopy_object.calculate_thermal_units(Canopy_object.ds, tmax, tmin, Day_Length, Canopy_object.BTemp_Phen, Canopy_object.OTemp_Phen, Canopy_object.CTemp_Phen, Canopy_object.TempCurve_Res)
    Canopy_object.calculate_developement_rate(Canopy_object.ds, Canopy_object.CropType_Photoperiod, Day_Length, Canopy_object.StartPhotoperiod_Phase, Canopy_object.EndPhotoperiod_Phase,
                                              Canopy_object.Photoperiod_Sensitivity, Canopy_object.MinThermal_Day_Veg, Canopy_object.MinThermal_Day_Rep, Canopy_object.DailyThermalUnit)
       
    Canopy_object.Initialize_Biomass_Formation(Canopy_object.LiveRoot_Carbon)
    
    Canopy_object.Initialize_Nitrogen_Accumulation()
    
    Canopy_object.Initialize_Carbon_Accumulation()
    
    
    
    # =============================================================================
    # Canopy processes for Respiration, Nitrogen demand and partitioning
    # =============================================================================
    
    Canopy_object.Calculate_Respiration()
    
    Canopy_object.Calculate_Nitrogen_Fixation()
    
    Canopy_object.Calculate_Photo_Assimilates()
    
    Canopy_object.Calculate_Crop_Nitrogen_Demand(Leaf_object.specific_leaf_n_output['Spec_Leaf_N'])
    
    Canopy_object.Calculate_Nitrogen_Partitioning(Leaf_object.specific_leaf_n_output['Spec_Leaf_N_Top_Increment'])
    
    Canopy_object.Calculate_Seed_Properties(Determinate, End_SeedFill_Indeterminate, Nitrogen_Removal_EarlyFlowering, Nitrogen_Removal_EndSeedFill, Fraction_SeedNitrogen_Remobilizable, Expected_SeedNitrogen_Concentration, Seed_Weight)
    
    Canopy_object.Calculate_Senescence(Leaf_object.laic, Leaf_object.leaf_area_output['N_determined_LAI'], Canopy_object.Model_TimeStep)
    
    
    # =============================================================================
    # Canopy processes for Carbon flow and partitioning
    # =============================================================================
    
    Canopy_object.Calculate_Carbon_Flow()

    Canopy_object.Calculate_Carbon_Flow_Seed_Filling()

    Canopy_object.Calculate_Carbon_Flow_Stem_Growth()

    Root_object.Calculate_Root_Senescence(Canopy_object.carbon_structural_roots,Canopy_object.CarbonFrac_Veg, Canopy_object.MinRootN_Conc, Canopy_object.ReserveRoot_Carbon, Canopy_object.Nitrogen_Root)


    Canopy_object.Calculate_Carbon_Partitioning(Leaf_object.leaf_area_output['N_determined_LAI'],Leaf_object.laic, Root_object.csrtn)
    
    Canopy_object.Update_Carbon_Reserve_Pools()
    
    Canopy_object.Calculate_Carbon_Production_Rate()
    
    Canopy_object.Calculate_Biomass_Change_Rate()
    
    Canopy_object.Calculate_Carbon_rates()
    
    
    # =============================================================================
    # Soil processes for caculating the state of carbon and Nitrogen for N uptake
    # =============================================================================
    
    Soil_object.Organic_Carbon_Composition()
    
    Soil_object.Organic_Nitrogen_Composition(rain)

    Soil_object.Calculate_Nitrogen_Uptake()
    
    
    # =============================================================================
    # Canopy processes for updating Nitrogen dynamic based on the nitrogen absorption from soil
    # =============================================================================
 
    Canopy_object.Calculate_Nitrogen_Dynamics()
    
    Canopy_object.Calculate_Seed_Number_Rate()
    
    Canopy_object.Calculate_Nitrogen_Accumulation_Rate()
    
    Canopy_object.Calculate_Leaf_Area_ChangeRate(Leaf_object.specific_leaf_n_output['Spec_Leaf_N_Base_Calc'],Leaf_object.leaf_area_output['Leaf_Nitro_Ext_Coeff'])

    Root_object.Calculate_Rooting_Depth()

    Canopy_object.Calculate_Respiration_Rates()
    
    Canopy_object.Calculate_Carbon_Nitrogen_Returns()
    
    Canopy_object.Check_Carbon_Nitrogen_Balance()
    

    # =============================================================================
    # updating Soil water balance, carbon, and Nitrogen
    # =============================================================================
    Soil_object.Calculate_Water_Balance()
    
    Soil_object.Calculate_Soil_Organic_Nitrogen_ChangeRate(dap,Root_object.Root_Depth)
    
    
    Soil_object.Calculate_Soil_Organic_Carbon_ChangeRate()
    
    
    # =============================================================================
    # updating state variables 
    # =============================================================================
    Canopy_object.Update_State_Variables()
    Leaf_object.Update_State_Variables(Canopy_object.Rate_LeafAreaIndex)
    Root_object.Update_State_Variables()
    Soil_object.Update_State_Variables()

        
    