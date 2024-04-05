
import numpy as np
import pandas as pd
import Root
import Soil
import Leaf
import Canopy
import AtmoWeather


weather_meteo = AtmoWeather.WeatherMeteo(f"{weatherfilename}.xlsx", lat, insp, sttime)
weather_data = weather_meteo.process_data()  # This contains weather data for all timesteps
Canopy_object=Canopy.Canopy(BTemp_Phen, OTemp_Phen, CTemp_Phen, TempCurve_Res, CropType_Photoperiod, StartPhotoperiod_Phase, EndPhotoperiod_Phase, Photoperiod_Sensitivity, MinThermal_Day_Veg, MinThermal_Day_Rep,
             CarbonAlloc_Shoot, NitrogenAlloc_Shoot, Plant_Density, Seed_Weight, 
             Crop_TypeDet, EndSeedNum_DetPeriod, SeedN_RemobFract, SLA_Const, MinSLN_Photosyn, MinRootN_Conc,
             CarbonCost_NFix, MaxN_Uptake, MinStemN_Conc,
             IniLeafN_Conc, MaxPlant_Height,
             legume,
             MaxStemGrowth_DS, MaxSeedGrowth_DS, StemDW_Height, Model_TimeStep,
             IniNConc_SeedFill, FinalNConc_SeedFill, tm)

Root_object=Root.Root(Canopy_object.cfv, root_to_shoot_ratio, max_root_depth, Canopy_object.Nitrogen_Root, Canopy_object.CarbonReserve_Root, Canopy_object.min_N_concentration_root, Canopy_object.Model_TimeStep, Soil_Depth_1)
Soil_object=Soil.Soil(Residual_Water_Content, Saturated_Water_Content, temperature_change_constant, lodging_condition,
             Soil_Depth_1, initial_soil_temp, initial_decomposable_plant_material, field_capacity_water_content, clay_percentage, residual_ammonium,
             residual_nitrate, dpm_decomposition_rate, rpm_decomposition_rate,
             biomass_incorporation_rate, humification_rate, total_organic_carbon, biochar_carbon_content,
             biochar_conversion_fraction, parameter_adjustment_multiplier, initial_ammonium_content, initial_nitrate_content,
             nitrogen_stress_water_index, water_supply_switch, daily_water_input, ammonium_nitrogen_input_rate,
             nitrate_nitrogen_input_rate, plant_material_dpm_rpm_ratio, soil_resistance_to_evaporation,
             sand_percentage)
Leaf_object=Leaf.Leaf(Veg_C_Fraction, Spec_Leaf_Area, LAI_ini, Leaf_Blade_Angle, Leaf_Width, Tot_Leaf_N, Min_Spec_Leaf_N, Pathway_C3C4, Ambient_CO2, Activation_Energy_JMAX, VCMAX_LeafN_Slope, JMAX_LeafN_Slope, Photosynthetic_Light_Response_Factor)
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
    dap=doy-planting_doy+1
    


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
    
    Soil_object.Calculate_Soil_Potential_Evaporation(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Sin_Beam_Exposure, tmax, tmin, Vapour_Pressure,Wind_Speed,rss,bld,
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
    
    


    Soil_object.Update_Evaporation_if_WaterStress(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,Soil_object.dwsup,Soil_Depth_1,Soil_object.Root_Depth,                                                  
                                                       Leaf_object.hourly_transpiration_SU, Leaf_object.hourly_transpiration_SH,
                                                       Soil_object.hourly_Soil_Evap,)
    
    Soil_object.Calculate_Soil_Temperature(tmin,tmax)


    
    # =============================================================================
    # Thermal units and initializing Canopy Carbon and Nitrogen Accumulation for Biomass formation
    # =============================================================================
    Canopy_object.calculate_thermal_units(Canopy_object.ds, tmax, tmin, Day_Length, BTemp_Phen, OTemp_Phen, CTemp_Phen, TempCurve_Res)
    Canopy_object.calculate_developement_rate(Canopy_object.ds, CropType_Photoperiod, Day_Length, StartPhotoperiod_Phase, EndPhotoperiod_Phase, Photoperiod_Sensitivity, MinThermal_Day_Veg, MinThermal_Day_Rep, Canopy_object.DailyThermalUnit)
       
    Canopy_object.Initialize_Biomass_Formation()
    
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
    
    Canopy_object.Calculate_Senescence(Leaf_object.laic, Leaf_object.leaf_area_output['N_determined_LAI'], Model_TimeStep)
    
    
    # =============================================================================
    # Canopy processes for Carbon flow and partitioning
    # =============================================================================
    
    Canopy_object.Calculate_Carbon_Flow()

    Canopy_object.Calculate_Carbon_Flow_Seed_Filling()

    Canopy_object.Calculate_Carbon_Flow_Stem_Growth()

    Root_object.Calculate_Root_Senescence()

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

        
    