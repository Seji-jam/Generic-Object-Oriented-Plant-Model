import numpy as np
import math
import Leaf

def SWITCH_FUN(x, y1, y2):
    return y1 if x < 0 else y2
wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
class Soil:
    def __init__(self, Residual_Water_Content, Saturated_Water_Content, temperature_change_constant,  field_capacity_water_content, 
                 Soil_Depth_1,clay_percentage,sand_percentage,
                  initial_soil_temp,  daily_water_input, soil_resistance_to_evaporation,
                 fraction_soil_greater_than_2mm,soil_bulk_density,organic_N_percentage,fraction_N_for_mineralization,
                 nitrate_concentration_ppm,ammonium_concentration_ppm,water,
                 Fertilizer_applications_count,Fertilizer_applications_amount,Fertilizer_applications_DAP,Fraction_volatilization,
                 Days_after_planting,):
        
        self.Residual_Water_Content = Residual_Water_Content  # Residual water content
        self.Saturated_Water_Content = Saturated_Water_Content  # Saturated water content
        self.temperature_change_constant = temperature_change_constant
        # self.lodging_condition = lodging_condition  # Lodging condition
        self.Soil_Depth_1 = Soil_Depth_1  # Soil bulk density
        self.initial_soil_temp = initial_soil_temp  # Initial soil temperature
        # self.initial_decomposable_plant_material = initial_decomposable_plant_material  # Initial decomposable plant material in soil
        self.field_capacity_water_content = field_capacity_water_content  # Water content at field capacity
        self.clay_percentage = clay_percentage  # Clay content in soil
        # self.residual_ammonium = residual_ammonium  # Residual ammonium-N in the soil
        # self.residual_nitrate = residual_nitrate  # Residual nitrate-N in the soil
        # self.dpm_decomposition_rate = dpm_decomposition_rate  # Decomposition rate of plant material (dpm rate)
        # self.rpm_decomposition_rate = rpm_decomposition_rate  # Decomposition rate of resistant plant material (rpm rate)
        # self.biomass_incorporation_rate = biomass_incorporation_rate  # Biomass incorporation rate into the soil organic matter pool
        # self.humification_rate = humification_rate  # Humification rate of organic matter
        # self.total_organic_carbon = total_organic_carbon  # Total organic carbon in the soil
        # self.biochar_carbon_content = biochar_carbon_content  # Biochar carbon content in the soil
        # self.biochar_conversion_fraction = biochar_conversion_fraction  # Fraction of biomass converted to biochar
        # self.parameter_adjustment_multiplier = parameter_adjustment_multiplier  # Multiplier factor for adjusting model parameters
        # self.initial_ammonium_content = initial_ammonium_content  # Initial Ammonium nitrogen content in soil
        # self.initial_nitrate_content = initial_nitrate_content  # Initial Nitrate nitrogen content in soil
        # self.nitrogen_stress_water_index = nitrogen_stress_water_index  # Nitrogen stress water index
        # self.water_supply_switch = water_supply_switch  # Switch variable for water supply to crop
        self.daily_water_input = daily_water_input  # Daily water input (precipitation + irrigation)
        # self.ammonium_nitrogen_input_rate = ammonium_nitrogen_input_rate  # Ammonium nitrogen input rate
        # self.nitrate_nitrogen_input_rate = nitrate_nitrogen_input_rate  # Nitrate nitrogen input rate
        # self.plant_material_dpm_rpm_ratio = plant_material_dpm_rpm_ratio  # Ratio dpm/rpm of added plant material
        self.soil_resistance_to_evaporation = soil_resistance_to_evaporation  # Soil resistance to evaporation
        self.sand_percentage = sand_percentage  # Percent of Sand in soil
        self.fraction_soil_greater_than_2mm = fraction_soil_greater_than_2mm/100  # Fraction soil > 2 mm (soil coarse fraction)
        self.soil_bulk_density = soil_bulk_density  # Soil bulk density (g.cm-3)
        self.organic_N_percentage = organic_N_percentage  # Organic N (%)
        self.fraction_N_for_mineralization = fraction_N_for_mineralization  # Fraction Organic N available for mineralization
        self.nitrate_concentration_ppm = nitrate_concentration_ppm  # Nitrate concentration (ppm = mg.kg−1)
        self.ammonium_concentration_ppm = ammonium_concentration_ppm  # Ammonium concentration (ppm = mg.kg−1)
        self.Soil_Depth_1 = Soil_Depth_1  # Depth, assuming necessary for calculation
        self.water = water  # Water, assuming necessary for calculation
        self.Fertilizer_applications_count = Fertilizer_applications_count  
        self.Fertilizer_applications_amount = Fertilizer_applications_amount  
        self.Fraction_volatilization = Fraction_volatilization  
       
        self.Fertilizer_applications_DAP=list(Fertilizer_applications_DAP)
        self.Fertilizer_applications_amount=list(Fertilizer_applications_amount)
        self.Fraction_volatilization=list(Fraction_volatilization)
       
       
        soil_mass = self.Soil_Depth_1 * soil_bulk_density * (1 - fraction_soil_greater_than_2mm) * 1000  # g.m-2
        organic_N_mass = organic_N_percentage * 0.01 * soil_mass  # g.m-2
        mineralizable_organic_N = organic_N_mass * fraction_N_for_mineralization  # gN.m-2 org. N avail. for miner.
        nitrate_mass = nitrate_concentration_ppm * (14 / 62) * 0.000001 * soil_mass  # from ppm to gN.m-2
        ammonium_mass = ammonium_concentration_ppm * (14 / 18) * 0.000001 * soil_mass  # from ppm to gN.m-2
        self.soluble_N = nitrate_mass + ammonium_mass  # gN.m-2
        self.N_concentration = self.soluble_N / (water * 1000)  # gN.g-1 H2O
        



        
        #initializing intermediate parameters
        self.decomposition_rate_decomposable_pm = 0  # Rate of decomposition for decomposable plant material
        self.decomposition_rate_resistant_pm = 0  # Rate of decomposition for resistant plant material
        self.decomposition_rate_decomposable_pn = 0  # Rate of decomposition for decomposable plant nitrogen
        self.decomposition_rate_resistant_pn = 0  # Rate of decomposition for resistant plant nitrogen
        self.decomposition_rate_biomass = 0  # Rate of decomposition for biomass
        self.decomposition_rate_humus = 0  # Rate of decomposition for humus
        self.mineralized_nitrogen = 0  # Total mineralized nitrogen
        self.mineralized_nitrogen_upper_layer = 0  # Mineralized nitrogen in the upper soil layer
        self.mineralized_nitrogen_lower_layer = 0  # Mineralized nitrogen in the lower soil layer
        
        self.nitrification_upper_layer = 0  # Nitrification rate in the upper soil layer
        self.nitrification_lower_layer = 0  # Nitrification rate in the lower soil layer
        self.mineralized_ammonium_upper_layer = 0  # Mineralized ammonium in the upper soil layer
        self.mineralized_ammonium_lower_layer = 0  # Mineralized ammonium in the lower soil layer
        # Note: Duplicate self.minaul removed
        self.water_stress_factor = 0  # Factor indicating water stress
        self.nitrogen_uptake = 0  # Nitrogen uptake rate


        self.hourly_Soil_Evap=[]
        self.hourly_Soil_Rad=[]
        

        # Irrigation and fertilization can be input in the weather file as well! 
        # Irrigation schedule: Days and corresponding irrigation amount (currently set to 0)
        self.irrigation_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), 
                                    (55.0, 0), (65.0, 0), (75.0, 0), (85.0, 0), (95.0, 0)]
        
        # # Fertilization schedules for ammonium  and nitrate nitrogen (currently set to 0)
        # self.ammonium_nitrogen_application_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), 
        #                                                (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]
        # self.nitrate_nitrogen_application_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), 
        #                                               (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]

        self.deniul = 0  # Denitrification in the upper layer
        self.denill = 0
        self.nuptn = 0  # Nitrogen uptake in the nitrogen form
        self.nupta= 0



        
        # Initial conditions and parameter calculations
        root_depth_ini = max(2.0, Soil_Depth_1)
        # biochar_initial_content = biochar_conversion_fraction * total_organic_carbon
        water_content_initial = field_capacity_water_content * parameter_adjustment_multiplier
        water_upper_soil_initial = 10.0 * (water_content_initial - Residual_Water_Content) * root_depth_ini
        water_lower_soil_initial = 10.0 * (water_content_initial - Residual_Water_Content) * (150.0 - root_depth_ini)
        # resistant_plant_material_initial = total_organic_carbon - biochar_initial_content - initial_decomposable_plant_material
        # humus_initial_content = biochar_initial_content - biochar_initial_content
        # decomposable_plant_nitrogen_initial = 1.0 / 40.0 * initial_decomposable_plant_material
        # ammonium_upper_layer_initial = (1.0 - np.exp(-0.065 * root_depth_ini)) * initial_ammonium_content + root_depth_ini / 150.0 * residual_ammonium
        # ammonium_lower_layer_initial = np.exp(-0.065 * root_depth_ini) * initial_ammonium_content + (1.0 - root_depth_ini / 150.0) * residual_ammonium
        # nitrate_uptake_limit_initial = (1.0 - np.exp(-0.065 * root_depth_ini)) * initial_nitrate_content + root_depth_ini / 150.0 * residual_nitrate
        # nitrate_loss_limit_initial = np.exp(-0.065 * root_depth_ini) * initial_nitrate_content + (1.0 - root_depth_ini / 150.0) * residual_nitrate
        # resistant_plant_nitrogen_initial = 1.0 / 100.0 * resistant_plant_material_initial

    
        self.average_soil_temperature = 0
        self.rate_of_change_soil_temp = 0
    
        self.soil_water_content_upper_limit = 0
        self.soil_water_content_lower_limit = 0
        self.nitrate_nitrogen_supply = 0
        self.ammonium_fertilizer_application = 0
        self.nitrate_fertilizer_application = 0

        self.ammonia_volatilization = 0
        # State variables
        self.total_nitrogen_uptake = 0
        # self.decomposable_plant_material = initial_decomposable_plant_material  
        # self.resistant_plant_material = resistant_plant_material_initial  
        # self.biomass_incorporated = biochar_initial_content  
        # self.humus_content = humus_initial_content  
        # self.decomposable_plant_nitrogen = decomposable_plant_nitrogen_initial  
        # self.resistant_plant_nitrogen = resistant_plant_nitrogen_initial   
        self.water_upper_soil = water_upper_soil_initial 
        self.water_lower_soil = water_lower_soil_initial 
        self.Root_Depth = root_depth_ini  
        self.initial_soil_temperature = initial_soil_temp 
        # self.ammonium_upper_layer = ammonium_upper_layer_initial
        # self.nitrate_upper_layer = nitrate_uptake_limit_initial
        # self.nitrate_lower_layer = nitrate_loss_limit_initial
        # self.ammonium_lower_layer = ammonium_lower_layer_initial
        self.total_nitrogen_leached = 0
        self.soil_ammonium_volatilization_rate = 0

        self.potential_evap_daily = 0
        self.potential_SoilRad_daily = 0
        self.Day_Air_Soil_temp_dif  = 0
        self.Actual_evap_daily =0
        self.hourly_rbhs=[]
        self.hourly_rts=[]


     


        
    
        def Calculate_Soil_Water_Content(self):
            # Calculate the water content in the upper soil layer, adjusting for root depth
            water_content_upper_layer = (self.water_content_upper_soil_layer + self.Residual_Water_Content * 10.0 * self.Root_Depth) / 10.0 / self.Root_Depth
            
            # Calculate the water content of the lower soil layer, considering the soil's root depth and making sure it doesn't exceed the saturated water content
            water_content_lower_layer = min(self.Saturated_Water_Content, (self.water_content_lower_soil_layer + self.Residual_Water_Content * 10.0 * (150.0 - self.Root_Depth)) / 10.0 / (150.0 - self.Root_Depth))
            
            # Determine the daily water supply available for evapotranspiration, with a safeguard for a minimum value
            daily_water_supply_for_et = self.switch_fun(-1, self.daily_water_input, max(0.1, self.water_content_upper_soil_layer / self.temperature_change_constant + 0.1))
            
            # Update the class attributes with the calculated values
            self.soil_water_content_upper_layer = water_content_upper_layer
            self.soil_water_content_lower_layer = water_content_lower_layer
            self.daily_water_supply_for_et = daily_water_supply_for_et
    





    
        
        def Calculate_Soil_Potential_Evaporation(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature, Min_Temperature, Vapour_Pressure_Deficit, Wind_Speed, Leaf_Blade_Angle, soil_resistance_to_evaporation, Total_Leaf_Area_Index, Light_Extinction_Coefficient, hourly_transpiration_Shaded, hourly_transpiration_Sunlit):
                
                # Convert daily climate data to hourly data
                hourly_climate_data, gaussian_weights = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature, Min_Temperature, Wind_Speed)
                
                hourly_Soil_Evap = []
                hourly_Air_Soil_temp_dif = []
                hourly_Soil_Rad = []
                hourly_rbhs = []
                hourly_rts = []
        
                for (hourly_Solar_Constant, hourly_temperature, hourly_sin_beam, hourly_diffuse_ratio, hourly_wind_speed), shaded_transpiration, sunlit_transpiration in zip(hourly_climate_data, hourly_transpiration_Shaded, hourly_transpiration_Sunlit):
                    
                    Incoming_PAR = 0.5 * Solar_Radiation  # Partition of incoming solar radiation to PAR
                    Incoming_NIR = 0.5 * Solar_Radiation  # Partition of incoming solar radiation to NIR
                    
                    # Atmospheric transmissivity based on incoming PAR and solar geometry
                    atmospheric_transmissivity = Incoming_PAR / (0.5 * hourly_Solar_Constant * hourly_sin_beam)
        
                    # Saturation vapor pressure and vapor pressure deficit calculations
                    saturation_vapor_pressure = 0.611 * math.exp(17.4 * hourly_temperature / (hourly_temperature + 239.))
                    vapor_pressure_deficit = max(0., saturation_vapor_pressure - Vapour_Pressure_Deficit)
                    
                    
                    slope_vapor_pressure_temperature = 4158.6 * saturation_vapor_pressure / (hourly_temperature + 239.) ** 2
    
                    # Turbulence resistance for soil using logarithmic wind profile
                    turbulence_resistance_soil = 0.74 * (np.log(56.)) ** 2 / (0.4 ** 2 * hourly_wind_speed)
                    
                    # Boundary layer resistance for soil considering wind speed and leaf area
                    boundary_layer_resistance_soil = 172. * np.sqrt(0.05 / max(0.1, hourly_wind_speed * np.exp(-Light_Extinction_Coefficient * Total_Leaf_Area_Index)))
                    
                
                
                
                    water_boundary_layer_resistance_soil = 0.93 * boundary_layer_resistance_soil
                    
                    # Calculation of the fractional diffuse light (frdf) based on atmospheric transmissivity
                    if atmospheric_transmissivity < 0.22:
                        fractional_diffuse_light = 1
                    elif 0.22 < atmospheric_transmissivity <= 0.35:
                        fractional_diffuse_light = 1 - 6.4 * (atmospheric_transmissivity - 0.22) ** 2
                    else:
                        fractional_diffuse_light = 1.47 - 1.66 * atmospheric_transmissivity
        
                    # Ensuring a minimum threshold for the fractional diffuse light
                    fractional_diffuse_light = max(fractional_diffuse_light, 0.15 + 0.85 * (1 - np.exp(-0.1 / hourly_sin_beam)))
            
                    # Incoming diffuse PAR (PARDF) and direct PAR (PARDR) based on the calculated fractional diffuse light
                    PAR_Diffuse = Incoming_PAR * fractional_diffuse_light
                    PAR_Direct = Incoming_PAR - PAR_Diffuse
        
                    # Extinction and RefLection coefficients
                    # Convert leaf blade angle from degrees to radians for calculations
                    Leaf_Blade_Angle_Radians = Leaf_Blade_Angle * np.pi / 180
                    
                    # Calculate the direct beam extinction coefficient for PAR using leaf blade angle
                    Direct_Beam_Extinction_Coefficient_PAR = Leaf.KDR_Coeff(hourly_sin_beam, Leaf_Blade_Angle_Radians)
                    
                    # Leaf scattering coefficients for PAR and NIR
                    Scattering_Coefficient_PAR = 0.2  # Scattering coefficient for Photosynthetically Active Radiation
                    Scattering_Coefficient_NIR = 0.8  # Scattering coefficient for Near-Infrared Radiation
                    
                    # Calculate diffuse extinction coefficients for PAR and NIR using total leaf area index and leaf blade angle
                    Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
                    Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
                    
                    # Calculate beam and canopy reflection coefficients for PAR and NIR
                    Beam_Plus_Canopy_Reflection_Coefficient_PAR, Canopy_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient_PAR)
                    Beam_Plus_Canopy_Reflection_Coefficient_NIR, Canopy_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient_PAR)
        
                    
                    # Incoming diffuse NIR (NIRDF) and direct NIR (NIRDR)
                    NIR_Diffuse = Incoming_NIR * fractional_diffuse_light
                    NIR_Direct = Incoming_NIR - NIR_Diffuse
                    
                    # Absorbed total radiation (PAR + NIR) by soil
                    soil_par_reflection = 0.1  # Soil PAR reflection coefficient
                    # Soil NIR reflection coefficient, varying with soil water content (wcul)
                    soil_nir_reflection = self.switch_function(self.soil_water_content_upper_layer - 0.5, 0.52 - 0.68 * self.soil_water_content_upper_layer, 0.18)
                    absorbed_total_radiation_by_soil = (1 - soil_par_reflection) * (PAR_Direct * np.exp(-Beam_Plus_Canopy_Reflection_Coefficient_PAR * Total_Leaf_Area_Index) + PAR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_PAR * Total_Leaf_Area_Index)) + \
                                                       (1 - soil_nir_reflection) * (NIR_Direct * np.exp(-Beam_Plus_Canopy_Reflection_Coefficient_NIR * Total_Leaf_Area_Index) + NIR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_NIR * Total_Leaf_Area_Index))
                    # Calculate potential evaporation and net radiation using Penman-Monteith equation
                    potential_evaporation, net_radiation = Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, hourly_temperature, vapor_pressure_deficit, slope_vapor_pressure_temperature, vapor_pressure_deficit)
        
                    # Proportional transpiration
                    proportional_transpiration = (shaded_transpiration + sunlit_transpiration) * self.sd1 / self.Root_Depth
        
        
                   # Daytime course of water supply
                    water_supply = self.daily_water_supply * (hourly_sin_beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
                    proportional_water_supply = water_supply * self.sd1 / self.Root_Depth
                    
                    potential_evaporation_soil = max(0, potential_evaporation)
                    actual_evaporation_soil = min(potential_evaporation_soil, potential_evaporation_soil / (proportional_transpiration + potential_evaporation_soil) * proportional_water_supply)
                    
                    Latent_Heat_of_Water_Vaporization = 2.4E6  # Latent heat of vaporization (J/kg)
                    Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity of air (J/m^3/°C)
                    
                    # Temperature difference driven by soil evaporation and net radiation
                    temperature_difference_soil = self.Limit_Function(-25., 25., (net_radiation - Latent_Heat_of_Water_Vaporization * actual_evaporation_soil) * (boundary_layer_resistance_soil + turbulence_resistance_soil) / Volumetric_Heat_Capacity_Air) 
                    
                    average_soil_temperature = hourly_temperature + temperature_difference_soil
                    
                    # Recalculate potential evaporation with updated soil temperature
                    saturation_vapor_pressure_soil = 0.611 * math.exp(17.4 * average_soil_temperature / (average_soil_temperature + 239.))
                    slope_vapor_pressure_curve_soil = (saturation_vapor_pressure_soil - saturation_vapor_pressure) / self.Avoid_Zero_Division(temperature_difference_soil)
                    potential_evaporation_second_estimate, net_radiation_second_estimate = Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, average_soil_temperature, Vapour_Pressure_Deficit, slope_vapor_pressure_curve_soil, vapor_pressure_deficit)
                    potential_soil_evaporation = max(0, potential_evaporation_second_estimate)
                    print(potential_soil_evaporation)
                    hourly_Soil_Evap.append(potential_soil_evaporation)
                    hourly_Soil_Rad.append(net_radiation_second_estimate)
                    hourly_Air_Soil_temp_dif.append(temperature_difference_soil)
                    hourly_rbhs.append(boundary_layer_resistance_soil)
                    hourly_rts.append(turbulence_resistance_soil)
    
        
        
                    
                self.hourly_Soil_Evap=hourly_Soil_Evap
                self.hourly_Soil_Rad=hourly_Soil_Rad
                self.hourly_Air_Soil_temp_dif=hourly_Air_Soil_temp_dif
                self.hourly_rbhs=hourly_rbhs
                self.hourly_rts=hourly_rts
                
                # Aggregate hourly data back to daily totals
                potential_evap_daily = Leaf.Leaf.aggregate_to_daily(hourly_Soil_Evap,  Day_Length)
                potential_SoilRad_daily = Leaf.Leaf.aggregate_to_daily(hourly_Soil_Rad,  Day_Length)
                self.potential_evap_daily=potential_evap_daily    
                self.potential_SoilRad_daily=potential_SoilRad_daily    
                self.Day_Air_Soil_temp_dif=(hourly_Air_Soil_temp_dif * wgauss).sum()
        
    
        def Update_Evaporation_if_WaterStress(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, daily_water_supply, soil_depth_1, Root_Depth, hourly_transpiration_Sunlit, hourly_transpiration_Shaded, hourly_Soil_Evap):
            gaussian_points = 5
            gaussian_weights_x = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
            gaussian_weights_w = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
            Latent_Heat_of_Vaporization = 2.4E6  # J/kg
            Volumetric_Heat_Capacity_Air = 1200  # J/m^3/°C
    
            hourly_Actual_Soil_Evap = []
            hourly_Air_Soil_Temperature_Difference = []
    
            for i, evap_potential, transp_sunlit, transp_shaded, solar_radiation, layer_resistance, turbulence_resistance in zip(range(gaussian_points), hourly_Soil_Evap, hourly_transpiration_Sunlit, hourly_transpiration_Shaded, self.hourly_Soil_Rad, self.hourly_Boundary_Layer_Resistance, self.hourly_Turbulence_Resistance):
                hour = 12 - 0.5 * Day_Length + Day_Length * gaussian_weights_x[i]
                Sin_Beam = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (hour - 12.) / 24.))
    
                # Diurnal availability of soil water supply
                water_supply_hourly = daily_water_supply * (Sin_Beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
    
                # Available water for evaporation from the top soil layer
                water_supply_for_evaporation = water_supply_hourly * (soil_depth_1 / Root_Depth)
    
                # Total canopy transpiration (sunlit + shaded)
                total_canopy_transpiration = transp_sunlit + transp_shaded
    
                # Transpiration from the top soil layer
                top_soil_layer_transpiration = total_canopy_transpiration * (soil_depth_1 / Root_Depth)
    
                # Maximum possible evaporation
                Max_possible_evap = evap_potential / (evap_potential + top_soil_layer_transpiration) * water_supply_for_evaporation
                Actual_Soil_Evap = Max_possible_evap if Max_possible_evap < evap_potential else evap_potential
    
                Actual_Air_Soil_Temperature_Difference = self.Limit_Function(-25., 25., (solar_radiation - Latent_Heat_of_Vaporization * Actual_Soil_Evap) * (layer_resistance + turbulence_resistance) / Volumetric_Heat_Capacity_Air)
                hourly_Air_Soil_Temperature_Difference.append(Actual_Air_Soil_Temperature_Difference)
                hourly_Actual_Soil_Evap.append(Actual_Soil_Evap)
    
            self.hourly_Actual_Soil_Evap = hourly_Actual_Soil_Evap
            # Aggregation to daily air-soil temperature difference and actual daily evaporation
            self.Day_Air_Soil_Temperature_Difference = np.dot(hourly_Air_Soil_Temperature_Difference, gaussian_weights_w).sum()
            # Assuming `aggregate_to_daily` is a method that sums or averages hourly data to a daily total
            self.Actual_Daily_Evaporation = np.dot(hourly_Actual_Soil_Evap, gaussian_weights_w).sum()
            print(self.Actual_Daily_Evaporation)
    
    
    
    
    
    
        def Calculate_Soil_Temperature(self, tmin, tmax):
            # Calculate the daily average temperature with a weighted average favoring tmax
            daily_avg_temperature = 0.29 * tmin + 0.71 * tmax
            
            # Calculate the nightly average temperature with a weighted average favoring tmin
            nightly_avg_temperature = 0.71 * tmin + 0.29 * tmax
            
            # Calculate the soil's average steady state temperature
            soil_steady_state_avg_temp = (daily_avg_temperature + self.Day_Air_Soil_temp_dif + nightly_avg_temperature) / 2.
            
            # Calculate the rate of change in soil temperature
            rate_of_change_soil_temp = (soil_steady_state_avg_temp - self.initial_soil_temp) / self.temperature_change_constant
            
            # Print the rate of change in soil temperature for debugging or information purposes
            print(rate_of_change_soil_temp)
            
            # Update class attributes with the new calculated values
            self.soil_steady_state_avg_temp = soil_steady_state_avg_temp
            self.rate_of_change_soil_temp = rate_of_change_soil_temp
    
    
    
        def Calculate_Water_Balance(self, days_after_planting, rainfall, total_canopy_transpiration):
                # Irrigation schedule: Days and corresponding irrigation amount
                irrigation_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0), (85.0, 0), (95.0, 0)]
                irrigation = sum([amount for day, amount in irrigation_schedule if day <= days_after_planting])
                
                # Sum of effective rainfall and irrigation
                rain_and_irrigation = rainfall + irrigation
                
                # Recharges to upper and lower soil layers
                recharge_to_upper_layer = min(10.0 * (self.Saturated_Water_Content - self.soil_water_content_upper_layer) * self.Root_Depth / self.temperature_change_constant, rain_and_irrigation)
                recharge_to_lower_layer = min(10.0 * (self.Saturated_Water_Content - self.soil_water_content_lower_layer) * (150.0 - self.Root_Depth) / self.temperature_change_constant, rain_and_irrigation - recharge_to_upper_layer)
                
                # Groundwater recharge and water loss
                groundwater_recharge = max(0., rain_and_irrigation - recharge_to_upper_layer - recharge_to_lower_layer)
                water_loss_lower_layer = recharge_to_lower_layer - 10.0 * (self.soil_water_content_lower_layer - self.Residual_Water_Content) * self.Root_Depth
                water_loss_upper_layer = recharge_to_upper_layer + 10.0 * (self.soil_water_content_lower_layer - self.Residual_Water_Content) * self.Root_Depth - self.switch_function(self.water_supply_switch, 0.0, total_canopy_transpiration) + 0.1
                
                # Update instance attributes with calculated values
                self.rain_and_irrigation = rain_and_irrigation
                self.recharge_to_upper_layer = recharge_to_upper_layer
                self.recharge_to_lower_layer = recharge_to_lower_layer
                self.groundwater_recharge = groundwater_recharge
                self.water_loss_lower_layer = water_loss_lower_layer
                self.water_loss_upper_layer = water_loss_upper_layer
    
       
        
        def Calculate_Soil_N_Dynamics(self):
            N_mineralization_rate_constant = 24 * (math.exp(17.753 - 6350.5 / (self.soil_temperature + 273))) / 168
            N_mineralization_factor = 1 - math.exp(-N_mineralization_rate_constant)
            
            # Adjust for soil water content
            if fraction_transpirable_soil_water < 0.9:
                relative_N = 1.111 * fraction_transpirable_soil_water
            else:
                relative_N = 10 - 10 * fraction_transpirable_soil_water
            if relative_N < 0:
                relative_N = 0
            
            # Calculate net N mineralization
            net_N_mineralization = mineralizable_organic_N * relative_N * N_mineralization_factor
            # Apply threshold of 200 mgN.L-1
            net_N_mineralization = net_N_mineralization * (0.0002 - self.N_concentration) / 0.0002
            if net_N_mineralization < 0:
                net_N_mineralization = 0
            
           
            
           
            
            if self.Days_after_planting == self.Fertilizer_applications_DAP[0]:  # Check if today matches any fertilization day
                
                N_fertilization = self.Fertilizer_applications_amount[0]  # gN.m-2, N to be applied
                volatilization_factor = self.Fraction_volatilization[0] / 100  # Convert percentage to decimal
                
                self.Fertilizer_applications_DAP=Fertilizer_applications_DAP[1:]
                self.Fertilizer_applications_amount=Fertilizer_applications_amount[1:]
                self.volatilization_factor=volatilization_factor[1:]
           
                
                N_volatilization = volatilization_factor * N_fertilization  # gN.m-2, N volatilized
           
            else:
               
               N_fertilization = N_volatilization = 0
                 
             
             
         
            # Calculate N leaching based on soluble N, water, and drainage
            #Assuming no water drainage otherwise         N_leaching = soluble_N * (DRAIN1 / (water + DRAIN1))  # gN.m-2
            N_leaching = self.soluble_N 
            
            # Apply a threshold of 1 mgN.L-1 to the N concentration
            if self.N_concentration <= 0.000001:
                N_leaching = 0
            
        
         
                 
        
            if fraction_transpirable_soil_water > 1:  # Condition for water saturation
                adjusted_N_concentration = self.N_concentration
                
            # Apply threshold of 400 mgN.L-1 to the N concentration
            if adjusted_N_concentration > 0.0004:
                adjusted_N_concentration = 0.0004
            
            # Calculate denitrification rate constant based on temperature
            denitrification_rate_constant = 6 * math.exp(0.07735 * self.soil_temperature - 6.593)
            
            # Calculate N denitrification per g of water
            N_denitrification = adjusted_N_concentration * (1 - math.exp(-denitrification_rate_constant))
            
            # Convert N denitrification to gN.m-2 based on water volume
            N_denitrification = N_denitrification * water * 1000
            
            
                    
            
            
            # Calculate the fraction of the top layer with roots
            fraction_top_layer_with_roots = self.Root_Depth / self.Soil_Depth_1
            if fraction_top_layer_with_roots > 1:
                fraction_top_layer_with_roots = 1
            
            # Update soluble N content
            self.soluble_N = (self.soluble_N + net_N_mineralization + 
                                N_fertilization - N_volatilization -
                                N_leaching - N_denitrification - self.total_nitrogen_uptake)
            
            # Recalculate N concentration in gN.g-1 H2O
            self.N_concentration = self.soluble_N / (water * 1000)
            
            # Calculate available soluble N in the root zone based on a threshold of 1 mgN.L-1
            available_soluble_N = (self.N_concentration - 0.000001) * ATSW1 * 1000 * fraction_top_layer_with_roots
            if available_soluble_N < 0:
                available_soluble_N = 0
            
    
    
    
    
    










    # def Organic_Carbon_Composition(self):


    #     temperature_factor = 47.9 / (1. + np.exp(106. / (self.initial_soil_temp + 18.3)))
    #     moisture_factor = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * (self.soil_water_content_upper_layer + self.soil_water_content_lower_layer) / 10. / 150. / (self.field_capacity_water_content - self.Residual_Water_Content))
        
    #     # Extract nitrogen levels
    #     nitrate_upper_layer, ammonium_upper_layer, nitrate_lower_layer, ammonium_lower_layer = self.nnuli, self.nauli, self.nnlli, self.nalli
        
    #     # Calculate decomposition rates for plant material
    #     decomposable_pm_rate_change = self.switch_function(nitrate_upper_layer + ammonium_upper_layer + nitrate_lower_layer + ammonium_lower_layer - self.residual_ammonium - self.residual_nitrate, 0., self.dpm_decomposition_rate)
    #     resistant_pm_rate_change = self.switch_function(nitrate_upper_layer + ammonium_upper_layer + nitrate_lower_layer + ammonium_lower_layer - self.residual_ammonium - self.residual_nitrate, 0., self.rpm_decomposition_rate)
        
    #     clay_bonus = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * self.clay_percentage))

    #     carbon_nitrogen_dpm_ratio = (self.initial_decomposable_plant_material + self.resistant_plant_material_initial) / self.Avoid_Zero_Division(self.decomposable_plant_nitrogen_initial + self.resistant_plant_nitrogen_initial)
    #     decomposed_biomass = self.biomass_incorporation_rate * (1 - math.exp(-temperature_factor * moisture_factor * self.biomass_incorporation_rate / 365)) / self.temperature_change_constant
    #     decomposed_humus = self.humus_initial_content * (1 - math.exp(-temperature_factor * moisture_factor * self.humification_rate / 365)) / self.temperature_change_constant

    #     dpm_rate_adjusted = self.switch_function(1.0 / self.Avoid_Zero_Division(carbon_nitrogen_dpm_ratio) - 1.0 / (8.5 * (1.0 + clay_bonus)), decomposable_pm_rate_change, self.dpm_decomposition_rate)
    #     rpm_rate_adjusted = self.switch_function(1.0 / self.Avoid_Zero_Division(carbon_nitrogen_dpm_ratio) - 1.0 / (8.5 * (1.0 + clay_bonus)), resistant_pm_rate_change, self.rpm_decomposition_rate)
        
    #     decomposed_dpm = self.initial_decomposable_plant_material * (1.0 - math.exp(-temperature_factor * moisture_factor * dpm_rate_adjusted / 365)) / self.temperature_change_constant
    #     decomposed_rpm = self.resistant_plant_material_initial * (1.0 - math.exp(-temperature_factor * moisture_factor * rpm_rate_adjusted / 365)) / self.temperature_change_constant
    #     decomposed_dpn = self.decomposable_plant_nitrogen_initial * (1.0 - math.exp(-temperature_factor * moisture_factor * dpm_rate_adjusted / 365)) / self.temperature_change_constant
    #     decomposed_rpn = self.resistant_plant_nitrogen_initial * (1.0 - math.exp(-temperature_factor * moisture_factor * rpm_rate_adjusted / 365)) / self.temperature_change_constant

    #     mineralized_nitrogen = 1.0 / 8.5 * (decomposed_biomass + decomposed_humus) + decomposed_dpn + decomposed_rpn - 1.0 / 8.5 / (1.0 + clay_bonus) * \
    #                           (decomposed_dpm + decomposed_rpm + decomposed_biomass + decomposed_humus)
    #     mineralized_nitrogen_upper_layer = (1. - math.exp(-0.065 * self.Root_Depth)) * mineralized_nitrogen
    #     mineralized_nitrogen_lower_layer = math.exp(-0.065 * self.Root_Depth) * mineralized_nitrogen
        
    #     # Update class attributes with the calculated values
    #     self.temperature_factor = temperature_factor
    #     self.clay_bonus = clay_bonus
    #     self.decomposed_dpm = decomposed_dpm
    #     self.decomposed_rpm = decomposed_rpm
    #     self.decomposed_dpn = decomposed_dpn
    #     self.decomposed_rpn = decomposed_rpn
    #     self.decomposed_biomass = decomposed_biomass
    #     self.decomposed_humus = decomposed_humus
    #     self.mineralized_nitrogen = mineralized_nitrogen
    #     self.mineralized_nitrogen_upper_layer = mineralized_nitrogen_upper_layer
    #     self.mineralized_nitrogen_lower_layer = mineralized_nitrogen_lower_layer

    
    

    def Calculate_Nitrogen_Uptake(self, nitrogen_demand, nitrogen_fixation_rate):
        # Calculate supply of ammonium nitrogen from the soil
        ammonium_nitrogen_supply_as = max(0., self.ammonium_upper_layer + (self.mineralized_ammonium_upper_layer - self.nitrate_upper_layer) * self.temperature_change_constant - self.Root_Depth / 150. * self.residual_ammonium) / self.temperature_change_constant
        
        # Calculate supply of nitrate nitrogen from the soil, factoring in water stress
        nitrate_nitrogen_supply_ns = max(0., self.nitrate_upper_layer + (self.mineralized_nitrate_upper_layer - self.denitrification_upper_layer) * self.temperature_change_constant - self.Root_Depth / 150. * self.residual_nitrate) / self.temperature_change_constant * self.nitrogen_stress_water_index
        
        # Determine actual nitrogen supply based on soil water index
        ammonium_nitrogen_supply = self.switch_function(self.water_supply_switch, self.ammonium_nitrogen_input_rate, ammonium_nitrogen_supply_as)
        nitrate_nitrogen_supply = self.switch_function(self.water_supply_switch, self.nitrate_nitrogen_input_rate, nitrate_nitrogen_supply_ns)
        
        # Total nitrogen supply
        total_nitrogen_supply = ammonium_nitrogen_supply + nitrate_nitrogen_supply
        
        # Calculate nitrogen uptake, respecting the demand and fixation rate
        ammonium_nitrogen_uptake = min(ammonium_nitrogen_supply, ammonium_nitrogen_supply / max(1e-10, total_nitrogen_supply) * max(0, nitrogen_demand - nitrogen_fixation_rate / self.temperature_change_constant))
        nitrate_nitrogen_uptake = min(nitrate_nitrogen_supply, nitrate_nitrogen_supply / max(1e-10, total_nitrogen_supply) * max(0, nitrogen_demand - nitrogen_fixation_rate / self.temperature_change_constant))
        
        # Total nitrogen uptake, ensuring it does not exceed demand
        total_nitrogen_uptake = max(0, ammonium_nitrogen_uptake + nitrate_nitrogen_uptake + min(nitrogen_demand, nitrogen_fixation_rate / self.temperature_change_constant))
        
        # Debugging print statements can be commented out or removed in production
        print(ammonium_nitrogen_uptake, nitrate_nitrogen_uptake, nitrogen_demand, nitrogen_fixation_rate, self.temperature_change_constant)
        
        # Update class attributes with the calculated values
        self.total_nitrogen_uptake = total_nitrogen_uptake
        self.nitrate_nitrogen_uptake = nitrate_nitrogen_uptake
        self.ammonium_nitrogen_uptake = ammonium_nitrogen_uptake
        self.total_nitrate_nitrogen_supply = nitrate_nitrogen_supply


    
    
    
    
    # def Organic_Nitrogen_Composition(self, rainfall):

    #     # Calculate mineralized ammonium in upper and lower layers
    #     mineralized_ammonium_upper_layer = self.switch_function(self.mineralized_nitrogen, -min((self.ammonium_lower_layer - self.Root_Depth / 150. * self.residual_ammonium) / self.temperature_change_constant, -self.mineralized_nitrogen_upper_layer), self.mineralized_nitrogen_upper_layer)
    #     mineralized_ammonium_lower_layer = self.switch_function(self.mineralized_nitrogen, -min((self.ammonium_lower_layer - (150. - self.Root_Depth) / 150. * self.residual_ammonium) / self.temperature_change_constant, -self.mineralized_nitrogen_lower_layer), self.mineralized_nitrogen_lower_layer)

    #     # Calculate mineralized nitrate in upper and lower layers
    #     mineralized_nitrate_upper_layer = self.switch_function(self.mineralized_nitrogen, -min(self.nitrate_upper_layer/self.temperature_change_constant, -self.mineralized_nitrogen_upper_layer + mineralized_ammonium_upper_layer), 0.)
    #     mineralized_nitrate_lower_layer = self.switch_function(self.mineralized_nitrogen, -min(self.nitrate_lower_layer/self.temperature_change_constant, -self.mineralized_nitrogen_lower_layer + mineralized_ammonium_lower_layer), 0.)

    #     # Calculate moisture factors for upper and lower layers
    #     moisture_factor_upper_layer = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * self.soil_water_content_upper_layer / 10. / self.Root_Depth / (self.field_capacity_water_content - self.Residual_Water_Content))
    #     moisture_factor_lower_layer = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * self.soil_water_content_lower_layer / 10. / (150. - self.Root_Depth) / (self.field_capacity_water_content - self.Residual_Water_Content))

    #     # Calculate nitrification in upper and lower layers
    #     nitrification_upper_layer = max(0., (self.ammonium_upper_layer + mineralized_ammonium_upper_layer * self.temperature_change_constant - self.Root_Depth / 150 * self.residual_ammonium)) * (1 - np.exp(-self.temperature_factor * moisture_factor_upper_layer * 0.6 / 7)) / self.temperature_change_constant
    #     nitrification_lower_layer = max(0., (self.ammonium_lower_layer + mineralized_ammonium_lower_layer * self.temperature_change_constant - (150 - self.Root_Depth) / 150 * self.residual_ammonium)) * (1 - np.exp(-self.temperature_factor * moisture_factor_lower_layer * 0.6 / 7)) / self.temperature_change_constant

    #     # Calculate CO2 respiration
    #     respiration_co2 = self.clay_bonus / (1.0 + self.clay_bonus) * (self.decomposed_dpm + self.decomposed_rpm + self.decomposed_biomass + self.decomposed_humus)

    #     # Calculate denitrification in upper and lower layers
    #     denitrification_upper_layer = .0005 * max(0., self.nitrate_upper_layer + mineralized_nitrate_upper_layer * self.temperature_change_constant - self.Root_Depth / 150. * self.residual_nitrate) * respiration_co2 * (1. - np.exp(-0.065 * self.Root_Depth))
    #     denitrification_lower_layer = .0005 * max(0., self.nitrate_lower_layer + mineralized_nitrate_lower_layer * self.temperature_change_constant - (150. - self.Root_Depth) / 150. * self.residual_nitrate) * respiration_co2 * np.exp(-0.065 * self.Root_Depth)

    #     # Calculate water stress factor
    #     water_stress_factor = min(1, self.soil_water_content_upper_layer / (self.Root_Depth * 10 * (self.field_capacity_water_content - self.Residual_Water_Content)))

    #     # Calculate total nitrogen and ammonia volatilization
    #     total_ammonium = self.ammonium_upper_layer + self.ammonium_lower_layer
    #     total_nitrate = self.nitrate_upper_layer + self.nitrate_lower_layer
    #     total_mineral_nitrogen = total_ammonium + total_nitrate
    #     ammonia_volatilization = self.switch_function(rainfall - 1., 0.15, 0.) * self.soil_resistance_to_evaporation

    #     # Update class attributes or return calculated values as needed
    #     self.nitrification_upper_layer = nitrification_upper_layer
    #     self.nitrification_lower_layer = nitrification_lower_layer
    #     self.mineralized_ammonium_upper_layer = mineralized_ammonium_upper_layer
    #     self.mineralized_ammonium_lower_layer = mineralized_ammonium_lower_layer
    #     self.mineralized_nitrate_upper_layer = mineralized_nitrate_upper_layer
    #     self.mineralized_nitrate_lower_layer = mineralized_nitrate_lower_layer
    #     self.denitrification_upper_layer = denitrification_upper_layer
    #     self.denitrification_lower_layer = denitrification_lower_layer
    #     self.water_stress_factor = water_stress_factor
    #     self.ammonia_volatilization = ammonia_volatilization




    # def Calculate_Soil_Organic_Nitrogen_ChangeRate(self, days_after_planting, root_depth_ratio):
    #     ammonium_fertilization_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]
    #     nitrate_fertilization_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]

    #     # Calculate ammonium and nitrate fertilization amounts based on schedule and days after planting
    #     ammonium_fertilizer_application = sum([self.Function_Adjust_Switch(time - days_after_planting, 0., amount, 0.) for time, amount in ammonium_fertilization_schedule])
    #     nitrate_fertilizer_application = sum([self.Function_Adjust_Switch(time - days_after_planting, 0., amount, 0.) for time, amount in nitrate_fertilization_schedule])

    #     # Residual soil fertilizer after application (assuming a certain percentage is immediately available)
    #     residual_soil_ammonium = ammonium_fertilizer_application - self.soil_ammonium_volatilization_rate / 3
        
    #     # Leaching losses for ammonium and nitrate in lower soil layer
    #     leaching_ammonium_lower_layer = max(0, self.nitrate_lower_layer + (self.mineralized_nitrate_lower_layer - self.denitrification_upper_layer) * self.temperature_change_constant - (150 - self.Root_Depth) / 150 * self.groundwater_recharge) * min(self.groundwater_recharge / self.Saturated_Water_Content / (150 - self.Root_Depth) / 10, 1)
    #     leaching_nitrate_upper_layer = max(0, (self.total_nitrate_nitrogen_supply - self.nitrate_nitrogen_uptake) * self.temperature_change_constant - self.Root_Depth / 150 * self.groundwater_recharge) * min((self.rain_and_irrigation - self.recharge_to_upper_layer) / self.Saturated_Water_Content / self.Root_Depth / 10, 1)
        
    #     # Nitrogen layer adjustments due to root depth ratio
    #     layer_adjustment_ammonium = root_depth_ratio / (150.0 - self.Root_Depth) * self.ammonium_lower_layer
    #     layer_adjustment_nitrate = root_depth_ratio / (150.0 - self.Root_Depth) * self.nitrate_lower_layer

    #     # Update class attributes with the calculated values
    #     self.residual_soil_ammonium = residual_soil_ammonium
    #     self.leaching_ammonium_lower_layer = leaching_ammonium_lower_layer
    #     self.leaching_nitrate_upper_layer = leaching_nitrate_upper_layer
    #     self.layer_adjustment_ammonium = layer_adjustment_ammonium
    #     self.layer_adjustment_nitrate = layer_adjustment_nitrate
    #     self.ammonium_fertilizer_application = ammonium_fertilizer_application
    #     self.nitrate_fertilizer_application = nitrate_fertilizer_application


    # def Calculate_Soil_Organic_Carbon_ChangeRate(self, litter_carbon, litter_nitrogen):
    #     # Decomposition and transformation rates for soil organic carbon
    #     decomposable_pm_change_rate = litter_carbon * self.plant_material_dpm_rpm_ratio / (1. + self.plant_material_dpm_rpm_ratio) - self.decomposed_dpm
    #     resistant_pm_change_rate = litter_carbon * 1. / (1. + self.plant_material_dpm_rpm_ratio) - self.decomposed_rpm
    #     decomposable_pn_change_rate = litter_nitrogen / (1. + 40. / self.plant_material_dpm_rpm_ratio / 100.) - self.decomposed_dpn
    #     resistant_pn_change_rate = litter_nitrogen / (1. + 100. * self.plant_material_dpm_rpm_ratio / 40.) - self.decomposed_rpn
        
    #     # Biochar and humus change rates based on decomposition
    #     biochar_change_rate = 0.46 / (1.0 + self.clay_bonus) * (self.decomposed_dpm + self.decomposed_rpm + self.decomposed_biomass + self.decomposed_humus) - self.decomposed_biomass
    #     humus_change_rate = 0.54 / (1.0 + self.clay_bonus) * (self.decomposed_dpm + self.decomposed_rpm + self.decomposed_biomass + self.decomposed_humus) - self.decomposed_humus
        
    #     # Ammonium and nitrate nitrogen in the upper and lower layers
    #     ammonium_upper_layer_change = self.ammonium_fertilizer_application + self.mineralized_ammonium_upper_layer + self.layer_adjustment_ammonium - self.switch_function(self.water_supply_switch, 0.0, self.ammonium_nitrogen_uptake) - self.nitrification_upper_layer - self.ammonia_volatilization
    #     ammonium_lower_layer_change = self.mineralized_ammonium_lower_layer - self.layer_adjustment_ammonium - self.nitrification_lower_layer
    #     nitrate_lower_layer_change = self.leaching_nitrate_upper_layer + self.mineralized_nitrate_lower_layer + self.nitrification_lower_layer - self.layer_adjustment_nitrate - self.denitrification_lower_layer - self.leaching_ammonium_lower_layer
    #     nitrate_upper_layer_change = self.nitrate_fertilizer_application + self.mineralized_nitrate_upper_layer + self.nitrification_upper_layer + self.layer_adjustment_nitrate - self.switch_function(self.water_supply_switch, 0.0, self.nitrate_nitrogen_uptake) - self.denitrification_upper_layer - self.leaching_nitrate_upper_layer
        
    #     # Update class attributes with the calculated values
    #     self.decomposable_pm_change_rate = decomposable_pm_change_rate
    #     self.resistant_pm_change_rate = resistant_pm_change_rate
    #     self.decomposable_pn_change_rate = decomposable_pn_change_rate
    #     self.resistant_pn_change_rate = resistant_pn_change_rate
    #     self.biochar_change_rate = biochar_change_rate
    #     self.humus_change_rate = humus_change_rate
    #     self.ammonium_upper_layer_change = ammonium_upper_layer_change
    #     self.ammonium_lower_layer_change = ammonium_lower_layer_change
    #     self.nitrate_lower_layer_change = nitrate_lower_layer_change
    #     self.nitrate_upper_layer_change = nitrate_upper_layer_change

        
        
        
        
    def switch_fun(self, x, y1, y2):
        return y1 if x < 0 else y2
    def Avoid_Zero_Division(self, x):
        
        if x != 0:
            out = x
        else:
            out = 1.0
        return out
    def Limit_Function(self,xl, xh, x):
       
        if xl <= x <= xh:
            out = x
        elif x > xh:
            out = xh
        else:
            out = xl
        return out

    
    def Update_State_Variables(self):
        # Update total nitrogen uptake
        self.total_nitrogen_uptake += self.nitrogen_uptake
        
        # Update soil temperature with rate of temperature change
        self.initial_soil_temperature += self.rate_of_change_soil_temp
    
        # Update water uptake and loss limits with respective rate changes
        self.water_uptake_limit += self.rate_of_change_water_uptake_limit
        self.water_loss_limit += self.rate_of_change_water_loss_limit
        
        # Update organic matter and nitrogen pools with their respective rate changes
        self.decomposable_plant_material += self.rate_of_change_decomposable_pm
        self.resistant_plant_material += self.rate_of_change_resistant_pm
        self.biomass_incorporated += self.rate_of_change_biomass
        self.humus_content += self.rate_of_change_humus
        self.decomposable_plant_nitrogen += self.rate_of_change_decomposable_pn
        self.resistant_plant_nitrogen += self.rate_of_change_resistant_pn
        
        # Update total nitrogen leached and nitrogen content in upper and lower soil layers
        self.total_nitrogen_leached += self.leached_ammonium_nitrogen_total
        self.ammonium_upper_layer += self.rate_of_change_ammonium_upper_layer
        self.ammonium_lower_layer += self.rate_of_change_ammonium_lower_layer
        self.nitrate_upper_layer += self.rate_of_change_nitrate_upper_layer
        self.nitrate_lower_layer += self.rate_of_change_nitrate_lower_layer
        
        # Update soil ammonium volatilization rate
        self.soil_ammonium_volatilization_rate += self.rate_of_change_soil_ammonium_volatilization


     