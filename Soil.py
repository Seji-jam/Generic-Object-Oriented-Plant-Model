import numpy as np
import math
import Leaf
               
                
Latent_Heat_of_Water_Vaporization = 2.4E6  # Latent heat of vaporization (J/kg)
Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity of air (J/m^3/°C)
wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
class Soil:
    def __init__(self, Residual_Soil_Moisture, Saturated_Soil_Moisture, Field_Capacity, Initial_Soil_Moisture,
                 Soil_Depth, Top_Layer_Depth, Soil_Evaporative_Depth,clay_percentage,sand_percentage,Drainage_Factor,
                  initial_soil_temp,  soil_resistance_to_evaporation,
                 fraction_soil_greater_than_2mm,soil_bulk_density,organic_N_percentage,fraction_N_for_mineralization,
                 nitrate_concentration_ppm,ammonium_concentration_ppm,
                 Fertilizer_applications_count,Fertilizer_applications_amount,Fertilizer_applications_DAP,Fraction_volatilization,
                 Soil_Dynamic_Temperature_Factor):
        
        self.Residual_Soil_Moisture = Residual_Soil_Moisture  # Residual water content
        self.Saturated_Soil_Moisture = Saturated_Soil_Moisture  # Saturated water content
        # self.lodging_condition = lodging_condition  # Lodging condition
        self.Soil_Evaporative_Depth = Soil_Evaporative_Depth  # Soil bulk density
        self.Soil_Depth=Soil_Depth
        self.Top_Layer_Depth=Top_Layer_Depth
        self.initial_soil_temp = initial_soil_temp  # Initial soil temperature
        # self.initial_decomposable_plant_material = initial_decomposable_plant_material  # Initial decomposable plant material in soil
        self.Field_Capacity = Field_Capacity  # Water content at field capacity
        self.Initial_Soil_Moisture=Initial_Soil_Moisture
        self.Drainage_Factor=Drainage_Factor
        self.clay_percentage = clay_percentage  # Clay content in soil
        self.soil_resistance_to_evaporation = soil_resistance_to_evaporation  # Soil resistance to evaporation
        self.sand_percentage = sand_percentage  # Percent of Sand in soil
        self.fraction_soil_greater_than_2mm = fraction_soil_greater_than_2mm/100  # Fraction soil > 2 mm (soil coarse fraction)
        self.soil_bulk_density = soil_bulk_density  # Soil bulk density (g.cm-3)
        self.organic_N_percentage = organic_N_percentage  # Organic N (%)
        self.fraction_N_for_mineralization = fraction_N_for_mineralization  # Fraction Organic N available for mineralization
        self.nitrate_concentration_ppm = nitrate_concentration_ppm  # Nitrate concentration (ppm = mg.kg−1)
        self.ammonium_concentration_ppm = ammonium_concentration_ppm  # Ammonium concentration (ppm = mg.kg−1)
        self.Soil_Evaporative_Depth = Soil_Evaporative_Depth  # Depth, assuming necessary for calculation
        # self.water = water  # Water, assuming necessary for calculation
        self.Fertilizer_applications_count = Fertilizer_applications_count  
        self.Fertilizer_applications_amount = Fertilizer_applications_amount  
        self.Soil_Dynamic_Temperature_Factor =Soil_Dynamic_Temperature_Factor
        self.Fertilizer_applications_DAP=[Fertilizer_applications_DAP]
        self.Fertilizer_applications_amount=[Fertilizer_applications_amount]
        self.Fraction_volatilization=[Fraction_volatilization]
       
        self.Current_Soil_Moisture_Top_Layer=Initial_Soil_Moisture
        self.Current_Soil_Moisture_Whole_Depth=Initial_Soil_Moisture

        # Available_Soil_Water_Initial= min(Initial_Soil_Moisture-Residual_Soil_Moisture, -Field_Capacity-Residual_Soil_Moisture)
        # self.Available_Soil_Water=Available_Soil_Water_Initial
        self.time_change_constant = 1 ########################## this is set as 1 day ... if the model changes to Hourly this constant needs to be changed to 1/24

        

        self.water_stress_factor = 0  # Factor indicating water stress
        self.Nitrogen_uptake = 0  # Nitrogen uptake rate


        self.Hourly_Soil_Evap=[]
        self.Hourly_Soil_Rad=[]



        
        # Initial conditions and parameter calculations
        root_depth_ini = max(2.0, Soil_Evaporative_Depth)
        # Root_zone_soil_moisture_initial = 10.0 * (self.Initial_Soil_Moisture - Residual_Soil_Moisture) * root_depth_ini
        # Below_Root_zone_soil_moisture_initial = 10.0 * (self.Initial_Soil_Moisture - Residual_Soil_Moisture) * (150.0 - root_depth_ini)


    
        self.average_soil_temperature = initial_soil_temp
        self.rate_of_change_soil_temp = 0
    
        self.soil_water_content_upper_limit = 0
        self.soil_water_content_lower_limit = 0
        self.nitrate_nitrogen_supply = 0
        self.ammonium_fertilizer_application = 0
        self.nitrate_fertilizer_application = 0

        self.ammonia_volatilization = 0
        # State variables
        self.Cumulative_Nitrogen_uptake = 0 
        self.Root_Depth = root_depth_ini  
        self.total_nitrogen_leached = 0
        self.soil_ammonium_volatilization_rate = 0

        self.potential_evap_daily = 0
        self.potential_SoilRad_daily = 0
        self.Day_Air_Soil_temp_dif  = 0
        self.Actual_evap_daily =0
        self.Hourly_boundary_layer_resistance_soil=[]
        self.Hourly_turbulence_resistance_soil=[]
        self.available_soluble_N=0
        
        self.Current_Soil_Water_Content_Top_Layer = 0
        self.Total_Soil_Water_Content_Top_Layer = 0
        self.Theoritical_Soil_Water_Content_Top_Layer = 0
        self.Soil_Moisture_Fraction_Whole_Depth=0
        self.Total_Soil_Water_Content_Whole_Depth=0
        self.Current_Soil_Water_Content_Whole_Depth=0

    def switch_function(self,x, y1, y2):
        return y1 if x < 0 else y2
        
    
    def Calculate_Soil_Water_Content(self):

        Current_Soil_Water_Content_Top_Layer=(self.Current_Soil_Moisture_Top_Layer-self.Residual_Soil_Moisture)*self.Top_Layer_Depth
        Total_Soil_Water_Content_Top_Layer = (self.Residual_Soil_Moisture * self.Top_Layer_Depth) + Current_Soil_Water_Content_Top_Layer
        
        Theoritical_Soil_Water_Content_Top_Layer=(self.Field_Capacity-self.Residual_Soil_Moisture)*self.Top_Layer_Depth

        Current_Soil_Water_Content_Whole_Depth=(self.Current_Soil_Moisture_Whole_Depth-self.Residual_Soil_Moisture)*self.Soil_Depth
        Total_Soil_Water_Content_Whole_Depth = (self.Residual_Soil_Moisture * self.Soil_Depth) + Current_Soil_Water_Content_Whole_Depth
        Soil_Moisture_Fraction_Whole_Depth=Current_Soil_Water_Content_Whole_Depth/Total_Soil_Water_Content_Whole_Depth
        
        print(self.Current_Soil_Moisture_Top_Layer,self.Residual_Soil_Moisture,Current_Soil_Water_Content_Top_Layer)
        # Update the class attributes with the calculated values
        self.Current_Soil_Water_Content_Top_Layer = Current_Soil_Water_Content_Top_Layer
        self.Total_Soil_Water_Content_Top_Layer = Total_Soil_Water_Content_Top_Layer
        self.Theoritical_Soil_Water_Content_Top_Layer=Theoritical_Soil_Water_Content_Top_Layer
        
        self.Soil_Moisture_Fraction_Whole_Depth=Soil_Moisture_Fraction_Whole_Depth
        self.Total_Soil_Water_Content_Whole_Depth=Total_Soil_Water_Content_Whole_Depth
        self.Current_Soil_Water_Content_Whole_Depth=Current_Soil_Water_Content_Whole_Depth



    
    def Calculate_Soil_Potential_Evaporation(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination,
                                             Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature,
                                             Min_Temperature, Vapour_Pressure_Deficit, Wind_Speed, 
                                             soil_resistance_to_evaporation,Leaf_Blade_Angle, Total_Leaf_Area_Index,
                                             Light_Extinction_Coefficient, Hourly_transpiration_Shaded, Hourly_transpiration_Sunlit):
            
            # Convert daily climate data to Hourly data
            Hourly_climate_data, gaussian_weights = Leaf.Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature, Min_Temperature, Wind_Speed)
            
            Hourly_Soil_Evap = []
            Hourly_Air_Soil_temp_dif = []
            Hourly_Soil_Rad = []
            Hourly_boundary_layer_resistance_soil=[]
            Hourly_turbulence_resistance_soil=[]

            for (Hourly_Solar_Constant, Hourly_temperature, Hourly_sin_beam, Hourly_Solar_Radiation, Hourly_wind_speed), shaded_transpiration, sunlit_transpiration in zip(Hourly_climate_data, Hourly_transpiration_Shaded, Hourly_transpiration_Sunlit):
                
                Incoming_PAR = 0.5 * Hourly_Solar_Radiation  # Partition of incoming solar radiation to PAR
                Incoming_NIR = 0.5 * Hourly_Solar_Radiation  # Partition of incoming solar radiation to NIR
                
                # Atmospheric transmissivity based on incoming PAR and solar geometry
                atmospheric_transmissivity = Incoming_PAR / (0.5 * Hourly_Solar_Constant * Hourly_sin_beam)
    
                # Saturation vapor pressure and vapor pressure deficit calculations
                saturation_vapor_pressure = 0.611 * math.exp(17.4 * Hourly_temperature / (Hourly_temperature + 239.))
                vapor_pressure_deficit = max(0., saturation_vapor_pressure - Vapour_Pressure_Deficit)
                
                
                slope_vapor_pressure_temperature = 4158.6 * saturation_vapor_pressure / (Hourly_temperature + 239.) ** 2

                # Turbulence resistance for soil using logarithmic wind profile
                turbulence_resistance_soil = 0.74 * (np.log(56.)) ** 2 / (0.4 ** 2 * Hourly_wind_speed)
                
                # Boundary layer resistance for soil considering wind speed and leaf area
                boundary_layer_resistance_soil = 172. * np.sqrt(0.05 / max(0.1, Hourly_wind_speed * np.exp(-Light_Extinction_Coefficient * Total_Leaf_Area_Index)))
                
            
            
            
                water_boundary_layer_resistance_soil = 0.93 * boundary_layer_resistance_soil
                
                # Calculation of the fractional diffuse light (frdf) based on atmospheric transmissivity
                if atmospheric_transmissivity < 0.22:
                    fractional_diffuse_light = 1
                elif 0.22 < atmospheric_transmissivity <= 0.35:
                    fractional_diffuse_light = 1 - 6.4 * (atmospheric_transmissivity - 0.22) ** 2
                else:
                    fractional_diffuse_light = 1.47 - 1.66 * atmospheric_transmissivity
    
                # Ensuring a minimum threshold for the fractional diffuse light
                fractional_diffuse_light = max(fractional_diffuse_light, 0.15 + 0.85 * (1 - np.exp(-0.1 / Hourly_sin_beam)))
        
                # Incoming diffuse PAR (PARDF) and direct PAR (PARDR) based on the calculated fractional diffuse light
                PAR_Diffuse = Incoming_PAR * fractional_diffuse_light
                PAR_Direct = Incoming_PAR - PAR_Diffuse
    
                # Extinction and RefLection coefficients
                # Convert leaf blade angle from degrees to radians for calculations
                Leaf_Blade_Angle_Radians = Leaf_Blade_Angle * np.pi / 180
                
                # Calculate the direct beam extinction coefficient for PAR using leaf blade angle
                Direct_Beam_Extinction_Coefficient_PAR = Leaf.Leaf.KDR_Coeff(Hourly_sin_beam, Leaf_Blade_Angle_Radians)
                
                # Leaf scattering coefficients for PAR and NIR
                Scattering_Coefficient_PAR = 0.2  # Scattering coefficient for Photosynthetically Active Radiation
                Scattering_Coefficient_NIR = 0.8  # Scattering coefficient for Near-Infrared Radiation
                
                # Calculate diffuse extinction coefficients for PAR and NIR using total leaf area index and leaf blade angle
                Diffuse_Extinction_Coefficient_PAR = Leaf.Leaf.KDF_Coeff(Total_Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
                Diffuse_Extinction_Coefficient_NIR = Leaf.Leaf.KDF_Coeff(Total_Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
                
                # Calculate beam and canopy reflection coefficients for PAR and NIR
                Beam_Reflection_Coefficient_PAR, Canopy_Reflection_Coefficient_PAR = Leaf.Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient_PAR)
                Beam_Reflection_Coefficient_NIR, Canopy_Reflection_Coefficient_NIR = Leaf.Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient_PAR)
    
                
                # Incoming diffuse NIR (NIRDF) and direct NIR (NIRDR)
                NIR_Diffuse = Incoming_NIR * fractional_diffuse_light
                NIR_Direct = Incoming_NIR - NIR_Diffuse
                
                # Absorbed total radiation (PAR + NIR) by soil
                soil_par_reflection = 0.1  # Soil PAR reflection coefficient
                # Soil NIR reflection coefficient, varying with soil water content (wcul)
                soil_nir_reflection = self.switch_function(self.Current_Soil_Water_Content_Top_Layer - 0.5, 0.52 - 0.68 * self.Current_Soil_Water_Content_Top_Layer, 0.18)
                absorbed_total_radiation_by_soil = (1 - soil_par_reflection) * (PAR_Direct * np.exp(-Beam_Reflection_Coefficient_PAR * Total_Leaf_Area_Index) + PAR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_PAR * Total_Leaf_Area_Index)) + \
                                                   (1 - soil_nir_reflection) * (NIR_Direct * np.exp(-Beam_Reflection_Coefficient_NIR * Total_Leaf_Area_Index) + NIR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_NIR * Total_Leaf_Area_Index))
                # Calculate potential evaporation and net radiation using Penman-Monteith equation
                potential_evaporation, soil_absorbed_radiation = Leaf.Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, Hourly_temperature, vapor_pressure_deficit, slope_vapor_pressure_temperature, vapor_pressure_deficit)
                #print(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil,)
                
                # Proportional transpiration
                potential_transpiration_from_evaporative_soil_layer  = (shaded_transpiration + sunlit_transpiration) * self.Soil_Evaporative_Depth / self.Root_Depth
                #print(potential_transpiration_from_evaporative_soil_layer)    
    
               # Daytime course of water supply
                water_supply = self.Total_Soil_Water_Content_Top_Layer * (Hourly_sin_beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
                proportional_water_supply = water_supply * self.Soil_Evaporative_Depth / self.Root_Depth
                
                potential_evaporation_soil = max(0, potential_evaporation)
                actual_evaporation_soil = min(potential_evaporation_soil, potential_evaporation_soil / (potential_transpiration_from_evaporative_soil_layer + potential_evaporation_soil) * proportional_water_supply)
                #print(potential_evaporation_soil, potential_transpiration_from_evaporative_soil_layer, proportional_water_supply)

                
                # Temperature difference driven by soil evaporation and net radiation
                temperature_difference_soil = self.Limit_Function(-25., 25., (soil_absorbed_radiation - Latent_Heat_of_Water_Vaporization * actual_evaporation_soil) * (boundary_layer_resistance_soil + turbulence_resistance_soil) / Volumetric_Heat_Capacity_Air) 
                #print(soil_absorbed_radiation ,  actual_evaporation_soil, boundary_layer_resistance_soil , turbulence_resistance_soil)
                average_soil_temperature = Hourly_temperature + temperature_difference_soil
                
                # Recalculate potential evaporation with updated soil temperature
                saturation_vapor_pressure_soil = 0.611 * math.exp(17.4 * average_soil_temperature / (average_soil_temperature + 239.))
                slope_vapor_pressure_curve_soil = (saturation_vapor_pressure_soil - saturation_vapor_pressure) / self.Avoid_Zero_Division(temperature_difference_soil)
                potential_evaporation_second_estimate, soil_absorbed_radiation_second_estimate = Leaf.Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, average_soil_temperature, Vapour_Pressure_Deficit, slope_vapor_pressure_curve_soil, vapor_pressure_deficit)
                potential_soil_evaporation = max(0, potential_evaporation_second_estimate)
                #print(average_soil_temperature)
                
                Hourly_Soil_Evap.append(potential_soil_evaporation)
                Hourly_Soil_Rad.append(soil_absorbed_radiation_second_estimate)
                Hourly_Air_Soil_temp_dif.append(temperature_difference_soil)
                Hourly_boundary_layer_resistance_soil.append(boundary_layer_resistance_soil)
                Hourly_turbulence_resistance_soil.append(turbulence_resistance_soil)

    
    
                
            self.Hourly_Soil_Evap=Hourly_Soil_Evap
            self.Hourly_Soil_Rad=Hourly_Soil_Rad
            self.Hourly_Air_Soil_temp_dif=Hourly_Air_Soil_temp_dif
            self.Hourly_boundary_layer_resistance_soil=Hourly_boundary_layer_resistance_soil
            self.Hourly_turbulence_resistance_soil=Hourly_turbulence_resistance_soil
            
            # Aggregate Hourly data back to daily totals
            potential_evap_daily = Leaf.Leaf.aggregate_to_daily(Hourly_Soil_Evap,  Day_Length)
            potential_SoilRad_daily = Leaf.Leaf.aggregate_to_daily(Hourly_Soil_Rad,  Day_Length)
            self.potential_evap_daily=potential_evap_daily    
            self.potential_SoilRad_daily=potential_SoilRad_daily    
            self.Day_Air_Soil_temp_dif=(Hourly_Air_Soil_temp_dif * wgauss).sum()
    

    def Update_Evaporation_if_WaterStress(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, daily_water_supply, Soil_Evaporative_Depth, Root_Depth, Hourly_transpiration_Sunlit, Hourly_transpiration_Shaded, Hourly_Soil_Evap):
        gaussian_points = 5
        gaussian_weights_x = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
        gaussian_weights_w = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
        Latent_Heat_of_Vaporization = 2.4E6  # J/kg
        Volumetric_Heat_Capacity_Air = 1200  # J/m^3/°C

        Hourly_Actual_Soil_Evap = []
        Hourly_Air_Soil_Temperature_Difference = []

        for i, evap_potential, transp_sunlit, transp_shaded, Hourly_Soil_Rad, layer_resistance, turbulence_resistance in zip(range(gaussian_points), Hourly_Soil_Evap, Hourly_transpiration_Sunlit, Hourly_transpiration_Shaded, self.Hourly_Soil_Rad, self.Hourly_boundary_layer_resistance_soil, self.Hourly_turbulence_resistance_soil):
            hour = 12 - 0.5 * Day_Length + Day_Length * gaussian_weights_x[i]
            Sin_Beam = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (hour - 12.) / 24.))

            # Diurnal availability of soil water supply
            water_supply_Hourly = daily_water_supply * (Sin_Beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
            #print(water_supply_Hourly)
            # Available water for evaporation from the top soil layer
            water_supply_for_evaporation = water_supply_Hourly * (Soil_Evaporative_Depth / Root_Depth)

            # Total canopy transpiration (sunlit + shaded)
            total_canopy_transpiration = transp_sunlit + transp_shaded

            # Transpiration from the top soil layer
            top_soil_layer_transpiration = total_canopy_transpiration * (Soil_Evaporative_Depth / Root_Depth)
            #print(top_soil_layer_transpiration)

            # Maximum possible evaporation
            Max_possible_evap = evap_potential / (evap_potential + top_soil_layer_transpiration) * water_supply_for_evaporation
            Actual_Soil_Evap = Max_possible_evap if Max_possible_evap < evap_potential else evap_potential
            #print(evap_potential)
            
            Actual_Air_Soil_Temperature_Difference = self.Limit_Function(-25., 25., (Hourly_Soil_Rad - Latent_Heat_of_Vaporization * Actual_Soil_Evap) * (layer_resistance + turbulence_resistance) / Volumetric_Heat_Capacity_Air)
            Hourly_Air_Soil_Temperature_Difference.append(Actual_Air_Soil_Temperature_Difference)
            Hourly_Actual_Soil_Evap.append(Actual_Soil_Evap)
            #print(Actual_Air_Soil_Temperature_Difference)

        self.Hourly_Actual_Soil_Evap = Hourly_Actual_Soil_Evap
        # Aggregation to daily air-soil temperature difference and actual daily evaporation
        self.Day_Air_Soil_temp_dif = (Hourly_Air_Soil_Temperature_Difference* gaussian_weights_w).sum()
        # Assuming `aggregate_to_daily` is a method that sums or averages Hourly data to a daily total
        self.Actual_Daily_Evaporation = Leaf.Leaf.aggregate_to_daily(Hourly_Actual_Soil_Evap, Day_Length)
        #print(self.Day_Air_Soil_temp_dif)


    def Calculate_Soil_Temperature(self, tmin, tmax):
        # Calculate the daily average temperature with a weighted average favoring tmax
        daily_avg_temperature = 0.29 * tmin + 0.71 * tmax
        
        # Calculate the nightly average temperature with a weighted average favoring tmin
        nightly_avg_temperature = 0.71 * tmin + 0.29 * tmax
        
        # Calculate the soil's average steady state temperature
        soil_steady_state_avg_temp = (daily_avg_temperature + self.Day_Air_Soil_temp_dif + nightly_avg_temperature) / 2.
        
        # Calculate the rate of change in soil temperature
        rate_of_change_soil_temp = (soil_steady_state_avg_temp - self.average_soil_temperature) / self.Soil_Dynamic_Temperature_Factor

        
        # Print the rate of change in soil temperature for debugging or information purposes
        print("rate_of_change_soil_temp")
        print(self.rate_of_change_soil_temp)
        
        # Update class attributes with the new calculated values
        self.soil_steady_state_avg_temp = soil_steady_state_avg_temp
        self.rate_of_change_soil_temp = rate_of_change_soil_temp






    def Soil_Water_Components(self,Rainfall,Root_Depth,Actual_canopy_transpiration):
        
        
        if self.Current_Soil_Water_Content_Top_Layer <= self.Theoritical_Soil_Water_Content_Top_Layer:
            Drainage_From_Top_Layer = 0
        else:
           Drainage_From_Top_Layer = (self.Current_Soil_Water_Content_Top_Layer - self.Theoritical_Soil_Water_Content_Top_Layer) * self.Drainage_Factor

        
        if self.Current_Soil_Water_Content_Whole_Depth <= self.Theoritical_Soil_Water_Content_Top_Layer:
            Drainage_Overall = 0
        else:
           Drainage_Overall = (self.Current_Soil_Water_Content_Whole_Depth - self.Theoritical_Soil_Water_Content_Top_Layer) * self.Drainage_Factor


        
        Runof = 0
        
        
        Current_Soil_Water_Content_Top_Layer = self.Current_Soil_Water_Content_Top_Layer + Rainfall  - Drainage_From_Top_Layer - Runof - Actual_canopy_transpiration - self.Actual_Daily_Evaporation
        if self.Current_Soil_Water_Content_Top_Layer < 0:
            Current_Soil_Water_Content_Top_Layer = 0
         


        self.Current_Soil_Water_Content_Whole_Depth = self.Current_Soil_Water_Content_Whole_Depth + Rainfall  - Drainage_Overall - Runof - Actual_canopy_transpiration - self.Actual_Daily_Evaporation
        if self.Current_Soil_Water_Content_Whole_Depth < 0:
            self.Current_Soil_Water_Content_Whole_Depth = 0
        self.Theoritical_Soil_Water_Content_Whole_Depth = Root_Depth * (self.Field_Capacity-self.Residual_Soil_Moisture)
        self.Total_Soil_Water_Content_Whole_Depth = self.Current_Soil_Water_Content_Whole_Depth / self.Theoritical_Soil_Water_Content_Top_Layer
        
        
        Current_Soil_Moisture_Top_Layer=(Current_Soil_Water_Content_Top_Layer/self.Top_Layer_Depth)-self.Residual_Soil_Moisture
        Soil_Moisture_Fraction_Top_Layer=Current_Soil_Water_Content_Top_Layer/self.Theoritical_Soil_Water_Content_Top_Layer
        # Total_Soil_Water_Content_Top_Layer = (self.Residual_Soil_Moisture * self.Top_Layer_Depth) + Current_Soil_Water_Content_Top_Layer
        
        self.Current_Soil_Moisture_Top_Layer=Current_Soil_Moisture_Top_Layer
        self.Current_Soil_Water_Content_Top_Layer=Current_Soil_Water_Content_Top_Layer
        self.Soil_Moisture_Fraction_Top_Layer=Soil_Moisture_Fraction_Top_Layer
       
        
        
        
    def Calculate_Soil_N_Dynamics(self,Days_after_planting):

        soil_mass = self.Top_Layer_Depth * self.soil_bulk_density * (1 - self.fraction_soil_greater_than_2mm) * 1000  # g.m-2
        organic_N_mass = self.organic_N_percentage * 0.01 * soil_mass  # g.m-2
        mineralizable_organic_N = organic_N_mass * self.fraction_N_for_mineralization  # gN.m-2 org. N avail. for miner.
        nitrate_mass = self.nitrate_concentration_ppm * (14 / 62) * 0.000001 * soil_mass  # from ppm to gN.m-2
        ammonium_mass = self.ammonium_concentration_ppm * (14 / 18) * 0.000001 * soil_mass  # from ppm to gN.m-2
        self.soluble_N = nitrate_mass + ammonium_mass  # gN.m-2
        self.N_concentration = self.soluble_N / (self.Total_Soil_Water_Content_Top_Layer * 1000)  # gN.g-1 H2O
        
        
        N_mineralization_rate_constant = 24 * (math.exp(17.753 - 6350.5 / (self.average_soil_temperature + 273))) / 168
        N_mineralization_factor = 1 - math.exp(-N_mineralization_rate_constant)
        
        # Adjust for soil water content
        if self.Soil_Moisture_Fraction_Top_Layer < 0.9:
            mineralization_sensitivity = 1.111 * self.Soil_Moisture_Fraction_Top_Layer
        else:
            mineralization_sensitivity = 10 - 10 * self.Soil_Moisture_Fraction_Top_Layer
        if mineralization_sensitivity < 0:
            mineralization_sensitivity = 0

        
        # Calculate net N mineralization
        net_N_mineralization = mineralizable_organic_N * mineralization_sensitivity * N_mineralization_factor
        # Apply threshold of 200 mgN.L-1
        net_N_mineralization = net_N_mineralization * (0.0002 - self.N_concentration) / 0.0002
        if net_N_mineralization < 0:
            net_N_mineralization = 0
        
       
        
        if Days_after_planting == self.Fertilizer_applications_DAP[0]:  # Check if today matches any fertilization day
            
            N_fertilization = self.Fertilizer_applications_amount[0]  # gN.m-2, N to be applied
            volatilization_factor = self.Fraction_volatilization[0] / 100  # Convert percentage to decimal
            
            self.Fertilizer_applications_DAP=self.Fertilizer_applications_DAP[1:]
            self.Fertilizer_applications_amount=self.Fertilizer_applications_amount[1:]
            self.Fraction_volatilization=self.Fraction_volatilization[1:]
       
            
            N_volatilization = volatilization_factor * N_fertilization  # gN.m-2, N volatilized
       
        else:
           
            N_fertilization = N_volatilization = 0
             
        
         
     
        # Calculate N leaching based on soluble N, water, and drainage
        #Assuming no water drainage otherwise         N_leaching = soluble_N * (DRAIN1 / (water + DRAIN1))  # gN.m-2
        N_leaching = self.soluble_N 
        
        # Apply a threshold of 1 mgN.L-1 to the N concentration
        if self.N_concentration <= 0.000001:
            N_leaching = 0
        
    
     
             
        N_denitrification=0
        if self.Current_Soil_Moisture_Top_Layer >=  self.Saturated_Soil_Moisture:  # Condition for water saturation
            adjusted_N_concentration = self.N_concentration
            
            # Apply threshold of 400 mgN.L-1 to the N concentration
            if adjusted_N_concentration > 0.0004:
                adjusted_N_concentration = 0.0004
        
            # Calculate denitrification rate constant based on temperature
            denitrification_rate_constant = 6 * math.exp(0.07735 * self.soil_temperature - 6.593)
            
            # Calculate N denitrification per g of water
            N_denitrification = adjusted_N_concentration * (1 - math.exp(-denitrification_rate_constant))
            
            # Convert N denitrification to gN.m-2 based on water volume
            N_denitrification = N_denitrification * self.Total_Soil_Water_Content_Top_Layer * 1000
                        
            
        # Calculate the fraction of the top layer with roots
        fraction_top_layer_with_roots = self.Root_Depth / self.Top_Layer_Depth
        if fraction_top_layer_with_roots > 1:
            fraction_top_layer_with_roots = 1
        
        # Update soluble N content
        self.soluble_N = (self.soluble_N + net_N_mineralization + 
                            N_fertilization - N_volatilization -
                            N_leaching - N_denitrification - self.Cumulative_Nitrogen_uptake)
        
        # Recalculate N concentration in gN.g-1 H2O
        self.N_concentration = self.soluble_N / (self.Total_Soil_Water_Content_Top_Layer * 1000)
        
        # Calculate available soluble N in the root zone based on a threshold of 1 mgN.L-1
        available_soluble_N = (self.N_concentration - 0.000001) * self.Current_Soil_Water_Content_Top_Layer * 1000 * fraction_top_layer_with_roots
        if available_soluble_N < 0:
            available_soluble_N = 0
        
        self.available_soluble_N=available_soluble_N
        

    
    

    def Calculate_Nitrogen_Uptake(self, nitrogen_demand):

        
        # Calculate nitrogen uptake, respecting the demand and fixation rate
        Nitrogen_uptake = min(self.available_soluble_N,  max(0, nitrogen_demand))


        # Update class attributes with the calculated values
        self.Nitrogen_uptake = Nitrogen_uptake



    
        
        
        
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
        self.Cumulative_Nitrogen_uptake += self.Nitrogen_uptake
        
        # Update soil temperature with rate of temperature change
        self.average_soil_temperature += self.rate_of_change_soil_temp
    



     