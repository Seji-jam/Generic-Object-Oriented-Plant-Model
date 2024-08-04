import numpy as np
import math
from collections import namedtuple
import Leaf
import pandas as pd
def SWITCH_FUN(x, y1, y2):
    return y1 if x < 0 else y2
wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
                
Latent_Heat_of_Water_Vaporization = 2.4E6  # Latent heat of vaporization (J/kg)
Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity of air (J/m^3/°C)
class Soil:
    def __init__(self,Soil_Layer_Property, planting_depth, num_intervals, number_of_layers):
        
        self.Soil_Layer_Property=Soil_Layer_Property
        # self.Residual_Soil_Moisture = self.soil_parameters['Residual_Soil_Moisture']
        # self.Saturated_Soil_Moisture= self.soil_parameters['Saturated_Soil_Moisture']
        # self.field_capacity_water_content = self.soil_parameters['Field_Capacity']  # Water content at field capacity
        
        
        total_water_FC_profile=0
        total_water_S_profile=0
        total_water_PWP_profile=0
        total_depth=0
        for j, layer in enumerate(Soil_Layer_Property['layers']):

            total_depth+=layer['Layer_Thickness']

            total_water_FC_profile+=layer['Field_Capacity_']*layer['Layer_Thickness']
            self.Field_Capacity_soil_profile=total_water_FC_profile/total_depth

            total_water_S_profile+=layer['Saturated_Soil_Moisture']*layer['Layer_Thickness']
            self.Saturated_Soil_Moisture_soil_profile=total_water_S_profile/total_depth
            
            total_water_PWP_profile+=layer['Permanent_Wilting_Point']*layer['Layer_Thickness']
            self.Permanent_Wilting_Point_soil_profile=total_water_PWP_profile/total_depth
        self.total_depth=total_depth
        self.clay_percentage = 23.4  # Clay content in soil

        self.Evaporative_Depth=Soil_Layer_Property['layers'][0]['Layer_Thickness']
        self.planting_depth=planting_depth
        self.num_intervals=num_intervals
        self.number_of_layers=number_of_layers

        self.temperature_change_constant = 1
        self.temperature_change_T = 4
        self.initial_soil_temp = 15  # Initial soil temperature
        self.initial_decomposable_plant_material = 0  # Initial decomposable plant material in soil
        self.residual_ammonium = 1  # Residual ammonium-N in the soil
        self.residual_nitrate = 1  # Residual nitrate-N in the soil
        self.dpm_decomposition_rate = 10  # Decomposition rate of plant material (dpm rate)
        self.rpm_decomposition_rate = 0.3  # Decomposition rate of resistant plant material (rpm rate)
        self.Microbial_Biomass_decomposition_rate = 0.66  # Biomass incorporation rate into the soil organic matter pool
        self.humification_rate = 0.02  # Humification rate of organic matter
        self.total_organic_carbon = 7193.0  # Total organic carbon in the soil
        self.Microbial_Biomass_carbon_content = 3500  # Microbial_Biomass carbon content in the soil
        self.Microbial_Biomass_conversion_fraction = 0.03  # Fraction of biomass converted to Microbial_Biomass
        self.parameter_adjustment_multiplier = 1  # Multiplier factor for adjusting model parameters
        self.initial_ammonium_content = 2  # Initial Ammonium nitrogen content in soil
        self.initial_nitrate_content = 2  # Initial Nitrate nitrogen content in soil
        self.nitrogen_stress_water_index = -1  # Nitrogen stress water index
        self.water_supply_switch = -1  # Switch variable for water supply to crop
        # self.daily_water_input = 15  # Daily water input (precipitation + irrigation)
        self.ammonium_nitrogen_input_rate = 0  # Ammonium nitrogen input rate
        self.nitrate_nitrogen_input_rate = 0.65  # Nitrate nitrogen input rate
        self.plant_material_dpm_rpm_ratio = 1.44  # Ratio dpm/rpm of added plant material
        self.soil_resistance_to_evaporation = 100  # Soil resistance to evaporation

        
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
        self.nitrogen_stress_water_index = 0  # Factor indicating water stress
        self.nitrogen_uptake = 0  # Nitrogen uptake rate


        # self.hourly_Soil_Evap=[]
        self.hourly_Soil_Rad=[]
        self.Hourly_Actual_Soil_Evap=[]

        # # Irrigation and fertilization can be input in the weather file as well! 
        # # Irrigation schedule: Days and corresponding irrigation amount (currently set to 0)
        # self.irrigation_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), 
        #                             (55.0, 0), (65.0, 0), (75.0, 0), (85.0, 0), (95.0, 0)]
        
        # Fertilization schedules for ammonium  and nitrate nitrogen (currently set to 0)
        self.ammonium_nitrogen_application_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), 
                                                       (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]
        self.nitrate_nitrogen_application_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), 
                                                      (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]

        self.deniul = 0  # Denitrification in the upper layer
        self.denill = 0
        self.nuptn = 0  # Nitrogen uptake in the nitrogen form
        self.nupta= 0


        # self.soil_parameters['current_soil_moisture'] = self.soil_parameters['initial_soil_moisture']
        # self.evaporation_depth=self.soil_parameters['evaporation_depth']
        # self.total_depth=self.soil_parameters['number_of_soil_layers']* self.soil_parameters['layer_depth']
        # self.Top_Layer_Depth = self.soil_parameters['layer_depth'][0]*100  # 


        # Initial conditions and parameter calculations
        root_depth_ini = 5 #cm
        Microbial_Biomass_initial_content =  self.Microbial_Biomass_conversion_fraction *  self.total_organic_carbon
        # soil_moisture_initial =  self.field_capacity_water_content *  self.parameter_adjustment_multiplier
       
        resistant_plant_material_initial =  self.total_organic_carbon -  Microbial_Biomass_initial_content -  self.initial_decomposable_plant_material
        humus_initial_content = Microbial_Biomass_initial_content - Microbial_Biomass_initial_content
        decomposable_plant_nitrogen_initial = 1.0 / 40.0 *  self.initial_decomposable_plant_material
        ammonium_upper_layer_initial = (1.0 - np.exp(-0.065 * root_depth_ini)) *  self.initial_ammonium_content + root_depth_ini / self.total_depth *  self.residual_ammonium
        ammonium_lower_layer_initial = np.exp(-0.065 * root_depth_ini) *  self.initial_ammonium_content + (1.0 - root_depth_ini / self.total_depth) *  self.residual_ammonium
        nitrate_upper_layer_initial = (1.0 - np.exp(-0.065 * root_depth_ini)) *  self.initial_nitrate_content + root_depth_ini / self.total_depth *  self.residual_nitrate
        nitrate_lower_layer_initial = np.exp(-0.065 * root_depth_ini) *  self.initial_nitrate_content + (1.0 - root_depth_ini / self.total_depth) *  self.residual_nitrate
        resistant_plant_nitrogen_initial = 1.0 / 100.0 * resistant_plant_material_initial

        self.average_root_zone_moisture=0
        self.average_below_root_zone_moisture=0
        self.average_root_zone_water_content=0
        self.average_below_root_zone_moisture=0
        
        self.water_supply_for_evaporation=0
        self.water_supply_for_Transpiration=0
        self.average_soil_moisture_soil_profile=0
        
        
        
        self.groundwater_recharge=0
        self.ammonium_upper_layer=ammonium_upper_layer_initial
        self.ammonium_lower_layer=ammonium_lower_layer_initial
        self.nitrate_lower_layer=nitrate_lower_layer_initial
        self.nitrate_upper_layer=nitrate_upper_layer_initial
        # self.water_content_upper_layer_ChangeRate=0
        # self.water_content_lower_layer_ChangeRate=0
        self.rain_and_irrigation=0
        self.recharge_to_upper_layer=0
        
        

        self.average_soil_temperature = 0
        self.rate_of_change_soil_temp = 0
    
        self.water_content_upper_limit = 0
        self.water_content_lower_limit = 0
        self.nitrate_nitrogen_supply = 0
        self.ammonium_fertilizer_application = 0
        self.nitrate_fertilizer_application = 0

        self.ammonia_volatilization = 0
        # State variables
        self.total_nitrogen_uptake = 0
        # self.water_content_upper_layer=self.water_upper_soil_initial
        # self.water_content_lower_layer=self.water_lower_soil_initial

        self.decomposable_plant_material =  self.initial_decomposable_plant_material  
        self.resistant_plant_material = resistant_plant_material_initial  
        self.Microbial_Biomass_content = Microbial_Biomass_initial_content  
        self.humus_content = humus_initial_content  
        self.decomposable_plant_nitrogen = decomposable_plant_nitrogen_initial  
        self.resistant_plant_nitrogen = resistant_plant_nitrogen_initial   
        # self.water_upper_soil = water_upper_soil_initial 
        # self.water_lower_soil = water_lower_soil_initial 
        # Root_Depth = root_depth_ini  
        self.Soil_Temperature =  self.initial_soil_temp 
        self.ammonium_upper_layer = ammonium_upper_layer_initial
        self.nitrate_upper_layer = nitrate_upper_layer_initial
        self.nitrate_lower_layer = nitrate_lower_layer_initial
        self.ammonium_lower_layer = ammonium_lower_layer_initial
        self.total_nitrogen_leached = 0
        self.N_volatilization = 0
        self.soil_moisure_evaporative_layer = 0

        self.potential_SoilRad_daily = 0
        self.Day_Air_Soil_temp_dif  = 0
        self.Actual_Daily_Evaporation =0
        self.hourly_rbhs=[]
        self.hourly_rts=[]

    def switch_function(self,x, y1, y2):
        return y1 if x < 0 else y2
        
    def Function_Adjust_Switch(self,x, y1, y2, y3):
        """
        Input switch. Returns y1 if x < 0, y2 if x = 0, and y3 if x > 0.
    
        Parameters:
        x (float): Control variable
        y1 (float): Returned value if x < 0
        y2 (float): Returned value if x = 0
        y3 (float): Returned value if x > 0
    
        Returns:
        float: The value of y based on the conditions.
        """
        if x < 0:
            return y1
        elif x == 0:
            return y2
        else:
            return y3

    def van_genuchten(self,soil_moisture, alpha, n, Residual_Soil_Moisture, Saturated_Soil_Moisture):
        m = 1 - 1 / n
        effective_saturation = (soil_moisture - Residual_Soil_Moisture) / (Saturated_Soil_Moisture - Residual_Soil_Moisture)
        effective_saturation = np.clip(effective_saturation, 0.0001, 1)  # prevent division by zero
        pressure_head = ((effective_saturation**(-1/m) - 1)**(1/n)) / alpha
        return pressure_head
    
    def hydraulic_conductivity(self,soil_moisture, Saturated_Hydraulic_Conductivity, alpha, n, Residual_Soil_Moisture, Saturated_Soil_Moisture):
        m = 1 - 1 / n
        effective_saturation = (soil_moisture - Residual_Soil_Moisture) / (Saturated_Soil_Moisture - Residual_Soil_Moisture)
        effective_saturation = np.clip(effective_saturation, 0.0001, 1)  # prevent division by zero
        conductivity = Saturated_Hydraulic_Conductivity * (effective_saturation**0.5) * (1 - (1 - effective_saturation**(1/m))**m)**2
        return conductivity
 
    # def Calculate_Soil_Water_Content(self,irrigation_amount,Root_Depth):

    #     number_of_layers = self.soil_parameters['number_of_soil_layers']   # number of soil layers
    #     layer_depth = self.soil_parameters['layer_depth']  # depth of each soil layer in cm
    #     # time_step = 1.0  # time step in days

    #     # # Soil parameters (Van Genuchten parameters)
    #     # alpha = self.soil_parameters['alpha']  # (1/cm)
    #     # n = self.soil_parameters['n']
    #     Residual_Soil_Moisture = self.soil_parameters['Residual_Soil_Moisture']  # residual soil moisture
    #     Saturated_Soil_Moisture = self.soil_parameters['Saturated_Soil_Moisture']   # saturated soil moisture
    #     # Saturated_Hydraulic_Conductivity = self.soil_parameters['Saturated_Hydraulic_Conductivity']  # saturated hydraulic conductivity (cm/day)


    #     # # Irrigation and evapotranspiration rates (cm/day)
    #     # irrigation_rate =irrigation_amount  # cm/day (only applied to the top layer)
   
        
    #     # # Initialize soil moisture array
    #     soil_moisture = np.zeros((1, number_of_layers))
    #     soil_moisture[0, :] = np.full( number_of_layers,self.soil_parameters['current_soil_moisture'])

    #     # timesteps=2
    #     # # Richards equation solver
    #     # for t in range(1, timesteps):
    #     #     conductivity = np.zeros(number_of_layers)
    #     #     pressure_head = np.zeros(number_of_layers)
    #     #     for z in range(number_of_layers):
    #     #         conductivity[z] = self.hydraulic_conductivity(soil_moisture[t-1, z], Saturated_Hydraulic_Conductivity, alpha, n, Residual_Soil_Moisture, Saturated_Soil_Moisture)
    #     #         pressure_head[z] = self.van_genuchten(soil_moisture[t-1, z], alpha, n, Residual_Soil_Moisture, Saturated_Soil_Moisture)
            
    #     #     # Calculate fluxes
    #     #     flux = np.zeros(number_of_layers + 1)
    #     #     flux[0] = irrigation_rate  # Irrigation at the surface
    #     #     for z in range(1, number_of_layers):
    #     #         flux[z] = -conductivity[z-1] * (pressure_head[z-1] - pressure_head[z]) / layer_depth
    #     #     flux[number_of_layers] = conductivity[number_of_layers-1] * (pressure_head[number_of_layers-1] - 0) / layer_depth  # Free drainage at the bottom
            
    #     #     # Update soil moisture
    #     #     for z in range(number_of_layers):
    #     #         soil_moisture[t, z] = soil_moisture[t-1, z] + (time_step / layer_depth) * (flux[z] - flux[z+1])

            
    #     #     # Apply boundary conditions
    #     #     soil_moisture[t, :] = np.clip(soil_moisture[t, :], Residual_Soil_Moisture, Saturated_Soil_Moisture)
    #     # # self.soil_parameters['current_soil_moisture']=soil_moisture[-1, :]



    #     water_stored_evaporation_depth = 0
    #     for z in range(number_of_layers):
    #         if z * layer_depth < self.evaporation_depth:
    #             effective_depth = min(self.evaporation_depth - z * layer_depth, layer_depth)
    #             # print(self.evaporation_depth)
    #             water_stored_evaporation_depth += (soil_moisture[0, z]-Residual_Soil_Moisture) * (effective_depth / layer_depth) * layer_depth * 1e-2  # Convert cm to m for volume
    #     self.soil_moisure_evaporative_layer=100* water_stored_evaporation_depth/self.evaporation_depth
    #     # Ensure the new soil moisture does not go below residual soil moisture
    #     # self.soil_parameters['current_soil_moisture'] = np.maximum(soil_moisture[-1, :], self.soil_parameters['Residual_Soil_Moisture'])
    #     # print(self.soil_parameters['current_soil_moisture'])
    #     self.water_supply_for_evaporation=water_stored_evaporation_depth*1000 #mm
    #     # print("258",self.water_supply_for_evaporation,soil_moisture)

    #     water_stored_root_zone = 0
    #     water_stored_below_root_zone = 0
    #     for z in range(number_of_layers):
    #         if (z) * layer_depth < Root_Depth:
    #             effective_depth = min(Root_Depth - z * layer_depth, layer_depth)
    #             water_stored_root_zone += (soil_moisture[0, z]-Residual_Soil_Moisture) * (effective_depth / layer_depth) * layer_depth * 1e-2  # Convert cm to m for volume
    #             water_stored_below_root_zone += (soil_moisture[0, z]-Residual_Soil_Moisture) * (layer_depth-effective_depth) * 1e-2  # Convert cm to m for volume
    #         else:
    #             water_stored_below_root_zone += (soil_moisture[0, z]-Residual_Soil_Moisture) * layer_depth * 1e-2  # Convert cm to m for volume
        
                
    #     self.average_root_zone_water_content=water_stored_root_zone #m
    #     self.average_below_root_zone_water_content=water_stored_below_root_zone #m

    #     # print (self.average_root_zone_water_content*1000)

    
    def soil_water_balance(self, duration,Soil_Layer_Property, planting_depth, num_intervals, number_of_layers, root_depth,   rain, potential_transpiration_canopy, potential_evaporation_soil):
        # Define namedtuples
        Soil_Water_Content = namedtuple('Soil_Water_Content', 'Saturated_Soil_Moisture Field_Capacity_ Permanent_Wilting_Point Pore_Size_Distribution Air_Entry')
        # RootZone = namedtuple('RootZone', 'Water_Content Soil_Moisture critical Saturated_Soil_Moisture Field_Capacity_ Permanent_Wilting_Point')
        ActualET = namedtuple('ActualET', 'crop soil')
        Fluxes = namedtuple('Fluxes', 'transpiration evaporation influx outflux netflux')
    
        def initialize_layer(idx, prev_layer, next_layer):
            nonlocal SoilWater__accumulated_layer_depth
            layer = layers[idx]
            layer['prev'] = prev_layer
            layer['next'] = next_layer
    
            prev_accumulated_layer_depth = prev_layer['Accumulated_Layer_Depth'] if prev_layer else 0.0
            layer['Accumulated_Layer_Depth'] = layer['Layer_Thickness'] + prev_accumulated_layer_depth
            prev_thickness = prev_layer['Layer_Thickness'] if prev_layer else 0.0
            depth = 0.5 * (prev_thickness + layer['Layer_Thickness'])
            layer['depth'] = SoilWater__accumulated_layer_depth + depth
            SoilWater__accumulated_layer_depth += depth
    
            if layer['Soil_Moisture'] < 0:
                Soil_Moisture = -layer['Soil_Moisture']
                Field_Capacity_ = layer['Soil_Water_Content'].Field_Capacity_
                if 1 <= Soil_Moisture <= 2:
                    Saturated_Soil_Moisture = layer['Soil_Water_Content'].Saturated_Soil_Moisture
                    Soil_Moisture = Saturated_Soil_Moisture - (Soil_Moisture - 1) * (Saturated_Soil_Moisture - Field_Capacity_)
                elif 2 < Soil_Moisture <= 3:
                    Permanent_Wilting_Point = layer['Soil_Water_Content'].Permanent_Wilting_Point
                    Soil_Moisture = Field_Capacity_ - (Soil_Moisture - 2) * (Field_Capacity_ - Permanent_Wilting_Point)
                else:
                    Soil_Moisture = Field_Capacity_
                layer['Soil_Moisture'] = Soil_Moisture
                layer['Water_Content'] = layer['Soil_Moisture'] * layer['Layer_Thickness'] * 1000
    
            update_heads_k(layer)
    
        def update_heads_k(layer):
            Field_Capacity_ = layer['Soil_Water_Content'].Field_Capacity_
            Soil_Moisture = layer['Soil_Moisture']
            if Soil_Moisture >= Field_Capacity_:
                delta_Field_Capacity_ = Soil_Moisture - Field_Capacity_
                head_matric = 33 - (33 - layer['Soil_Water_Content'].Air_Entry) * delta_Field_Capacity_ / (layer['Soil_Water_Content'].Saturated_Soil_Moisture - Field_Capacity_)
                head_matric /= 10
            else:
                exponent_b = 1 / layer['Soil_Water_Content'].Pore_Size_Distribution
                a_coefficient = math.exp(3.496508 + exponent_b * math.log(Field_Capacity_))
                head_matric = (a_coefficient * max(0.05, Soil_Moisture) ** (-exponent_b)) / 10
            layer['matric_potential'] = max(0.0, head_matric)
            layer['gravity_potential'] = layer['depth']
            Air_Entry = layer['Soil_Water_Content'].Air_Entry / 10
            head_matric = layer['matric_potential']
            moisture_ratio = layer['Soil_Moisture'] / layer['Soil_Water_Content'].Saturated_Soil_Moisture
            if head_matric > Air_Entry:
                layer['Unaturated_hydraulic_conductivity'] = layer['Saturated_Hydraulic_Conductivity'] * moisture_ratio ** (3 + 2 / layer['Soil_Water_Content'].Pore_Size_Distribution)
            else:
                layer['Unaturated_hydraulic_conductivity'] = layer['Saturated_Hydraulic_Conductivity']
    
        def tothead(layer):
            return layer['matric_potential'] + layer['gravity_potential']
    
    
        def calculate_water_fluxes(cumulative_fluxes, potential_transpiration_canopy, potential_evaporation_soil, net_rain):
            # calculate_root_zone_water()
            nonlocal actual_evapotranspiration
            actual_evapotranspiration = ActualET(potential_transpiration_canopy, potential_evaporation_soil)
    
            previous_soil_moisture_potential = 0.0
            for index in range(number_of_layers):
                current_layer = layers[index]
                previous_layer = current_layer['prev']
    
                update_heads_k(current_layer)
    
                evaporation_influx = 0.0 if previous_layer is not None else actual_evapotranspiration.soil / 1000
                root_fraction = min(1.0, current_layer['Accumulated_Layer_Depth'] / root_depth)
                current_soil_moisture_potential = 1.8 * root_fraction - 0.8 * root_fraction ** 2
                transpiration_influx = actual_evapotranspiration.crop * (current_soil_moisture_potential - previous_soil_moisture_potential) / 1000
                previous_soil_moisture_potential = current_soil_moisture_potential
    
                if previous_layer is not None:
                    log_k_ratio = math.log(current_layer['Unaturated_hydraulic_conductivity']) - math.log(previous_layer['Unaturated_hydraulic_conductivity'])
                    k_average = (current_layer['Unaturated_hydraulic_conductivity'] - previous_layer['Unaturated_hydraulic_conductivity']) / log_k_ratio if log_k_ratio != 0.0 else current_layer['Unaturated_hydraulic_conductivity']
                    gradient = (tothead(current_layer) - tothead(previous_layer)) / (current_layer['depth'] - previous_layer['depth'])
                    current_influx = k_average * gradient - transpiration_influx
                else:
                    net_rain_mm = net_rain / 1000
                    current_influx = min(net_rain_mm, current_layer['Saturated_Hydraulic_Conductivity']) - evaporation_influx - transpiration_influx
    
                cumulative_fluxes[index]['transpiration'] = transpiration_influx
                cumulative_fluxes[index]['evaporation'] = evaporation_influx
                cumulative_fluxes[index]['influx'] = current_influx
    
            for index in range(number_of_layers):
                current_layer = layers[index]
                next_layer = current_layer['next']
    
                Water_Content = current_layer['Soil_Moisture'] * current_layer['Layer_Thickness']
                influx = cumulative_fluxes[index]['influx']
                if next_layer is not None:
                    outflux = cumulative_fluxes[index + 1]['influx']
                else:
                    outflux = current_layer['Unaturated_hydraulic_conductivity']
    
                next_Water_Content = influx + Water_Content - outflux
                dry_limit = current_layer['Layer_Thickness'] * 0.005
                saturated_limit = current_layer['Layer_Thickness'] * current_layer['Soil_Water_Content'].Saturated_Soil_Moisture
                if next_Water_Content < dry_limit:
                    outflux = influx + Water_Content - dry_limit
                elif next_Water_Content > saturated_limit:
                    outflux = influx + Water_Content - saturated_limit
    
                if next_layer is not None:
                    cumulative_fluxes[index + 1]['influx'] = outflux
    
                cumulative_fluxes[index]['outflux'] = outflux
                cumulative_fluxes[index]['netflux'] = cumulative_fluxes[index]['influx'] - cumulative_fluxes[index]['outflux']
    
                Water_Content += cumulative_fluxes[index]['netflux'] / num_intervals
                current_layer['Soil_Moisture'] = max(0.005, min(current_layer['Soil_Water_Content'].Saturated_Soil_Moisture, Water_Content / current_layer['Layer_Thickness']))
                current_layer['Water_Content'] = current_layer['Soil_Moisture'] * current_layer['Layer_Thickness'] * 1000
    
                for field in Fluxes._fields:
                    flux_increment = cumulative_fluxes[index][field] / num_intervals
                    cumulative_fluxes[index][field] += flux_increment
    
        def daily_water_balance(rain, potential_transpiration_canopy, potential_evaporation_soil):
            nonlocal net_rain, root_depth
            net_rain = rain
            root_depth = root_depth
    
            cumulative_fluxes = [{field: 0.0 for field in Fluxes._fields} for _ in range(number_of_layers)]
    
            for interval in range(num_intervals):
                calculate_water_fluxes(cumulative_fluxes, potential_transpiration_canopy, potential_evaporation_soil, net_rain)
    
            # root_water = RootZone(*[root_zone[field] for field in RootZone._fields])
            # print('root_water------------->',root_water)
            for index, layer in enumerate(layers):
                layer['fluxes'] = Fluxes(*[cumulative_fluxes[index][field] for field in Fluxes._fields])
    
        def calculate_root_zone_stat(layers, planting_depth, root_depth):
            # total_moisture = 0.0
            Root_zone_water_content = 0.0
            total_thickness = 0.0
            Root_zone_plant_avail_water=0.0
            for layer in layers:
    
                top_of_layer = layer['Accumulated_Layer_Depth'] - layer['Layer_Thickness']
                top_of_layer=int(round(top_of_layer,2)*100)/100
                # top_of_layer=top_of_layer
                # print(top_of_layer,layer['Layer_Thickness'],layer['Accumulated_Layer_Depth'])
                bottom_of_layer = layer['Accumulated_Layer_Depth']
                if top_of_layer >= planting_depth and bottom_of_layer <= planting_depth + root_depth:
                    # print('c1')
                    thickness_within_zone = layer['Layer_Thickness']
                elif top_of_layer < planting_depth and bottom_of_layer > planting_depth:
                    # print('c2')
                    # print(top_of_layer)
                    thickness_within_zone = bottom_of_layer - planting_depth
                elif top_of_layer < planting_depth + root_depth and bottom_of_layer > planting_depth + root_depth:
                    thickness_within_zone = planting_depth + root_depth - top_of_layer
                    # print('c3')
    
                else:
                    # print('c4')
    
                    continue
    
                Root_zone_water_content += layer['Soil_Moisture'] * thickness_within_zone
                Root_zone_plant_avail_water +=max(0,(layer['Soil_Moisture']-layer['Soil_Water_Content'].Permanent_Wilting_Point)* thickness_within_zone)
               
                total_thickness += thickness_within_zone
    
            if total_thickness > 0:
                avg_soil_moisture = Root_zone_water_content / total_thickness
                avg_water_content = Root_zone_water_content*1000
                Root_zone_plant_avail_water =Root_zone_plant_avail_water*1000
            else:
                avg_soil_moisture = 0.0
                avg_water_content = 0.0
                Root_zone_plant_avail_water =0.0
    
            return avg_soil_moisture, avg_water_content,Root_zone_plant_avail_water
    
        layers = []
        net_rain = 0.0
        actual_evapotranspiration = ActualET(0.0, 0.0)
        cumulative_fluxes = [{field: 0.0 for field in Fluxes._fields} for _ in range(number_of_layers)]
        # root_zone = {field: 0.0 for field in RootZone._fields}
        SoilWater__accumulated_layer_depth = 0.0
    
        layers_config = Soil_Layer_Property['layers']
        for i in range(number_of_layers):
            layer = {
                'Layer_Thickness': layers_config[i]['Layer_Thickness'],
                'Soil_Moisture': layers_config[i]['Soil_Moisture'],
                'Water_Content': layers_config[i]['Soil_Moisture'] * layers_config[i]['Layer_Thickness'] * 1000,
                'Soil_Water_Content': Soil_Water_Content(layers_config[i]['Saturated_Soil_Moisture'],
                                                         layers_config[i]['Field_Capacity_'],
                                                         layers_config[i]['Permanent_Wilting_Point'],
                                                         layers_config[i]['Pore_Size_Distribution'],
                                                         layers_config[i]['Air_Entry']),
                'Saturated_Hydraulic_Conductivity': layers_config[i]['Saturated_Hydraulic_Conductivity'],
                'Unaturated_hydraulic_conductivity': 0.0,
                'matric_potential': 0.0,
                'gravity_potential': 0.0,
                'fluxes': Fluxes(0.0, 0.0, 0.0, 0.0, 0.0),
                'prev': None,
                'next': None
            }
            layers.append(layer)
    
        for i in range(number_of_layers):
            prev_layer = layers[i - 1] if i > 0 else None
            next_layer = layers[i + 1] if i < number_of_layers - 1 else None
            initialize_layer(i, prev_layer, next_layer)
    
        # root_zone_stat=calculate_root_zone_water()
    
        for day in range(1, duration + 1):
            # print(f'Day {day}')
            daily_water_balance(rain, potential_transpiration_canopy, potential_evaporation_soil)
    
        soil_moisture_array = [layer['Soil_Moisture'] for layer in layers]
        water_content_array = [layer['Water_Content'] for layer in layers]
    
        Root_zone_avg_soil_moisture, Root_zone_avg_water_content,Root_zone_plant_avail_water = calculate_root_zone_stat(layers, planting_depth, root_depth)
    
        print('\ndone.')
        
        return soil_moisture_array, water_content_array, Root_zone_avg_soil_moisture, Root_zone_avg_water_content,Root_zone_plant_avail_water

    def Current_Soil_Water_Status(self, irrigation_amount,Root_Depth):
        soil_moisture_array, _, _, _,Root_zone_plant_avail_water=self.soil_water_balance(1,self.Soil_Layer_Property, self.planting_depth, self.num_intervals, self.number_of_layers, Root_Depth, irrigation_amount, 0, 0 )
        self.water_supply_for_evaporation=soil_moisture_array[0]*self.Evaporative_Depth
        self.water_supply_for_Transpiration=Root_zone_plant_avail_water
        
    def Calculate_Soil_Evaporation(self,irrigation_amount, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination,
                                             Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature,
                                             Min_Temperature, Vapour_Pressure_Deficit, Wind_Speed, 
                                             soil_resistance_to_evaporation,Root_Depth,Leaf_Blade_Angle, Total_Leaf_Area_Index,
                                             Light_Extinction_Coefficient, Potential_Canopy_Transpiration,Hourly_transpiration_Shaded, Hourly_transpiration_Sunlit):

            
            # Convert daily climate data to Hourly data
            Hourly_climate_data, gaussian_weights = Leaf.Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature, Min_Temperature, Wind_Speed)
            potential_evaporation_soil_first_estimate=[]
            for gaussian_weight,(Hourly_Solar_Constant, Hourly_temperature, Hourly_sin_beam, Hourly_Solar_Radiation, Hourly_wind_speed), shaded_transpiration, sunlit_transpiration in zip(gaussian_weights,Hourly_climate_data, Hourly_transpiration_Shaded, Hourly_transpiration_Sunlit):
                
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
                soil_nir_reflection = self.switch_function(self.soil_moisure_evaporative_layer - 0.5, 0.52 - 0.68 * self.soil_moisure_evaporative_layer, 0.18)
                absorbed_total_radiation_by_soil = (1 - soil_par_reflection) * (PAR_Direct * np.exp(-Beam_Reflection_Coefficient_PAR * Total_Leaf_Area_Index) + PAR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_PAR * Total_Leaf_Area_Index)) + \
                                                   (1 - soil_nir_reflection) * (NIR_Direct * np.exp(-Beam_Reflection_Coefficient_NIR * Total_Leaf_Area_Index) + NIR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_NIR * Total_Leaf_Area_Index))
                # Calculate potential evaporation and net radiation using Penman-Monteith equation
                potential_evaporation_soil, soil_absorbed_radiation = Leaf.Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, Hourly_temperature, vapor_pressure_deficit, slope_vapor_pressure_temperature, vapor_pressure_deficit)
                #print(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil,)
                potential_evaporation_soil = max(0, potential_evaporation_soil)
                potential_evaporation_soil_first_estimate.append(potential_evaporation_soil)
            potential_evaporation_soil_first_estimate_daily = Leaf.Leaf.aggregate_to_daily(potential_evaporation_soil_first_estimate, Day_Length)
            # print(potential_evaporation_soil_first_estimate)
            Maximum_possible_Evaporation=self.water_supply_for_evaporation
            if potential_evaporation_soil_first_estimate_daily > Maximum_possible_Evaporation :
                Water_Stress_Fraction=Maximum_possible_Evaporation/potential_evaporation_soil_first_estimate_daily
                potential_evaporation_soil_estimate=np.array(potential_evaporation_soil_first_estimate)*Water_Stress_Fraction
            else:
                potential_evaporation_soil_estimate=np.array(potential_evaporation_soil_first_estimate)
            # print('372',potential_evaporation_soil_estimate,Maximum_possible_Evaporation)

            Hourly_Actual_Soil_Evap = []
            Hourly_Air_Soil_temp_dif = []
            Hourly_Soil_Rad = []
            Hourly_boundary_layer_resistance_soil=[]
            Hourly_turbulence_resistance_soil=[]

            for potential_evaporation_soil_estimate_instant,(Hourly_Solar_Constant, Hourly_temperature, Hourly_sin_beam, Hourly_Solar_Radiation, Hourly_wind_speed), shaded_transpiration, sunlit_transpiration in zip(potential_evaporation_soil_estimate,Hourly_climate_data, Hourly_transpiration_Shaded, Hourly_transpiration_Sunlit):
                
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
                soil_nir_reflection = self.switch_function(self.soil_moisure_evaporative_layer - 0.5, 0.52 - 0.68 * self.soil_moisure_evaporative_layer, 0.18)
                absorbed_total_radiation_by_soil = (1 - soil_par_reflection) * (PAR_Direct * np.exp(-Beam_Reflection_Coefficient_PAR * Total_Leaf_Area_Index) + PAR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_PAR * Total_Leaf_Area_Index)) + \
                                                   (1 - soil_nir_reflection) * (NIR_Direct * np.exp(-Beam_Reflection_Coefficient_NIR * Total_Leaf_Area_Index) + NIR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_NIR * Total_Leaf_Area_Index))
               
           
                
                # Temperature difference driven by soil evaporation and net radiation
                temperature_difference_soil = self.Limit_Function(-25., 25., (soil_absorbed_radiation - Latent_Heat_of_Water_Vaporization * potential_evaporation_soil_estimate_instant) * (boundary_layer_resistance_soil + turbulence_resistance_soil) / Volumetric_Heat_Capacity_Air) 
                #print(soil_absorbed_radiation ,  actual_evaporation_soil, boundary_layer_resistance_soil , turbulence_resistance_soil)
                average_soil_temperature = Hourly_temperature + temperature_difference_soil
                
                # Recalculate potential evaporation with updated soil temperature
                saturation_vapor_pressure_soil = 0.611 * math.exp(17.4 * average_soil_temperature / (average_soil_temperature + 239.))
                slope_vapor_pressure_curve_soil = (saturation_vapor_pressure_soil - saturation_vapor_pressure) / self.Avoid_Zero_Division(temperature_difference_soil)
                potential_evaporation_second_estimate, soil_absorbed_radiation_second_estimate = Leaf.Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, average_soil_temperature, Vapour_Pressure_Deficit, slope_vapor_pressure_curve_soil, vapor_pressure_deficit)
                Actual_Soil_Evap = max(0, potential_evaporation_second_estimate)
                
                
                Hourly_Actual_Soil_Evap.append(Actual_Soil_Evap)
                Hourly_Soil_Rad.append(soil_absorbed_radiation_second_estimate)
                Hourly_Air_Soil_temp_dif.append(temperature_difference_soil)
                Hourly_boundary_layer_resistance_soil.append(boundary_layer_resistance_soil)
                Hourly_turbulence_resistance_soil.append(turbulence_resistance_soil)
            Actual_Daily_Evaporation = Leaf.Leaf.aggregate_to_daily(Hourly_Actual_Soil_Evap, Day_Length)

            if Actual_Daily_Evaporation > Maximum_possible_Evaporation :
                Water_Stress_Fraction=Maximum_possible_Evaporation/Actual_Daily_Evaporation
                potential_evaporation_soil_estimate=np.array(Hourly_Actual_Soil_Evap)*Water_Stress_Fraction
                self.Hourly_Actual_Soil_Evap=potential_evaporation_soil_estimate
                self.Actual_Daily_Evaporation = Leaf.Leaf.aggregate_to_daily(potential_evaporation_soil_estimate, Day_Length)

                # if (Water_Stress_Fraction <0.95) and (Water_Stress_Fraction>0):
                #     print(Water_Stress_Fraction)
                #     print('warning ..... Soil evaporation convergence issue!!!!!!!!!')
            else:
                self.Hourly_Actual_Soil_Evap=Hourly_Actual_Soil_Evap
                self.Actual_Daily_Evaporation = Actual_Daily_Evaporation
            # print('484',self.Actual_Daily_Evaporation,Maximum_possible_Evaporation,irrigation_amount)

            # print(Hourly_Actual_Soil_Evap)    

            self.Hourly_Soil_Rad=Hourly_Soil_Rad
            self.Hourly_Air_Soil_temp_dif=Hourly_Air_Soil_temp_dif
            self.Hourly_boundary_layer_resistance_soil=Hourly_boundary_layer_resistance_soil
            self.Hourly_turbulence_resistance_soil=Hourly_turbulence_resistance_soil
            
            # Aggregate Hourly data back to daily totals
            self.potential_SoilRad_daily = Leaf.Leaf.aggregate_to_daily(Hourly_Soil_Rad,  Day_Length)
            self.Day_Air_Soil_temp_dif=(Hourly_Air_Soil_temp_dif * wgauss).sum()
    

  
    # def Update_Evaporation_if_WaterStress(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, daily_water_supply, Root_Depth, Hourly_transpiration_Sunlit, Hourly_transpiration_Shaded, Hourly_Soil_Evap):
    #     gaussian_points = 5
    #     gaussian_weights_x = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
    #     gaussian_weights_w = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
    #     Latent_Heat_of_Vaporization = 2.4E6  # J/kg
    #     Volumetric_Heat_Capacity_Air = 1200  # J/m^3/°C

    #     Hourly_Actual_Soil_Evap = []
    #     Hourly_Air_Soil_Temperature_Difference = []
    #     # water_supply_for_evaporation = self.water_supply_for_evaporation
    #     # print('next')
    #     for i, potential_evaporation_soil, sunlit_transpiration, shaded_transpiration, Hourly_Soil_Rad, layer_resistance, turbulence_resistance in zip(range(gaussian_points), Hourly_Soil_Evap, Hourly_transpiration_Sunlit, Hourly_transpiration_Shaded, self.Hourly_Soil_Rad, self.Hourly_boundary_layer_resistance_soil, self.Hourly_turbulence_resistance_soil):
    #         hour = 12 - 0.5 * Day_Length + Day_Length * gaussian_weights_x[i]
    #         Hourly_sin_beam = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (hour - 12.) / 24.))

    #         # Diurnal availability of soil water supply
    #         water_supply_Hourly = self.water_supply_for_evaporation * (Hourly_sin_beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
    #         # water_supply_for_evaporation = water_supply_Hourly * self.Top_Layer_Depth / Root_Depth
    #         water_supply_for_evaporation = water_supply_Hourly *self.evaporation_depth/ Root_Depth

    #         # print( Day_Length*3600*(Hourly_sin_beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure)
    #         # print((Hourly_Solar_Radiation/Solar_Radiation)*Day_Length*3600)

    #         # Transpiration from the top soil layer
    #         potential_transpiration_from_evaporative_soil_layer  = (shaded_transpiration + sunlit_transpiration) * self.evaporation_depth / Root_Depth
    #         # print(water_supply_for_evaporation ,potential_evaporation_soil,potential_evaporation_soil / (potential_transpiration_from_evaporative_soil_layer + potential_evaporation_soil) * water_supply_for_evaporation)

    #         # Maximum possible evaporation
    #         Actual_Soil_Evap = min(potential_evaporation_soil, potential_evaporation_soil / (potential_transpiration_from_evaporative_soil_layer + potential_evaporation_soil) * water_supply_for_evaporation)
    #         # water_supply_for_evaporation-=Actual_Soil_Evap
    #         # water_supply_for_evaporation=max(0,water_supply_for_evaporation)
    #         # print(water_supply_for_evaporation ,self.water_supply_for_evaporation)

    #         # print(Actual_Soil_Evap)

    #         Actual_Air_Soil_Temperature_Difference = self.Limit_Function(-25., 25., (Hourly_Soil_Rad - Latent_Heat_of_Vaporization * Actual_Soil_Evap) * (layer_resistance + turbulence_resistance) / Volumetric_Heat_Capacity_Air)
    #         Hourly_Air_Soil_Temperature_Difference.append(Actual_Air_Soil_Temperature_Difference)
    #         Hourly_Actual_Soil_Evap.append(Actual_Soil_Evap)
    #         #print(Actual_Air_Soil_Temperature_Difference)
    #     # print('Next day')
    #     self.Hourly_Actual_Soil_Evap = Hourly_Actual_Soil_Evap
    #     # Aggregation to daily air-soil temperature difference and actual daily evaporation
    #     self.Day_Air_Soil_temp_dif = (Hourly_Air_Soil_Temperature_Difference* gaussian_weights_w).sum()
    #     # Assuming `aggregate_to_daily` is a method that sums or averages Hourly data to a daily total
    #     self.Actual_Daily_Evaporation = Leaf.Leaf.aggregate_to_daily(Hourly_Actual_Soil_Evap, Day_Length)
    #     #print(self.Day_Air_Soil_temp_dif)







    def Calculate_Soil_Temperature(self, tmin, tmax):
        # Calculate the daily average temperature with a weighted average favoring tmax
        daily_avg_temperature = 0.29 * tmin + 0.71 * tmax
        
        # Calculate the nightly average temperature with a weighted average favoring tmin
        nightly_avg_temperature = 0.71 * tmin + 0.29 * tmax
        
        # Calculate the soil's average steady state temperature
        average_soil_temperature = (daily_avg_temperature + self.Day_Air_Soil_temp_dif + nightly_avg_temperature) / 2.
        
        # Calculate the rate of change in soil temperature
        rate_of_change_soil_temp = (average_soil_temperature - self.Soil_Temperature) / self.temperature_change_T
        
        # Print the rate of change in soil temperature for debugging or information purposes
        # print(rate_of_change_soil_temp)
        
        # Update class attributes with the new calculated values
        self.average_soil_temperature = average_soil_temperature
        self.rate_of_change_soil_temp = rate_of_change_soil_temp




    def Calculate_Water_Balance(self,irrigation_amount,Actual_Canopy_Transpiration,Root_Depth):

      
        soil_moisture_array, water_content_array, Root_zone_avg_soil_moisture, Root_zone_avg_water_content,Root_zone_plant_avail_water=self.soil_water_balance(1,self.Soil_Layer_Property, self.planting_depth, self.num_intervals, self.number_of_layers,
                                                               Root_Depth,   irrigation_amount, Actual_Canopy_Transpiration, self.Actual_Daily_Evaporation )

        total_water_soil_profile=0
        total_depth=0
        for j, layer in enumerate(self.Soil_Layer_Property['layers']):
            if j==0:
                layer['Soil_Moisture']=soil_moisture_array[j]
                continue
            else:
                layer['Soil_Moisture']=soil_moisture_array[j]
                total_water_soil_profile+=layer['Soil_Moisture']*layer['Layer_Thickness']
                total_depth+=layer['Layer_Thickness']
                average_soil_moisture_soil_profile=total_water_soil_profile/total_depth

        below_root_zone_water_content=total_water_soil_profile-Root_zone_avg_water_content/1000
        self.average_below_root_zone_moisture=below_root_zone_water_content/(self.total_depth-Root_Depth-self.Evaporative_Depth)
        self.average_root_zone_water_content=Root_zone_avg_water_content #m
        self.water_supply_for_Transpiration=Root_zone_plant_avail_water #m
        self.average_root_zone_moisture=Root_zone_avg_soil_moisture
        self.average_soil_moisture_soil_profile=average_soil_moisture_soil_profile
        # self.average_below_root_zone_moisture=100*self.average_below_root_zone_water_content/(self.total_depth-Root_Depth)


    






   
    def Organic_Carbon_Composition(self,Root_Depth):

        temperature_factor = 47.9 / (1. + np.exp(106. / (self.Soil_Temperature + 18.3)))
        # moisture_factor = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * (self.average_root_zone_water_content + self.average_below_root_zone_water_content) * 100. / self.total_depth / (self.field_capacity_water_content - self.Residual_Soil_Moisture))
        moisture_factor = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * (self.average_soil_moisture_soil_profile -self.Permanent_Wilting_Point_soil_profile ) / (self.Field_Capacity_soil_profile - self.Permanent_Wilting_Point_soil_profile))

        # Extract nitrogen levels
        nitrate_upper_layer, ammonium_upper_layer, nitrate_lower_layer, ammonium_lower_layer =self.nitrate_upper_layer, self.ammonium_upper_layer,self.nitrate_lower_layer,self.ammonium_lower_layer
        
        # Calculate decomposition rates for plant material
        decomposable_pm_rate_change = self.switch_function(nitrate_upper_layer + ammonium_upper_layer + nitrate_lower_layer + ammonium_lower_layer - self.residual_ammonium - self.residual_nitrate, 0., self.dpm_decomposition_rate)
        resistant_pm_rate_change = self.switch_function(nitrate_upper_layer + ammonium_upper_layer + nitrate_lower_layer + ammonium_lower_layer - self.residual_ammonium - self.residual_nitrate, 0., self.rpm_decomposition_rate)
        
        clay_bonus = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * self.clay_percentage))

        carbon_nitrogen_dpm_ratio = (self.initial_decomposable_plant_material + self.resistant_plant_material) / self.Avoid_Zero_Division(self.decomposable_plant_nitrogen + self.resistant_plant_nitrogen)
        decomposed_biomass = self.Microbial_Biomass_content * (1 - math.exp(-temperature_factor * moisture_factor * self.Microbial_Biomass_decomposition_rate / 365)) / self.temperature_change_constant
        decomposed_humus = self.humus_content * (1 - math.exp(-temperature_factor * moisture_factor * self.humification_rate / 365)) / self.temperature_change_constant

        dpm_rate_adjusted = self.switch_function(1.0 / self.Avoid_Zero_Division(carbon_nitrogen_dpm_ratio) - 1.0 / (8.5 * (1.0 + clay_bonus)), decomposable_pm_rate_change, self.dpm_decomposition_rate)
        rpm_rate_adjusted = self.switch_function(1.0 / self.Avoid_Zero_Division(carbon_nitrogen_dpm_ratio) - 1.0 / (8.5 * (1.0 + clay_bonus)), resistant_pm_rate_change, self.rpm_decomposition_rate)
        
        decomposed_dpm = self.initial_decomposable_plant_material * (1.0 - math.exp(-temperature_factor * moisture_factor * dpm_rate_adjusted / 365)) / self.temperature_change_constant
        decomposed_rpm = self.resistant_plant_material * (1.0 - math.exp(-temperature_factor * moisture_factor * rpm_rate_adjusted / 365)) / self.temperature_change_constant
        decomposed_dpn = self.decomposable_plant_nitrogen * (1.0 - math.exp(-temperature_factor * moisture_factor * dpm_rate_adjusted / 365)) / self.temperature_change_constant
        decomposed_rpn = self.resistant_plant_nitrogen * (1.0 - math.exp(-temperature_factor * moisture_factor * rpm_rate_adjusted / 365)) / self.temperature_change_constant

        mineralized_nitrogen = 1.0 / 8.5 * (decomposed_biomass + decomposed_humus) + decomposed_dpn + decomposed_rpn - 1.0 / 8.5 / (1.0 + clay_bonus) * \
                              (decomposed_dpm + decomposed_rpm + decomposed_biomass + decomposed_humus)
        mineralized_nitrogen_upper_layer = (1. - math.exp(-0.065 * Root_Depth)) * mineralized_nitrogen
        mineralized_nitrogen_lower_layer = math.exp(-0.065 * Root_Depth) * mineralized_nitrogen
        
        # Update class attributes with the calculated values
        self.temperature_factor = temperature_factor
        self.clay_bonus = clay_bonus
        self.decomposed_dpm = decomposed_dpm
        self.decomposed_rpm = decomposed_rpm
        self.decomposed_dpn = decomposed_dpn
        self.decomposed_rpn = decomposed_rpn
        self.decomposed_biomass = decomposed_biomass
        self.decomposed_humus = decomposed_humus
        self.mineralized_nitrogen = mineralized_nitrogen
        self.mineralized_nitrogen_upper_layer = mineralized_nitrogen_upper_layer
        self.mineralized_nitrogen_lower_layer = mineralized_nitrogen_lower_layer

    
    

    def Calculate_Nitrogen_Uptake(self, nitrogen_demand, NitrogenFixation_Reserve_Pool_ChangeRate,Root_Depth):
        

        # Calculate supply of ammonium nitrogen from the soil
        ammonium_nitrogen_supply_as = max(0., self.ammonium_upper_layer + (self.mineralized_ammonium_upper_layer - self.nitrate_upper_layer) * self.temperature_change_constant - Root_Depth / self.total_depth * self.residual_ammonium) / self.temperature_change_constant
        
        # Calculate supply of nitrate nitrogen from the soil, factoring in water stress
        nitrate_nitrogen_supply_ns = max(0., self.nitrate_upper_layer + (self.mineralized_nitrate_upper_layer - self.denitrification_upper_layer) * self.temperature_change_constant - Root_Depth / self.total_depth * self.residual_nitrate) / self.temperature_change_constant * self.nitrogen_stress_water_index
        
        # Determine actual nitrogen supply based on soil water index
        ammonium_nitrogen_supply = self.switch_function(self.water_supply_switch, self.ammonium_nitrogen_input_rate, ammonium_nitrogen_supply_as)
        nitrate_nitrogen_supply = self.switch_function(self.water_supply_switch, self.nitrate_nitrogen_input_rate, nitrate_nitrogen_supply_ns)
        
        # Total nitrogen supply
        total_nitrogen_supply = ammonium_nitrogen_supply + nitrate_nitrogen_supply
        
        # Calculate nitrogen uptake, respecting the demand and fixation rate
        ammonium_nitrogen_uptake = min(ammonium_nitrogen_supply, ammonium_nitrogen_supply / max(1e-10, total_nitrogen_supply) * max(0, nitrogen_demand - NitrogenFixation_Reserve_Pool_ChangeRate / self.temperature_change_constant))
        nitrate_nitrogen_uptake = min(nitrate_nitrogen_supply, nitrate_nitrogen_supply / max(1e-10, total_nitrogen_supply) * max(0, nitrogen_demand - NitrogenFixation_Reserve_Pool_ChangeRate / self.temperature_change_constant))
        
        # Total nitrogen uptake, ensuring it does not exceed demand
        nitrogen_uptake = max(0, ammonium_nitrogen_uptake + nitrate_nitrogen_uptake + min(nitrogen_demand, NitrogenFixation_Reserve_Pool_ChangeRate / self.temperature_change_constant))
        
        # Debugging print statements can be commented out or removed in production
        # print(ammonium_nitrogen_uptake, nitrate_nitrogen_uptake, nitrogen_demand, NitrogenFixation_Reserve_Pool_ChangeRate, self.temperature_change_constant)
        
        # Update class attributes with the calculated values
        self.nitrogen_uptake = nitrogen_uptake
        self.nitrate_nitrogen_uptake = nitrate_nitrogen_uptake
        self.ammonium_nitrogen_uptake = ammonium_nitrogen_uptake
        self.total_nitrate_nitrogen_supply = nitrate_nitrogen_supply


    
    
    
    
    def Organic_Nitrogen_Composition(self, rainfall,Root_Depth):

        # Calculate mineralized ammonium in upper and lower layers
        mineralized_ammonium_upper_layer = self.switch_function(self.mineralized_nitrogen, -min((self.ammonium_lower_layer - Root_Depth / self.total_depth * self.residual_ammonium) / self.temperature_change_constant, -self.mineralized_nitrogen_upper_layer), self.mineralized_nitrogen_upper_layer)
        mineralized_ammonium_lower_layer = self.switch_function(self.mineralized_nitrogen, -min((self.ammonium_lower_layer - (self.total_depth - Root_Depth) / self.total_depth * self.residual_ammonium) / self.temperature_change_constant, -self.mineralized_nitrogen_lower_layer), self.mineralized_nitrogen_lower_layer)

        # Calculate mineralized nitrate in upper and lower layers
        mineralized_nitrate_upper_layer = self.switch_function(self.mineralized_nitrogen, -min(self.nitrate_upper_layer/self.temperature_change_constant, -self.mineralized_nitrogen_upper_layer + mineralized_ammonium_upper_layer), 0.)
        mineralized_nitrate_lower_layer = self.switch_function(self.mineralized_nitrogen, -min(self.nitrate_lower_layer/self.temperature_change_constant, -self.mineralized_nitrogen_lower_layer + mineralized_ammonium_lower_layer), 0.)

        # Calculate moisture factors for upper and lower layers
        moisture_factor_upper_layer = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * (self.average_root_zone_moisture - self.Permanent_Wilting_Point_soil_profile) / (self.Field_Capacity_soil_profile - self.Permanent_Wilting_Point_soil_profile))
        moisture_factor_lower_layer = self.Limit_Function(0.2, 1.0, 0.2 + 0.8 * (self.average_below_root_zone_moisture - self.Permanent_Wilting_Point_soil_profile)/  (self.Field_Capacity_soil_profile - self.Permanent_Wilting_Point_soil_profile))

                                                          
                                                          

        # Calculate nitrification in upper and lower layers
        nitrification_upper_layer = max(0., (self.ammonium_upper_layer + mineralized_ammonium_upper_layer * self.temperature_change_constant - Root_Depth / self.total_depth * self.residual_ammonium)) * (1 - np.exp(-self.temperature_factor * moisture_factor_upper_layer * 0.6 / 7)) / self.temperature_change_constant
        nitrification_lower_layer = max(0., (self.ammonium_lower_layer + mineralized_ammonium_lower_layer * self.temperature_change_constant - (self.total_depth - Root_Depth) / self.total_depth * self.residual_ammonium)) * (1 - np.exp(-self.temperature_factor * moisture_factor_lower_layer * 0.6 / 7)) / self.temperature_change_constant

        # Calculate CO2 respiration
        respiration_co2 = self.clay_bonus / (1.0 + self.clay_bonus) * (self.decomposed_dpm + self.decomposed_rpm + self.decomposed_biomass + self.decomposed_humus)

        # Calculate denitrification in upper and lower layers
        denitrification_upper_layer = .0005 * max(0., self.nitrate_upper_layer + mineralized_nitrate_upper_layer * self.temperature_change_constant - Root_Depth / self.total_depth * self.residual_nitrate) * respiration_co2 * (1. - np.exp(-0.065 * Root_Depth))
        denitrification_lower_layer = .0005 * max(0., self.nitrate_lower_layer + mineralized_nitrate_lower_layer * self.temperature_change_constant - (self.total_depth - Root_Depth) / self.total_depth * self.residual_nitrate) * respiration_co2 * np.exp(-0.065 * Root_Depth)

        # Calculate water stress factor
        nitrogen_stress_water_index = min(1,  (self.average_root_zone_moisture - self.Permanent_Wilting_Point_soil_profile) / (self.Field_Capacity_soil_profile - self.Permanent_Wilting_Point_soil_profile))

        # Calculate total nitrogen and ammonia volatilization
        total_ammonium = self.ammonium_upper_layer + self.ammonium_lower_layer
        total_nitrate = self.nitrate_upper_layer + self.nitrate_lower_layer
        total_mineral_nitrogen = total_ammonium + total_nitrate
        ammonia_volatilization = self.switch_function(rainfall - 1., 0.15, 0.) * self.N_volatilization

        # Update class attributes or return calculated values as needed
        self.nitrification_upper_layer = nitrification_upper_layer
        self.nitrification_lower_layer = nitrification_lower_layer
        self.mineralized_ammonium_upper_layer = mineralized_ammonium_upper_layer
        self.mineralized_ammonium_lower_layer = mineralized_ammonium_lower_layer
        self.mineralized_nitrate_upper_layer = mineralized_nitrate_upper_layer
        self.mineralized_nitrate_lower_layer = mineralized_nitrate_lower_layer
        self.denitrification_upper_layer = denitrification_upper_layer
        self.denitrification_lower_layer = denitrification_lower_layer
        self.nitrogen_stress_water_index = nitrogen_stress_water_index
        self.ammonia_volatilization = ammonia_volatilization




    def Calculate_Soil_Organic_Nitrogen_ChangeRate(self, days_after_planting, root_depth_growth_rate,Root_Depth):

        ammonium_fertilization_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]
        nitrate_fertilization_schedule = [(5.0, 0), (15.0, 0), (25.0, 0), (35.0, 0), (45.0, 0), (55.0, 0), (65.0, 0), (75.0, 0)]

        # Calculate ammonium and nitrate fertilization amounts based on schedule and days after planting
        ammonium_fertilizer_application = sum([self.Function_Adjust_Switch(time - days_after_planting, 0., amount, 0.) for time, amount in ammonium_fertilization_schedule])
        nitrate_fertilizer_application = sum([self.Function_Adjust_Switch(time - days_after_planting, 0., amount, 0.) for time, amount in nitrate_fertilization_schedule])

        # Residual soil fertilizer after application (assuming a certain percentage is immediately available)
        N_volatilization_change_rate = ammonium_fertilizer_application - self.N_volatilization / 3
        
        # Leaching losses for ammonium and nitrate in lower soil layer
        leaching_nitrate_lower_layer = max(0, self.nitrate_lower_layer + (self.mineralized_nitrate_lower_layer - self.denitrification_upper_layer) * self.temperature_change_constant - (self.total_depth - Root_Depth) / self.total_depth * self.groundwater_recharge) * min(self.groundwater_recharge / self.Saturated_Soil_Moisture_soil_profile / (self.total_depth - Root_Depth) / 10, 1)
        leaching_nitrate_upper_layer = max(0, (self.total_nitrate_nitrogen_supply - self.nitrate_nitrogen_uptake) * self.temperature_change_constant - Root_Depth / self.total_depth * self.groundwater_recharge) * min((self.rain_and_irrigation - self.recharge_to_upper_layer) / self.Saturated_Soil_Moisture_soil_profile / Root_Depth / 10, 1)
        
        # Nitrogen layer adjustments due to root depth ratio
        layer_adjustment_ammonium = root_depth_growth_rate / (self.total_depth - Root_Depth) * self.ammonium_lower_layer
        layer_adjustment_nitrate = root_depth_growth_rate / (self.total_depth - Root_Depth) * self.nitrate_lower_layer

        # Update class attributes with the calculated values
        self.N_volatilization_change_rate = N_volatilization_change_rate
        self.leaching_nitrate_lower_layer = leaching_nitrate_lower_layer
        self.leaching_nitrate_upper_layer = leaching_nitrate_upper_layer
        self.layer_adjustment_ammonium = layer_adjustment_ammonium
        self.layer_adjustment_nitrate = layer_adjustment_nitrate
        self.ammonium_fertilizer_application = ammonium_fertilizer_application
        self.nitrate_fertilizer_application = nitrate_fertilizer_application


    def Calculate_Soil_Organic_Carbon_ChangeRate(self, litter_carbon, litter_nitrogen):
        # Decomposition and transformation rates for soil organic carbon
        decomposable_pm_change_rate = litter_carbon * self.plant_material_dpm_rpm_ratio / (1. + self.plant_material_dpm_rpm_ratio) - self.decomposed_dpm
        resistant_pm_change_rate = litter_carbon * 1. / (1. + self.plant_material_dpm_rpm_ratio) - self.decomposed_rpm
        decomposable_pn_change_rate = litter_nitrogen / (1. + 40. / self.plant_material_dpm_rpm_ratio / 100.) - self.decomposed_dpn
        resistant_pn_change_rate = litter_nitrogen / (1. + 100. * self.plant_material_dpm_rpm_ratio / 40.) - self.decomposed_rpn
        
        # Microbial_Biomass and humus change rates based on decomposition
        Microbial_Biomass_change_rate = 0.46 / (1.0 + self.clay_bonus) * (self.decomposed_dpm + self.decomposed_rpm + self.decomposed_biomass + self.decomposed_humus) - self.decomposed_biomass
        humus_change_rate = 0.54 / (1.0 + self.clay_bonus) * (self.decomposed_dpm + self.decomposed_rpm + self.decomposed_biomass + self.decomposed_humus) - self.decomposed_humus
        
        # Ammonium and nitrate nitrogen in the upper and lower layers
        ammonium_upper_layer_change = self.ammonium_fertilizer_application + self.mineralized_ammonium_upper_layer + self.layer_adjustment_ammonium - self.switch_function(self.water_supply_switch, 0.0, self.ammonium_nitrogen_uptake) - self.nitrification_upper_layer - self.ammonia_volatilization
        ammonium_lower_layer_change = self.mineralized_ammonium_lower_layer - self.layer_adjustment_ammonium - self.nitrification_lower_layer
        nitrate_lower_layer_change = self.leaching_nitrate_upper_layer + self.mineralized_nitrate_lower_layer + self.nitrification_lower_layer - self.layer_adjustment_nitrate - self.denitrification_lower_layer - self.leaching_nitrate_lower_layer
        nitrate_upper_layer_change = self.nitrate_fertilizer_application + self.mineralized_nitrate_upper_layer + self.nitrification_upper_layer + self.layer_adjustment_nitrate - self.switch_function(self.water_supply_switch, 0.0, self.nitrate_nitrogen_uptake) - self.denitrification_upper_layer - self.leaching_nitrate_upper_layer
        
        # Update class attributes with the calculated values
        self.decomposable_pm_change_rate = decomposable_pm_change_rate
        self.resistant_pm_change_rate = resistant_pm_change_rate
        self.decomposable_pn_change_rate = decomposable_pn_change_rate
        self.resistant_pn_change_rate = resistant_pn_change_rate
        self.Microbial_Biomass_change_rate = Microbial_Biomass_change_rate
        self.humus_change_rate = humus_change_rate
        self.ammonium_upper_layer_change = ammonium_upper_layer_change
        self.ammonium_lower_layer_change = ammonium_lower_layer_change
        self.nitrate_lower_layer_change = nitrate_lower_layer_change
        self.nitrate_upper_layer_change = nitrate_upper_layer_change

        
        
        
        
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
        self.Soil_Temperature += self.rate_of_change_soil_temp
    
        # Update water uptake and loss limits with respective rate changes
        # self.water_content_upper_layer += self.water_content_upper_layer_ChangeRate 
        # self.water_content_lower_layer += self.water_content_lower_layer_ChangeRate 
        
        # Update organic matter and nitrogen pools with their respective rate changes
        self.decomposable_plant_material += self.decomposable_pm_change_rate
        self.resistant_plant_material += self.resistant_pm_change_rate
        self.humus_content += self.humus_change_rate
        self.decomposable_plant_nitrogen += self.decomposable_pn_change_rate
        self.resistant_plant_nitrogen += self.resistant_pn_change_rate
        
        self.Microbial_Biomass_content += self.Microbial_Biomass_change_rate

        # Update total nitrogen leached and nitrogen content in upper and lower soil layers
        self.total_nitrogen_leached += self.leaching_nitrate_lower_layer
        
        self.ammonium_upper_layer += self.ammonium_upper_layer_change
        self.ammonium_lower_layer += self.ammonium_lower_layer_change
        self.nitrate_upper_layer += self.nitrate_upper_layer_change
        self.nitrate_lower_layer += self.nitrate_lower_layer_change
        
        # Update soil ammonium volatilization rate
        self.N_volatilization += self.N_volatilization_change_rate

     