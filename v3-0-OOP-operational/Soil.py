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

import json
import math
from collections import namedtuple
# Soil water characteristics.
SWC = namedtuple('SWC', 'sat fc pwp psd porosity airentry')

# Soil texture.
Texture = namedtuple('Texture', ['clay', 'sand', 'om'])

# Soil water characteristics in the rooting zone.
RootZone = namedtuple('RootZone', 'wc vwc Root_Zone_Plant_Avail sat fc pwp Evaporative_layer_water_content')

# Actual evapotranspiration (ET).
ActualET = namedtuple('ActualET', 'crop soil')

# Water fluxes into a given soil layer.
Fluxes = namedtuple('Fluxes', 't e influx outflux netflux')

class Soil:
        
    __accdepth = 0.0    # internal use: used to determine a layer's depth
    
    def __init__(self, Soil_Layer_Property):
        self.num_intervals = Soil_Layer_Property['num_intervals']
        self.planting_depth = Soil_Layer_Property['planting_depth']
        self.rootdepth = Soil_Layer_Property['rootdepth']
        self.has_watertable = Soil_Layer_Property['has_watertable']
        self.number_of_layers = Soil_Layer_Property['number_of_layers']
        self.layers = []
   
        layers = Soil_Layer_Property['layers']
        self.alpha= Soil_Layer_Property['alpha']
        self.layer_depths = [layer['thick'] * 100 for layer in Soil_Layer_Property['layers']]  # Convert thickness to cm

        
        
        for i in range(self.number_of_layers):
            layer = self._create_layer(layers[i])
            self.layers.append(layer)
   
        for i in range(self.number_of_layers):
            prevlayer = self.layers[i - 1] if i > 0 else None
            nextlayer = self.layers[i + 1] if i < self.number_of_layers - 1 else None
            self._initialize_layer(self.layers[i], prevlayer, nextlayer)
   
        self.__pf = [{field: 0.0 for field in Fluxes._fields} for _ in range(self.number_of_layers)]
        self.__prz = {field: 0.0 for field in RootZone._fields}
        self._rootzone_water()
        self.rootwater = RootZone(*[0.0] * len(RootZone._fields))
        self.netrain = 0.0
        # self.Evaporative_Depth=Soil_Layer_Property['layers'][0]['thick']
        self.Evaporative_Depth=Soil_Layer_Property['Evaporative_Depth']
        self.clay_percentage = 12.0  # Clay content in soil
        self.soil_moisure_evaporative_layer = layers[0]['vwc']

        root_depth_ini = Soil_Layer_Property['rootdepth']*100 #cm


        
        total_water_FC_profile=0
        total_water_S_profile=0
        total_water_PWP_profile=0
        total_depth=0
        for j, layer in enumerate(Soil_Layer_Property['layers']):

            total_depth+=layer['thick']

            total_water_FC_profile+=layer['ThetaFC']*layer['thick']
            self.Field_Capacity_soil_profile=total_water_FC_profile/total_depth

            total_water_S_profile+=layer['ThetaS']*layer['thick']
            self.Saturated_Soil_Moisture_soil_profile=total_water_S_profile/total_depth
            
            total_water_PWP_profile+=layer['ThetaPWP']*layer['thick']
            self.Permanent_Wilting_Point_soil_profile=total_water_PWP_profile/total_depth
        self.total_depth=total_depth



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

        self.Soil_Temperature =  self.initial_soil_temp 
        self.ammonium_upper_layer = ammonium_upper_layer_initial
        self.nitrate_upper_layer = nitrate_upper_layer_initial
        self.nitrate_lower_layer = nitrate_lower_layer_initial
        self.ammonium_lower_layer = ammonium_lower_layer_initial
        self.total_nitrogen_leached = 0
        self.N_volatilization = 0

        self.potential_SoilRad_daily = 0
        self.Day_Air_Soil_temp_dif  = 0
        self.Actual_Daily_Evaporation =0
        
        # self.hourly_Soil_Evap=[]
        self.hourly_Soil_Rad=[]
        self.Hourly_Actual_Soil_Evap=[]
        self.potential_evaporation=0
        self.Actual_Daily_Transpiration=0
        
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



    def _create_layer(self, layer_data):
        layer = {
            'thick': layer_data['thick'],
            'texture': Texture(layer_data['texture']['clay'], layer_data['texture']['sand'], layer_data['texture']['om']),
            'vwc': layer_data['vwc'],
            'wc': layer_data['vwc'] * layer_data['thick'] * 1000,
            'ThetaS': layer_data['ThetaS'],
            'ThetaFC': layer_data['ThetaFC'],
            'ThetaPWP': layer_data['ThetaPWP'],
            'ksat': layer_data['KSAT'],
            'Pore_Size_Distribution': layer_data['Pore_Size_Distribution'],
            'PORS': layer_data['PORS'],
            'swc': SWC(layer_data['ThetaS'], layer_data['ThetaFC'], layer_data['ThetaPWP'], 
                       layer_data['Pore_Size_Distribution'], layer_data['PORS'], 5.0),
            'accthick': 0.0,
            'depth': 0.0,
            'matric': 0.0,
            'gravity': 0.0,
            'fluxes': Fluxes(0.0, 0.0, 0.0, 0.0, 0.0),
            'prev': None,
            'next': None
        }
        return layer

    def _initialize_layer(self, layer, prevlayer, nextlayer):
        layer['prev'] = prevlayer
        layer['next'] = nextlayer

        prevaccthick = prevlayer['accthick'] if prevlayer else 0.0
        layer['accthick'] = layer['thick'] + prevaccthick
        prevthick = prevlayer['thick'] if prevlayer else 0.0
        d = 0.5 * (prevthick + layer['thick'])
        layer['depth'] = Soil.__accdepth + d
        Soil.__accdepth += d


        c, s, om = layer['texture']
        s /= 100
        c /= 100
        
        
        awc = layer['ThetaS'] - layer['ThetaFC']
        n1 = -21.674 * s - 27.932 * c - 81.975 * awc + 71.121 * (s * awc)
        n2 = 8.294 * (c * awc) + 14.05 * (s * c) + 27.161
        aet = n1 + n2
        ae = max(0.0, aet + (0.02 * aet ** 2 - 0.113 * aet - 0.7))
    

        layer['swc'] = SWC(layer['ThetaS'], layer['ThetaFC'], layer['ThetaPWP'], 
                           layer['Pore_Size_Distribution'], layer['PORS'], ae)

        
        if layer['vwc'] < 0:
            vwc = -layer['vwc']
            fc = layer['swc'].fc
            if 1 <= vwc <= 2:
                sat = layer['swc'].sat
                vwc = sat - (vwc - 1) * (sat - fc)
            elif 2 < vwc <= 3:
                pwp = layer['swc'].pwp
                vwc = fc - (vwc - 2) * (fc - pwp)
            else:
                vwc = fc
            layer['vwc'] = vwc
            layer['wc'] = layer['vwc'] * layer['thick'] * 1000

        self._update_heads_k(layer)

    def _update_heads_k(self, layer):
        fc = layer['swc'].fc
        vwc = layer['vwc']
        if vwc >= fc:
            df = vwc - fc
            hm = 33 - (33 - layer['swc'].airentry) * df / (layer['swc'].sat - fc)
            hm /= 10
        else:
            b = 1 / layer['swc'].psd
            a = math.exp(3.496508 + b * math.log(fc))
            hm = (a * max(0.05, vwc) ** (-b)) / 10
        layer['matric'] = max(0.0, hm)
        layer['gravity'] = layer['depth']
        ae = layer['swc'].airentry / 10
        hm = layer['matric']
        ratio = layer['vwc'] / layer['swc'].sat
        if hm > ae:
            layer['k'] = layer['ksat'] * ratio ** (3 + 2 / layer['swc'].psd)
        else:
            layer['k'] = layer['ksat']



    def _rootzone_water(self):
        Root_zone_water_content = 0.0
        total_thickness = 0.0
        Root_zone_plant_avail_water = 0.0
        Evaporative_layer_water_content=0        
        for layer_indx, layer in enumerate(self.layers):
            if layer['accthick'] <0.3:
                Evaporative_layer_water_content+= (layer['vwc'] - 0.005) * layer['thick']*1000
            # if    evaporative_depth_remaining > 0 :
            #       Evaporative_layer_water_content+=(layer['vwc'] - 0.005) *min(evaporative_depth_remaining, layer['thick'])*1000 
            #       evaporative_depth_remaining-=min(evaporative_depth_remaining, layer['thick'])
            top_of_layer = layer['accthick'] - layer['thick']
            bottom_of_layer = layer['accthick']
            if top_of_layer >= self.planting_depth and bottom_of_layer <= self.planting_depth + self.rootdepth:
                thickness_within_zone = layer['thick']
            elif top_of_layer < self.planting_depth and bottom_of_layer > self.planting_depth:
                thickness_within_zone = bottom_of_layer - self.planting_depth
            elif top_of_layer < self.planting_depth + self.rootdepth and bottom_of_layer > self.planting_depth + self.rootdepth:
                thickness_within_zone = self.planting_depth + self.rootdepth - top_of_layer
            else:
                continue

            Root_zone_water_content += (layer['vwc'] - 0.005) * thickness_within_zone
            Root_zone_plant_avail_water += max(0, (layer['vwc'] - layer['swc'].pwp) * thickness_within_zone)
            total_thickness += thickness_within_zone

        if total_thickness > 0:
            avg_soil_moisture = Root_zone_water_content / total_thickness
            avg_water_content = Root_zone_water_content * 1000
            Root_zone_plant_avail_water = Root_zone_plant_avail_water * 1000
        else:
            avg_soil_moisture = 0.0
            avg_water_content = 0.0
            Root_zone_plant_avail_water = 0.0

        self.__prz['wc'] = avg_water_content
        self.__prz['vwc'] = avg_soil_moisture
        self.__prz['Root_Zone_Plant_Avail'] = Root_zone_plant_avail_water
        self.__prz['Evaporative_layer_water_content'] = Evaporative_layer_water_content

    # Function to calculate evaporation as a function of depth
    def evaporation_profile(self,petsoil, layer_depths, alpha):
        depths = np.cumsum(layer_depths)  # Cumulative depth at the bottom of each layer
        evaporation_rates = petsoil * np.exp(-alpha * depths)  # Evaporation at each layer
        return evaporation_rates


    def _calc_water_fluxes(self, cummfluxes, petcrop, petsoil):
        prvpsi = 0.0
        
        evaporation_per_layer = self.evaporation_profile(petsoil, self.layer_depths, self.alpha)
        # print('evaporation_per_layer',evaporation_per_layer)
        for idx in range(self.number_of_layers):
            cur = self.layers[idx]
            prv = cur['prev']

            self._update_heads_k(cur)

            if prv is None:
                ei = evaporation_per_layer[idx] / 1000
                ti = 0
                prvpsi = 0.0
            else:
                ei = evaporation_per_layer[idx]/ 1000
                cj = min(1.0, (cur['accthick'] - self.planting_depth) / self.rootdepth)
                curpsi = 1.8 * cj - 0.8 * cj ** 2
                ti = petcrop * (curpsi - prvpsi) / 1000
                prvpsi = curpsi

            if prv is not None:
                n = math.log(cur['k']) - math.log(prv['k'])
                k = (cur['k'] - prv['k']) / n if n != 0.0 else cur['k']
                grad = (cur['matric'] + cur['gravity']) - (prv['matric'] + prv['gravity'])
                grad /= (cur['depth'] - prv['depth'])
                grad = max(0, grad)
                curinflux = k * grad
            else:
                netrain = self.netrain / 1000
                curinflux = min(netrain, cur['ksat'])

            self.__pf[idx]['t'] = ti
            self.__pf[idx]['e'] = ei
            self.__pf[idx]['influx'] = curinflux
        total_ei=0
        total_ti=0
        for idx in range(self.number_of_layers):
            cur = self.layers[idx]
            nxt = cur['next']

            wc = cur['vwc'] * cur['thick']
            wc_e = (cur['vwc'] - 0.005) * cur['thick']  # available water for evaporation
            wc_ti = (cur['vwc'] - cur['swc'].pwp) * cur['thick']  # available water for transpiration

            influx = self.__pf[idx]['influx']

            if nxt is not None:
                outflux = self.__pf[idx + 1]['influx']
            elif not self.has_watertable:
                outflux = cur['k']
            else:
                outflux = self._influx_from_watertable()

            ti = self.__pf[idx]['t']
            ei = self.__pf[idx]['e']
            

            if idx == 0:
                if wc_e + influx - ei <= 0:
                    # print(f'E is not sufficient in layer {idx}')
                    ei = wc_e + influx
                    outflux = 0
                else:
                    
                    if (wc_e + influx - ei)>0 and (wc_e + influx - ei <= outflux):
                        outflux = wc_e + influx - ei
                if nxt is not None:
                    self.__pf[idx + 1]['influx'] = outflux

                netflux = influx - ei - outflux
                
                
                
            elif cur['accthick'] <0.3:
                if wc_e + influx - ei <= 0:
                    # print(f'E is not sufficient in layer {idx}')
                    ei = wc_e + influx
                    ti=0
                    outflux = 0
                
                else:
                    if wc_ti + influx - ti -ei <= 0:
                        # print(f'T is not sufficient in layer {idx}')
                        ti = wc_ti + influx- ei
                        outflux = 0
                    if (wc_ti + influx - ti - ei>0) & (wc_ti + influx - ti - ei<= outflux):
                        outflux = wc_ti + influx - ti - ei
                if nxt is not None:
                    self.__pf[idx + 1]['influx'] = outflux

                netflux = influx - ti - ei- outflux

            else :
                if wc_e + influx - ti <= 0:
                    # print(f'E is not sufficient in layer {idx}')
                    ti = wc_e + influx
                    outflux = 0
                
        
                if (wc_ti + influx - ti >0) & (wc_ti + influx - ti <= outflux):
                    outflux = wc_ti + influx - ti - ei
                if nxt is not None:
                    self.__pf[idx + 1]['influx'] = outflux

                netflux = influx - ti -  outflux
            total_ei+=ei
            total_ti+=ti
            self.__pf[idx]['outflux'] = outflux
            self.__pf[idx]['netflux'] = netflux
            self.__pf[idx]['influx'] = influx

            wc += self.__pf[idx]['netflux'] / self.num_intervals
        
            
            cur['vwc'] = max(0.005, min(cur['swc'].sat, wc / cur['thick']))
            
    
            
            cur['wc'] = cur['vwc'] * cur['thick'] * 1000

            for field in Fluxes._fields:
                n1 = self.__pf[idx][field] / self.num_intervals
                cummfluxes[idx][field] += n1

        self._rootzone_water()
        self.Actual_Daily_Evaporation=total_ei*1000
        self.Actual_Daily_Transpiration=total_ti*1000
        if self.Actual_Daily_Transpiration <= 0:
            self.Actual_Daily_Transpiration = 1e-10
        # print('Soil 535','A_evap',total_ei,'A_tranpiration',total_ti, 'PT',petcrop)

    def daily_water_balance(self, rain, petcrop, petsoil,rootdepth):
        self.netrain = rain
        self.rootdepth = rootdepth

        cummfluxes = [{field: 0.0 for field in Fluxes._fields} for _ in range(self.number_of_layers)]

        for i in range(self.num_intervals):
            self._calc_water_fluxes(cummfluxes, petcrop, petsoil)

        self.rootwater = RootZone(*[self.__prz[field] for field in RootZone._fields])
        for i, layer in enumerate(self.layers):
            layer['fluxes'] = Fluxes(*[cummfluxes[i][field] for field in Fluxes._fields])

    def _influx_from_watertable(self):
        last = self.layers[-1]
        k = (last['ksat'] - last['k']) / (math.log(last['ksat']) - math.log(last['k']))
        hm = (33 - (33 - last['swc'].airentry)) / 10
        hg = last['accthick']
        tothead = hm + hg
        return k * (tothead - (last['matric'] + last['gravity'])) / (last['thick'] * 0.5)

    def run_soil_water_model(self, duration, rain,  transpiration, evaporation,rootdepth):
        nlayers = self.number_of_layers
        res = {
                    'wc': [[] for _ in range(nlayers)],
                    'vwc': [[] for _ in range(nlayers)]
                }
        for i in range(duration):
            # day = i + 1
            self.daily_water_balance(rain,  transpiration, evaporation,rootdepth)

            plantwc = self.rootwater.Root_Zone_Plant_Avail
            evapowc = self.rootwater.Evaporative_layer_water_content
            Root_zone_avg_water_content = self.rootwater.wc
            Root_zone_avg_soil_moisture = self.rootwater.vwc
        wc_result=[]    
        vwc_result=[]
        for j, layer in enumerate(self.layers):
            wc_result.append(layer['wc'])
            vwc_result.append(layer['vwc'])

        self.results = res  # Store the results in an attribute
        self.soil_moisure_evaporative_layer=vwc_result[0]
        return plantwc, evapowc,Root_zone_avg_water_content,Root_zone_avg_soil_moisture, wc_result, vwc_result
        

    # def Current_Soil_Water_Status(self, irrigation_amount,Root_Depth):
    #     plantwc, evapowc, _, _,_,soil_moisture_array=self.run_soil_water_model(1,irrigation_amount,  0, 0, Root_Depth,  )
    #     self.water_supply_for_evaporation=evapowc
    #     self.soil_moisure_evaporative_layer=soil_moisture_array[0]
    #     self.water_supply_for_Transpiration=plantwc

    def Calculate_Soil_Evaporation(self,irrigation_amount, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination,
                                             Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature,
                                             Min_Temperature, Vapour_Pressure_Deficit, Wind_Speed, 
                                             soil_resistance_to_evaporation,Root_Depth,Leaf_Blade_Angle, Total_Leaf_Area_Index,
                                             Light_Extinction_Coefficient, Potential_Canopy_Transpiration):

            
            # Convert daily climate data to Hourly data
            Hourly_climate_data, gaussian_weights = Leaf.Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temperature, Min_Temperature, Wind_Speed)
            potential_evaporation_soil_first_estimate=[]
            for gaussian_weight,(Hourly_Solar_Constant, Hourly_temperature, Hourly_sin_beam, Hourly_Solar_Radiation, Hourly_wind_speed) in zip(gaussian_weights,Hourly_climate_data):
                
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
                potential_evaporation_soil = max(0, potential_evaporation_soil)
                potential_evaporation_soil_first_estimate.append(potential_evaporation_soil)
            potential_evaporation_soil_first_estimate_daily = Leaf.Leaf.aggregate_to_daily(potential_evaporation_soil_first_estimate, Day_Length)
            self.potential_evaporation=potential_evaporation_soil_first_estimate_daily
            # Maximum_possible_Evaporation=self.water_supply_for_evaporation

            # if potential_evaporation_soil_first_estimate_daily > Maximum_possible_Evaporation :
            #     Water_Stress_Fraction=Maximum_possible_Evaporation/potential_evaporation_soil_first_estimate_daily
            #     potential_evaporation_soil_estimate=np.array(potential_evaporation_soil_first_estimate)*Water_Stress_Fraction
            # else:
            #     potential_evaporation_soil_estimate=np.array(potential_evaporation_soil_first_estimate)

            # Hourly_Actual_Soil_Evap = []
            # Hourly_Air_Soil_temp_dif = []
            # Hourly_Soil_Rad = []
            # Hourly_boundary_layer_resistance_soil=[]
            # Hourly_turbulence_resistance_soil=[]

            # for potential_evaporation_soil_estimate_instant,(Hourly_Solar_Constant, Hourly_temperature, Hourly_sin_beam, Hourly_Solar_Radiation, Hourly_wind_speed), shaded_transpiration, sunlit_transpiration in zip(potential_evaporation_soil_estimate,Hourly_climate_data, Hourly_transpiration_Shaded, Hourly_transpiration_Sunlit):
                
            #     Incoming_PAR = 0.5 * Hourly_Solar_Radiation  # Partition of incoming solar radiation to PAR
            #     Incoming_NIR = 0.5 * Hourly_Solar_Radiation  # Partition of incoming solar radiation to NIR
            #     # Atmospheric transmissivity based on incoming PAR and solar geometry
            #     atmospheric_transmissivity = Incoming_PAR / (0.5 * Hourly_Solar_Constant * Hourly_sin_beam)
           
            #     # Saturation vapor pressure and vapor pressure deficit calculations
            #     saturation_vapor_pressure = 0.611 * math.exp(17.4 * Hourly_temperature / (Hourly_temperature + 239.))
            #     vapor_pressure_deficit = max(0., saturation_vapor_pressure - Vapour_Pressure_Deficit)
                
                
            #     slope_vapor_pressure_temperature = 4158.6 * saturation_vapor_pressure / (Hourly_temperature + 239.) ** 2
           
            #     # Turbulence resistance for soil using logarithmic wind profile
            #     turbulence_resistance_soil = 0.74 * (np.log(56.)) ** 2 / (0.4 ** 2 * Hourly_wind_speed)
                
            #     # Boundary layer resistance for soil considering wind speed and leaf area
            #     boundary_layer_resistance_soil = 172. * np.sqrt(0.05 / max(0.1, Hourly_wind_speed * np.exp(-Light_Extinction_Coefficient * Total_Leaf_Area_Index)))
                
            
            
            
            #     water_boundary_layer_resistance_soil = 0.93 * boundary_layer_resistance_soil
                
            #     # Calculation of the fractional diffuse light (frdf) based on atmospheric transmissivity
            #     if atmospheric_transmissivity < 0.22:
            #         fractional_diffuse_light = 1
            #     elif 0.22 < atmospheric_transmissivity <= 0.35:
            #         fractional_diffuse_light = 1 - 6.4 * (atmospheric_transmissivity - 0.22) ** 2
            #     else:
            #         fractional_diffuse_light = 1.47 - 1.66 * atmospheric_transmissivity
           
            #     # Ensuring a minimum threshold for the fractional diffuse light
            #     fractional_diffuse_light = max(fractional_diffuse_light, 0.15 + 0.85 * (1 - np.exp(-0.1 / Hourly_sin_beam)))
           
            #     # Incoming diffuse PAR (PARDF) and direct PAR (PARDR) based on the calculated fractional diffuse light
            #     PAR_Diffuse = Incoming_PAR * fractional_diffuse_light
            #     PAR_Direct = Incoming_PAR - PAR_Diffuse
           
            #     # Extinction and RefLection coefficients
            #     # Convert leaf blade angle from degrees to radians for calculations
            #     Leaf_Blade_Angle_Radians = Leaf_Blade_Angle * np.pi / 180
                
            #     # Calculate the direct beam extinction coefficient for PAR using leaf blade angle
            #     Direct_Beam_Extinction_Coefficient_PAR = Leaf.Leaf.KDR_Coeff(Hourly_sin_beam, Leaf_Blade_Angle_Radians)
                
            #     # Leaf scattering coefficients for PAR and NIR
            #     Scattering_Coefficient_PAR = 0.2  # Scattering coefficient for Photosynthetically Active Radiation
            #     Scattering_Coefficient_NIR = 0.8  # Scattering coefficient for Near-Infrared Radiation
                
            #     # Calculate diffuse extinction coefficients for PAR and NIR using total leaf area index and leaf blade angle
            #     Diffuse_Extinction_Coefficient_PAR = Leaf.Leaf.KDF_Coeff(Total_Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            #     Diffuse_Extinction_Coefficient_NIR = Leaf.Leaf.KDF_Coeff(Total_Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
                
            #     # Calculate beam and canopy reflection coefficients for PAR and NIR
            #     Beam_Reflection_Coefficient_PAR, Canopy_Reflection_Coefficient_PAR = Leaf.Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient_PAR)
            #     Beam_Reflection_Coefficient_NIR, Canopy_Reflection_Coefficient_NIR = Leaf.Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient_PAR)
           
                
            #     # Incoming diffuse NIR (NIRDF) and direct NIR (NIRDR)
            #     NIR_Diffuse = Incoming_NIR * fractional_diffuse_light
            #     NIR_Direct = Incoming_NIR - NIR_Diffuse
                
            #     # Absorbed total radiation (PAR + NIR) by soil
            #     soil_par_reflection = 0.1  # Soil PAR reflection coefficient
            #     # Soil NIR reflection coefficient, varying with soil water content (wcul)
            #     soil_nir_reflection = self.switch_function(self.soil_moisure_evaporative_layer - 0.5, 0.52 - 0.68 * self.soil_moisure_evaporative_layer, 0.18)
            #     absorbed_total_radiation_by_soil = (1 - soil_par_reflection) * (PAR_Direct * np.exp(-Beam_Reflection_Coefficient_PAR * Total_Leaf_Area_Index) + PAR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_PAR * Total_Leaf_Area_Index)) + \
            #                                        (1 - soil_nir_reflection) * (NIR_Direct * np.exp(-Beam_Reflection_Coefficient_NIR * Total_Leaf_Area_Index) + NIR_Diffuse * np.exp(-Diffuse_Extinction_Coefficient_NIR * Total_Leaf_Area_Index))
               
           
                
            #     # Temperature difference driven by soil evaporation and net radiation
            #     temperature_difference_soil = self.Limit_Function(-25., 25., (soil_absorbed_radiation - Latent_Heat_of_Water_Vaporization * potential_evaporation_soil_estimate_instant) * (boundary_layer_resistance_soil + turbulence_resistance_soil) / Volumetric_Heat_Capacity_Air) 
            #     average_soil_temperature = Hourly_temperature + temperature_difference_soil
                
            #     # Recalculate potential evaporation with updated soil temperature
            #     saturation_vapor_pressure_soil = 0.611 * math.exp(17.4 * average_soil_temperature / (average_soil_temperature + 239.))
            #     slope_vapor_pressure_curve_soil = (saturation_vapor_pressure_soil - saturation_vapor_pressure) / self.Avoid_Zero_Division(temperature_difference_soil)
            #     potential_evaporation_second_estimate, soil_absorbed_radiation_second_estimate = Leaf.Leaf.Penman_Monteith(soil_resistance_to_evaporation, turbulence_resistance_soil, water_boundary_layer_resistance_soil, boundary_layer_resistance_soil, absorbed_total_radiation_by_soil, atmospheric_transmissivity, 1, average_soil_temperature, Vapour_Pressure_Deficit, slope_vapor_pressure_curve_soil, vapor_pressure_deficit)
            #     Actual_Soil_Evap = max(0, potential_evaporation_second_estimate)
                
                
            #     Hourly_Actual_Soil_Evap.append(Actual_Soil_Evap)
            #     Hourly_Soil_Rad.append(soil_absorbed_radiation_second_estimate)
            #     Hourly_Air_Soil_temp_dif.append(temperature_difference_soil)
            #     Hourly_boundary_layer_resistance_soil.append(boundary_layer_resistance_soil)
            #     Hourly_turbulence_resistance_soil.append(turbulence_resistance_soil)
            # Actual_Daily_Evaporation = Leaf.Leaf.aggregate_to_daily(Hourly_Actual_Soil_Evap, Day_Length)

            # if Actual_Daily_Evaporation > Maximum_possible_Evaporation :
            #     Water_Stress_Fraction=Maximum_possible_Evaporation/Actual_Daily_Evaporation
            #     potential_evaporation_soil_estimate=np.array(Hourly_Actual_Soil_Evap)*Water_Stress_Fraction
            #     self.Hourly_Actual_Soil_Evap=potential_evaporation_soil_estimate
            #     self.Actual_Daily_Evaporation = Leaf.Leaf.aggregate_to_daily(potential_evaporation_soil_estimate, Day_Length)

            #     # if (Water_Stress_Fraction <0.95) and (Water_Stress_Fraction>0):
            #     #     print(Water_Stress_Fraction)
            #     #     print('warning ..... Soil evaporation convergence issue!!!!!!!!!')
            # else:
            #     self.Hourly_Actual_Soil_Evap=Hourly_Actual_Soil_Evap
            #     self.Actual_Daily_Evaporation = Actual_Daily_Evaporation


            # self.Hourly_Soil_Rad=Hourly_Soil_Rad
            # self.Hourly_Air_Soil_temp_dif=Hourly_Air_Soil_temp_dif
            # self.Hourly_boundary_layer_resistance_soil=Hourly_boundary_layer_resistance_soil
            # self.Hourly_turbulence_resistance_soil=Hourly_turbulence_resistance_soil
            
            # # Aggregate Hourly data back to daily totals
            # self.potential_SoilRad_daily = Leaf.Leaf.aggregate_to_daily(Hourly_Soil_Rad,  Day_Length)
            # self.Day_Air_Soil_temp_dif=(Hourly_Air_Soil_temp_dif * wgauss).sum()
    

  





    def Calculate_Soil_Temperature(self, tmin, tmax):
        # Calculate the daily average temperature with a weighted average favoring tmax
        daily_avg_temperature = 0.29 * tmin + 0.71 * tmax
        
        # Calculate the nightly average temperature with a weighted average favoring tmin
        nightly_avg_temperature = 0.71 * tmin + 0.29 * tmax
        
        # Calculate the soil's average steady state temperature
        average_soil_temperature = (daily_avg_temperature + self.Day_Air_Soil_temp_dif + nightly_avg_temperature) / 2.
        
        # Calculate the rate of change in soil temperature
        rate_of_change_soil_temp = (average_soil_temperature - self.Soil_Temperature) / self.temperature_change_T
        
        
        # Update class attributes with the new calculated values
        self.average_soil_temperature = average_soil_temperature
        self.rate_of_change_soil_temp = rate_of_change_soil_temp




    def Calculate_Water_Balance(self,irrigation_amount,Potential_Canopy_Transpiration,Potential_evaporation,Root_Depth):
        _, _, _, _,_,soil_moisture_array=self.run_soil_water_model(1,irrigation_amount,  0, 0, Root_Depth  )

        # print(soil_moisture_array)
        plantwc, evapowc, Root_zone_avg_water_content, Root_zone_avg_soil_moisture,water_content_array,soil_moisture_array=self.run_soil_water_model(1,irrigation_amount,  Potential_Canopy_Transpiration, Potential_evaporation, Root_Depth,  )
        # print(784,soil_moisture_array,irrigation_amount,
        #       self.Actual_Daily_Evaporation, self.water_supply_for_evaporation,
        #       Actual_Canopy_Transpiration,self.water_supply_for_Transpiration,Root_Depth)


        total_water_soil_profile=0
        total_depth=0
        
        for j, layer in enumerate(self.layers):
            if j==0:
                layer['vwc']=soil_moisture_array[j]
            else:
                layer['vwc']=soil_moisture_array[j]
                total_water_soil_profile+=layer['vwc']*layer['thick']
                total_depth+=layer['thick']
                
        average_soil_moisture_soil_profile=total_water_soil_profile/total_depth
        below_root_zone_water_content=total_water_soil_profile-Root_zone_avg_water_content/1000

        
        self.average_below_root_zone_moisture=below_root_zone_water_content/(self.total_depth-Root_Depth-self.Evaporative_Depth)
        self.average_root_zone_water_content=Root_zone_avg_water_content #mm
        self.water_supply_for_Transpiration=plantwc #mm
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

     