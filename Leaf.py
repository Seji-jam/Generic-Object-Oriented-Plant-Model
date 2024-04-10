import numpy as np
import math
import Canopy

def SWITCH_FUN(x, y1, y2):
    return y1 if x < 0 else y2


##self.Plant_Height ... take care of Plant_Height because it's should be comming from canopy
class Leaf: 
    def __init__(self, Veg_C_Fraction, Spec_Leaf_Area, LAI_ini, Leaf_Blade_Angle, Leaf_Width, Tot_Leaf_N, Min_Spec_Leaf_N, Pathway_C3C4, Ambient_CO2, Activation_Energy_JMAX, VCMAX_LeafN_Slope, JMAX_LeafN_Slope, Photosynthetic_Light_Response_Factor):
        self.Veg_C_Fraction = Veg_C_Fraction
        self.Spec_Leaf_Area = Spec_Leaf_Area
        self.LAI_ini = LAI_ini
        self.Leaf_Blade_Angle = Leaf_Blade_Angle
        self.Leaf_Width = Leaf_Width
        self.Tot_Leaf_N = Tot_Leaf_N
        self.Min_Spec_Leaf_N = Min_Spec_Leaf_N
        self.Pathway_C3C4 = Pathway_C3C4
        self.Ambient_CO2 = Ambient_CO2
        self.Activation_Energy_JMAX = Activation_Energy_JMAX
        self.VCMAX_LeafN_Slope = VCMAX_LeafN_Slope
        self.JMAX_LeafN_Slope = JMAX_LeafN_Slope
        self.Photosynthetic_Light_Response_Factor = Photosynthetic_Light_Response_Factor

        # Initial conditions
        self.Carbon_Dead_Leaves=self.Carbon_Litters_Soil = self.Dif_Ait_leaf_T=0         

        # Initialize attributes to store outputs
        self.leaf_area_output = {}
        self.specific_leaf_n_output = {}
               
        
        self.hourly_Air_SU_leaf_T_diff=[]
        self.hourly_Air_SH_leaf_T_diff=[]

        
        self.hourly_SU_leaf_T=[]
        self.hourly_SH_leaf_T=[]
        
        self.hourly_photosyn_SU=[]
        self.hourly_darkresp_SU=[]
        self.hourly_transpiration_SU=[]
        self.hourly_Aradiation_SU=[]
        self.hourly_rsw_SU=[]
        self.hourly_slope_SU=[]
        self.hourly_apar_SU=[]

        self.Hourly_Actual_Photosynthesis_SU=[]
        self.Hourly_Actual_Transpiration_SU=[]
        self.hourly_Actual_Air_SU_leaf_T_diff=[]
        self.hourly_Actual_SU_leaf_T=[]
        
        
        self.hourly_photosyn_SH=[]
        self.hourly_darkresp_SH=[]
        self.hourly_transpiration_SH=[]
        self.hourly_Aradiation_SH=[]
        self.hourly_rsw_SH=[]
        self.hourly_slope_SH=[]
        self.hourly_apar_SH=[]

        
        self.Hourly_Actual_Photosynthesis_SH=[]
        self.Hourly_Actual_Transpiration_SH=[]
        self.hourly_Actual_Air_SH_leaf_T_diff=[]
        self.hourly_Actual_SH_leaf_T=[]

        self.potential_photosyn_SU_daily=0
        self.dark_rsp_SU_daily=0
        self.potential_photosyn_SH_daily=0
        self.dark_rsp_SH_daily=0
        self.potential_transpiration_SU_daily=0
        self.potential_transpiration_SH_daily=0
        self.Apar_SH_daily=0  
        self.Apar_SU_daily=0  


    
    

    
    def KDR_Coeff(Solar_Elev_Sin, Leaf_Blade_Angle):

        # Convert solar elevation sine to angle in radians
        Solar_Elev_Angle = np.arcsin(Solar_Elev_Sin)

        # Calculate average projection of leaves in the direction of solar beam
        if Solar_Elev_Sin >= np.sin(Leaf_Blade_Angle):
            Leaf_Orientation_Avg = Solar_Elev_Sin * np.cos(Leaf_Blade_Angle)
        else:
            Leaf_Orientation_Avg = (2 / np.pi) * (Solar_Elev_Sin * np.cos(Leaf_Blade_Angle) * np.arcsin(np.tan(Solar_Elev_Angle) / np.tan(Leaf_Blade_Angle)) + ((np.sin(Leaf_Blade_Angle))**2 - Solar_Elev_Sin**2)**0.5)

        # Calculate beam radiation extinction coefficient
        Direct_Beam_Ext_Coeff = Leaf_Orientation_Avg / Solar_Elev_Sin

        return Direct_Beam_Ext_Coeff

    def KDF_Coeff(Initial_LAI, Leaf_Blade_Angle, Scattering_Coeff):

        # Extinction coefficients for direct light at 15, 45, and 75 degrees elevation
        Beam_Ext_Coeff_15 = Leaf.KDR_Coeff(np.sin(15. * np.pi / 180.), Leaf_Blade_Angle)
        Beam_Ext_Coeff_45 = Leaf.KDR_Coeff(np.sin(45. * np.pi / 180.), Leaf_Blade_Angle)
        Beam_Ext_Coeff_75 = Leaf.KDR_Coeff(np.sin(75. * np.pi / 180.), Leaf_Blade_Angle)

        # Calculate diffuse light extinction coefficient
        Diffuse_Ext_Coeff = -1 / Initial_LAI * np.log(0.178 * np.exp(-Beam_Ext_Coeff_15 * (1.0 - Scattering_Coeff)**0.5 * Initial_LAI) +
                            0.514 * np.exp(-Beam_Ext_Coeff_45 * (1.0 - Scattering_Coeff)**0.5 * Initial_LAI) +
                            0.308 * np.exp(-Beam_Ext_Coeff_75 * (1.0 - Scattering_Coeff)**0.5 * Initial_LAI))

        return Diffuse_Ext_Coeff


    def Update_Leaf_Area(self):

        # Calculating delta and total leaf area index (LAI)
        Delta_LAI = (self.Carb_Dead_Leaves - self.Carb_Leaves_to_Litters) / self.Veg_C_Fraction * self.Spec_Leaf_Area
        Total_LAI = self.LAI_ini + self.Carb_Dead_Leaves / self.Veg_C_Fraction * self.Spec_Leaf_Area

        # Nitrogen and light extinction coefficients
        Light_Ext_Coeff = Leaf.KDF_Coeff(Total_LAI, self.Leaf_Blade_Angle * np.pi / 180., 0.2)
        Nitro_Ext_Coeff = Light_Ext_Coeff * (self.Tot_Leaf_N - self.Min_Spec_Leaf_N * Total_LAI)
        Nitro_Base_K = self.Min_Spec_Leaf_N * (1.0 - np.exp(-Light_Ext_Coeff * Total_LAI))
        
        # Assuming wind extinction coefficient is similar to light for simplicity
        Wind_Ext_Coeff = Light_Ext_Coeff  
        Leaf_Nitro_Ext_Coeff = 1.0 / Total_LAI * math.log((Nitro_Ext_Coeff + Nitro_Base_K) / (Nitro_Ext_Coeff * math.exp(-Light_Ext_Coeff * Total_LAI) + Nitro_Base_K))

        # Integrating LAI considering nitrogen effect
        N_determined_LAI = math.log(1. + Leaf_Nitro_Ext_Coeff * max(0., self.Tot_Leaf_N) / self.Min_Spec_Leaf_N) / Leaf_Nitro_Ext_Coeff
        Final_LAI = min(N_determined_LAI, self.LAI_ini)
        
        # Updating leaf area outputs
        self.leaf_area_output = {
            'Delta_LAI': Delta_LAI,
            'Total_LAI': Total_LAI,
            'Light_Ext_Coeff': Light_Ext_Coeff,
            'Nitro_Ext_Coeff': Nitro_Ext_Coeff,
            'Nitro_Base_K': Nitro_Base_K,
            'Wind_Ext_Coeff': Wind_Ext_Coeff,
            'Leaf_Nitro_Ext_Coeff': Leaf_Nitro_Ext_Coeff,
            'Final_LAI': Final_LAI,
            'N_determined_LAI': N_determined_LAI
        }

    def Update_Specific_Leaf_N(self):
        # Extracting final LAI and leaf nitrogen extinction coefficient
        Final_LAI = self.leaf_area_output['Final_LAI']
        Leaf_Nitro_Ext_Coeff = self.leaf_area_output['Leaf_Nitro_Ext_Coeff']
        
        # Calculating specific leaf nitrogen
        Spec_Leaf_N = self.Tot_Leaf_N / Final_LAI  # Average specific leaf nitrogen
        Spec_Leaf_N_Top = self.Tot_Leaf_N * Leaf_Nitro_Ext_Coeff / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Final_LAI))  # For top leaves
        Spec_Leaf_N_Base_Calc = self.Tot_Leaf_N * Leaf_Nitro_Ext_Coeff * np.exp(-Leaf_Nitro_Ext_Coeff * Final_LAI) / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Final_LAI))  # Exponential nitrogen profile
        Spec_Leaf_N_Top_Increment = (self.Tot_Leaf_N + 0.001 * self.Tot_Leaf_N) * Leaf_Nitro_Ext_Coeff / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Final_LAI))  # With small nitrogen increment

        # Updating specific leaf nitrogen outputs
        self.specific_leaf_n_output = {
            'Spec_Leaf_N': Spec_Leaf_N,
            'Spec_Leaf_N_Top': Spec_Leaf_N_Top,
            'Spec_Leaf_N_Base_Calc': Spec_Leaf_N_Base_Calc,
            'Spec_Leaf_N_Top_Increment': Spec_Leaf_N_Top_Increment
        }

    def aggregate_to_daily(hourly_results, Day_Length):
       
        wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])

        var_daily = 0
        
        for i, var_hourly in enumerate(hourly_results):
            var_daily += var_hourly * wgauss[i]
            # print(var_daily)
        # Convert to daily totals
        var_daily *= Day_Length * 3600
    
        return var_daily
    
    def convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind):

         Gaussian_Points = 5
         x_Gauss = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
         w_Gauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
         
         hourly_data = []
         for i in range(Gaussian_Points):
             hour = 12 - 0.5 * Day_Length + Day_Length * x_Gauss[i]
             Solar_Elevation_Sin = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (hour - 12.) / 24.))
             Diffuse_Ratio = Solar_Radiation * (Solar_Elevation_Sin * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
             Hourly_Temp = Min_Temp + (Max_Temp - Min_Temp) * np.sin(np.pi * (hour + Day_Length / 2 - 12) / (Day_Length + 3))
             
             hourly_data.append((Solar_Constant, Hourly_Temp, Solar_Elevation_Sin, Diffuse_Ratio, Wind))
         
         return hourly_data, w_Gauss
    
    def INTERNAL_CO2(Leaf_Temp, VPD, VPD_Slope, Ambient_CO2, Pathway_C3C4):

        # Air-to-leaf vapor pressure deficit
        Saturated_Vapor_Pressure_Leaf = 0.611 * np.exp(17.4 * Leaf_Temp / (Leaf_Temp + 239.))
        Vapor_Pressure_Deficit_Leaf = max(0, Saturated_Vapor_Pressure_Leaf - VPD)
        
        # Constants based on crop type
        Michaelis_Menten_CO2_25C = 650 if Pathway_C3C4 == 1 else 404.9
        Michaelis_Menten_O2_25C = 450 if Pathway_C3C4 == 1 else 278.4
        
        # Adjustment for temperature
        KMC = Michaelis_Menten_CO2_25C * np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 79430 / 8.314)
        KMO = Michaelis_Menten_O2_25C * np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 36380 / 8.314)
        
        # CO2 compensation point without dark respiration
        GAMMAX = 0.5 * np.exp(-3.3801 + 5220./(Leaf_Temp + 273.) / 8.314) * 210 * KMC / KMO
        RDVCX = 0.0089 * np.exp((1/298 - 1/(Leaf_Temp + 273)) * (46390 - 65330) / 8.314)
        GAMMA = (GAMMAX + RDVCX * KMC * (1 + 210 / KMO)) / (1 - RDVCX) / (10 if Pathway_C3C4 == 1 else 1)
        Intercellular_CO2_Ratio = 1 - (1 - GAMMA / Ambient_CO2) * (0.14 + VPD_Slope * Vapor_Pressure_Deficit_Leaf)
        
        # Intercellular CO2 concentration
        Intercellular_CO2 = Intercellular_CO2_Ratio * Ambient_CO2
        
        return Saturated_Vapor_Pressure_Leaf, Intercellular_CO2
      
        
            
        
    def LIGHT_ABSORB(Scattering_Coeff, Direct_Beam_Ext_Coeff, Scattered_Beam_Ext_Coeff, Diffuse_Ext_Coeff, Canopy_Beam_Reflect_Coeff, Canopy_Diffuse_Reflect_Coeff, Incident_Direct_Beam_Rad, Incident_Diffuse_Rad, Leaf_Area_Index):

        Total_Canopy_Absorbed_Light = (1. - Canopy_Beam_Reflect_Coeff) * Incident_Direct_Beam_Rad * (1. - np.exp(-Scattered_Beam_Ext_Coeff * Leaf_Area_Index)) + (1. - Canopy_Diffuse_Reflect_Coeff) * Incident_Diffuse_Rad * (1. - np.exp(-Diffuse_Ext_Coeff * Leaf_Area_Index))
        
        Absorbed_Sunlit_Rad = (1 - Scattering_Coeff) * Incident_Direct_Beam_Rad * (1 - np.exp(-Direct_Beam_Ext_Coeff * Leaf_Area_Index)) \
            + (1 - Canopy_Diffuse_Reflect_Coeff) * Incident_Diffuse_Rad / (Diffuse_Ext_Coeff + Direct_Beam_Ext_Coeff) * Diffuse_Ext_Coeff * (1 - np.exp(-(Diffuse_Ext_Coeff + Direct_Beam_Ext_Coeff) * Leaf_Area_Index)) \
            + Incident_Direct_Beam_Rad * ((1 - Canopy_Beam_Reflect_Coeff) / (Scattered_Beam_Ext_Coeff + Direct_Beam_Ext_Coeff) * Scattered_Beam_Ext_Coeff * (1 - np.exp(-(Scattered_Beam_Ext_Coeff + Direct_Beam_Ext_Coeff) * Leaf_Area_Index)) \
                                          - (1 - Scattering_Coeff) * (1 - np.exp(-2 * Direct_Beam_Ext_Coeff * Leaf_Area_Index)) / 2)
        
        Absorbed_Shaded_Rad = Total_Canopy_Absorbed_Light - Absorbed_Sunlit_Rad
        
        return Absorbed_Sunlit_Rad, Absorbed_Shaded_Rad
    
    def PHOTOSYN(Pathway_C3C4, Absorbed_PAR, Leaf_Temp, Intercellular_CO2, Leaf_N_Content, Activation_Energy_Jmax, Vcmax_N_Slope, Jmax_N_Slope, Light_Response_Convexity):
       
        # Constants and initial calculations for both C3 and C4 plants, adapted for temperature
        O2_Concentration = 210  # Oxygen concentration (mmol/mol)
        Kmc_25C = 650 if Pathway_C3C4 < 0 else 404.9  # Michaelis-Menten constant for CO2 at 25°C, adjusted for C3/C4
        Kmo_25C = 450 if Pathway_C3C4 < 0 else 278.4  # Michaelis-Menten constant for O2 at 25°C, adjusted for C3/C4
        
        # Calculating temperature-dependent parameters
        Kmc = Kmc_25C * np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 79430 / 8.314)
        Kmo = Kmo_25C * np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 36380 / 8.314)
        CO2_Compensation_Point = 0.5 * np.exp(-3.3801 + 5220. / (Leaf_Temp + 273.) / 8.314) * O2_Concentration * Kmc / Kmo
        
        Vcmax_Temperature_Adjustment = np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 65330 / 8.314)
        Electron_Transport_Temperature_Adjustment = np.exp((1./298. - 1./(Leaf_Temp + 273.)) * Activation_Energy_Jmax / 8.314)
        
        Vcmax = Vcmax_N_Slope * Vcmax_Temperature_Adjustment * Leaf_N_Content
        Jmax = Jmax_N_Slope * Electron_Transport_Temperature_Adjustment * Leaf_N_Content
        
        # Calculating photosynthesis rates based on available light and nitrogen content
        PAR_Photon_Flux = 4.56 * Absorbed_PAR  # Conversion factor from J to umol for PAR
        
        # Further calculations for photosynthesis rates, dark respiration, and adjustments for plant type
        # Omitted for brevity, but would include detailed biochemical model calculations as hinted above
        
        # Placeholder return values; real implementation would include the actual calculations
        Gross_Leaf_Photosynthesis = 0  # This would be calculated based on the model
        Leaf_Dark_Respiration = 0  # This would be calculated based on temperature adjustments and leaf N content
        
        return Gross_Leaf_Photosynthesis, Leaf_Dark_Respiration

    def Penman_Monteith(Stomatal_Resist_Water, Turbulence_Resist, Boundary_Layer_Resist_Water, 
                        Boundary_Layer_Resist_Heat, Absorbed_Global_Radiation, Atmospheric_Transmissivity, 
                        Fraction_Leaf_Classes, Leaf_Temp, Vapour_Pressure, Saturated_Vapour_Pressure_Slope, VPD):

     # Constants
     Stefan_Boltzmann_Const = 5.668E-8  # Stefan-Boltzmann constant (J/m2/s/K4)
     Latent_Heat_Vaporization = 2.4E6  # Latent heat of water vaporization (J/kg)
     Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity (J/m3/°C)
     Psychrometric_Constant = 0.067  # Psychrometric constant (kPa/°C)
    
     # Net absorbed radiation calculation
     Sky_Clearness = max(0, min(1, (Atmospheric_Transmissivity - 0.25) / 0.45))  # Sky clearness
     Black_Body_Radiation = Stefan_Boltzmann_Const * (Leaf_Temp + 273)**4
     Longwave_Net_Radiation = Black_Body_Radiation * (0.56 - 0.079 * np.sqrt(Vapour_Pressure * 10)) * (0.1 + 0.9 * Sky_Clearness) * Fraction_Leaf_Classes
    
     Net_Leaf_Absorbed_Radiation = Absorbed_Global_Radiation - Longwave_Net_Radiation
    
     # Intermediate variable for resistances
     Resistance_Sum = Psychrometric_Constant * (Boundary_Layer_Resist_Water + Turbulence_Resist + Stomatal_Resist_Water) / (Boundary_Layer_Resist_Heat + Turbulence_Resist)
    
     # Radiation-determined transpiration component
     Radiation_Component = Net_Leaf_Absorbed_Radiation * Saturated_Vapour_Pressure_Slope / (Saturated_Vapour_Pressure_Slope + Resistance_Sum) / Latent_Heat_Vaporization
    
     # Vapour pressure-determined transpiration component
     Vapour_Pressure_Component = (Volumetric_Heat_Capacity_Air * VPD / (Boundary_Layer_Resist_Heat + Turbulence_Resist)) / (Saturated_Vapour_Pressure_Slope + Resistance_Sum) / Latent_Heat_Vaporization
    
     # Potential leaf transpiration calculation
     Potential_Leaf_Transpiration = max(1E-10, Radiation_Component + Vapour_Pressure_Component)
     
     return Potential_Leaf_Transpiration, Net_Leaf_Absorbed_Radiation
    
    def REFLECTION_Coeff(Leaf_Scattering_Coeff, Direct_Beam_Ext_Coeff):
        """
        Calculates reflection coefficients for beam radiation.
    
        Parameters:
        - Leaf_Scattering_Coeff: Leaf scattering coefficient.
        - Direct_Beam_Ext_Coeff: Direct beam radiation extinction coefficient.
    
        Returns:
        - Scattered_Beam_Ext_Coeff: Scattered beam radiation extinction coefficient.
        - Canopy_Beam_Reflect_Coeff: Canopy beam radiation reflection coefficient.
        """
        # Scattered beam radiation extinction coefficient calculation
        Scattered_Beam_Ext_Coeff = Direct_Beam_Ext_Coeff * (1 - Leaf_Scattering_Coeff)**0.5
    
        # Canopy reflection coefficient for horizontal leaves calculation
        Horizontal_Leaf_Phase_Function = (1 - (1 - Leaf_Scattering_Coeff)**0.5) / (1 + (1 - Leaf_Scattering_Coeff)**0.5)
        
        # Canopy beam radiation reflection coefficient calculation
        Canopy_Beam_Reflect_Coeff = 1 - np.exp(-2 * Horizontal_Leaf_Phase_Function * Direct_Beam_Ext_Coeff / (1 + Direct_Beam_Ext_Coeff))
        
        return Scattered_Beam_Ext_Coeff, Canopy_Beam_Reflect_Coeff


    def Avoid_Zero_Division(self,x):

        if x != 0:
            out = x
        else:
            out = 1.0
        return out
    
    def Update_State_Variables(self,Rate_LeafAreaIndex):
        self.LAI_ini += Rate_LeafAreaIndex


    def Limit_Function(xl, xh, x):
        if xl <= x <= xh:
            out = x
        elif x > xh:
            out = xh
        else:
            out = xl
        return out





class Leaf_sunlit(Leaf):
    def __init__(self, leaf_object):
        self.leaf_object = leaf_object
        
    def Calculate_leaf_temp(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind, Plant_Height):
        # Using updated attributes
        Total_LAI = self.leaf_object.leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.leaf_object.leaf_area_output['Final_LAI']
        Wind_Ext_Coeff = self.leaf_object.leaf_area_output['Ext_Coeff_Wind']
        Leaf_Nitro_Ext_Coeff = self.leaf_object.leaf_area_output['Leaf_Nitro_Ext_Coeff']
        Spec_Leaf_N_Top = self.leaf_object.specific_leaf_n_output['Spec_Leaf_N_Top']

        hourly_data, w_Gauss = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind)
        hourly_Sunlit_Leaf_Temp = []
        hourly_Air_Leaf_Temp_Diff = []

        for Solar_Constant, Hourly_Temp, Sin_Beam, Diffuse_Ratio, Wind_Speed in hourly_data:
            
            Incoming_PAR = 0.5 * Solar_Radiation
            Incoming_NIR = 0.5 * Solar_Radiation

            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)

            # Calculation of diffuse light fraction based on atmospheric transmissivity
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity

            # Ensuring a minimum threshold for the diffuse light fraction
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
            
            # Calculation for Diffuse and Direct PAR (Photosynthetically Active Radiation)
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR
            
            # Calculation for Diffuse and Direct NIR (Near-Infrared Radiation)
            Diffuse_NIR = Incoming_NIR * Diffuse_Light_Fraction
            Direct_NIR = Incoming_NIR - Diffuse_NIR
            
            # Calculating Extinction and Reflection Coefficients
            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
            Scattering_Coefficient_NIR = 0.8  # Leaf scattering coefficient for NIR
            
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
            
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            Scattered_Beam_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient)
            
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057  # Canopy diffuse PAR reflection coefficient
            Canopy_Diffuse_Reflection_Coefficient_NIR = 0.389  # Canopy diffuse NIR reflection coefficient
            
                        
            # Calculating photosynthetic nitrogen for sunlit canopy parts
            Spec_Leaf_N_Sunlit = Spec_Leaf_N_Top * (1. - np.exp(-(Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) - self.leaf_object.Min_Spec_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            
            # Absorption of PAR and NIR by sunlit leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            Absorbed_NIR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_NIR, Diffuse_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR, Canopy_Diffuse_Reflection_Coefficient_NIR, Direct_NIR, Diffuse_NIR, Leaf_Area_Index)
            
            # Calculating potential photosynthesis and dark respiration for sunlit leaves
            Potential_Photosynthesis_Sunlit, Dark_Respiration_Sunlit = Leaf.PHOTOSYN(self.leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Hourly_Temp, Intercellular_CO2, Spec_Leaf_N_Sunlit, self.leaf_object.Activation_Energy_Jmax, self.leaf_object.Vcmax_LeafN_Slope, self.leaf_object.Jmax_LeafN_Slope, self.leaf_object.Light_Response_Factor)
            
            # Total absorbed radiation (PAR + NIR) by sunlit leaves
            Total_Absorbed_Radiation_Sunlit = Absorbed_PAR_Sunlit + Absorbed_NIR_Sunlit
            
            # Calculating the fraction of sunlit and shaded canopy components
            Sunlit_Fraction = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            
            # Boundary layer resistance for heat and water vapor for sunlit leaves
            Boundary_Layer_Resistance_Heat_Sunlit = 1. / (0.01 * np.sqrt(Wind_Speed / self.leaf_object.Leaf_Width))
            Boundary_Layer_Resistance_Heat_Sunlit_Adjusted = (1 - np.exp(-(0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Resistance_Heat_Sunlit
            Boundary_Layer_Resistance_Water_Sunlit = 0.93 * Boundary_Layer_Resistance_Heat_Sunlit_Adjusted
            
            # Calculating potential CO2 conductance for sunlit leaves
            Potential_CO2_Conductance_Sunlit = (Potential_Photosynthesis_Sunlit - Dark_Respiration_Sunlit) * (273.15 + Hourly_Temp) / 0.53717 / (self.leaf_object.Ambient_CO2 - Intercellular_CO2)
            
            # Calculating Turbulence Resistance for the canopy
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            Stomatal_Resistance_Sunlit = max(1E-30, 1 / Potential_CO2_Conductance_Sunlit - Boundary_Layer_Resistance_Water_Sunlit * 1.3 - Turbulence_Resistance) / 1.6
            
            # Applying Penman-Monteith equation to calculate transpiration and net radiation absorbed
            Transpiration_Sunlit, Net_Radiation_Absorbed_Sunlit = Leaf.Penman_Monteith(Stomatal_Resistance_Sunlit, Turbulence_Resistance, Boundary_Layer_Resistance_Water_Sunlit, Boundary_Layer_Resistance_Heat_Sunlit_Adjusted, Total_Absorbed_Radiation_Sunlit, Atmospheric_Transmissivity, Sunlit_Fraction, Hourly_Temp, Vapour_Pressure, SVP_Slope, Vapour_Pressure_Deficit)
            
            # Calculating leaf temperature difference
            Latent_Heat_Vaporization = 2.4E6  # Latent heat of vaporization (J/kg)
            Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity of air (J/m^3/°C)
            
            Air_Leaf_Temp_Diff_Sunlit = (Net_Radiation_Absorbed_Sunlit - Latent_Heat_Vaporization * Transpiration_Sunlit) * (Boundary_Layer_Resistance_Heat_Sunlit_Adjusted + (Turbulence_Resistance * Sunlit_Fraction)) / Volumetric_Heat_Capacity_Air
            
            # Appending calculated temperature difference and adjusted leaf temperature to lists
            hourly_Air_Leaf_Temp_Diff.append(Air_Leaf_Temp_Diff_Sunlit)
            Adjusted_Leaf_Temp_Sunlit = Hourly_Temp + Air_Leaf_Temp_Diff_Sunlit
            hourly_Sunlit_Leaf_Temp.append(Adjusted_Leaf_Temp_Sunlit)
        self.leaf_object.hourly_SU_leaf_T=hourly_Sunlit_Leaf_Temp
        self.leaf_object.hourly_Air_SU_leaf_T_diff=hourly_Air_Leaf_Temp_Diff             
        
    def Calculate_Potential_Photosynthesis(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Noon_Max):
        # Use updated attributes
        Total_LAI_ini = self.leaf_object.leaf_area_output['Total_LAI_ini']
        Final_LAI = self.leaf_object.leaf_area_output['Final_LAI']
        Leaf_Nitro_Ext_Coeff = self.leaf_object.leaf_area_output['KN']
        Specific_Leaf_N_Top = self.leaf_object.specific_leaf_n_output['Spec_Leaf_N_Top']
    
        # Convert daily data to hourly
        hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Noon_Max)
        hourly_Sunlit_Leaf_Temp = self.leaf_object.hourly_Sunlit_Leaf_Temp
        
        hourly_Photosynthesis = []
        hourly_Dark_Respiration = []
    
        for (_, Hourly_Temp, Sin_Beam, Diffuse_Ratio, Wind_Speed), Sunlit_Leaf_Temp in zip(hourly_data, hourly_Sunlit_Leaf_Temp):
            # Use sunlit leaf temperature if it differs from ambient

    
            Incoming_PAR = 0.5 * Solar_Radiation
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)
    
            # Calculate diffuse light fraction
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity
    
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
    
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR   
            
            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
            
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI_ini, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057  # Canopy diffuse PAR reflection coefficient
            
            # Adjusting for vapor pressure deficit influence on intercellular CO2 concentration
            Vapor_Pressure_Deficit_Response =  0.116214  # Slope for linear effect of VPD on Ci/Ca
            Sat_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Sunlit_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            
            # Calculating photosynthetic nitrogen availability for sunlit canopy parts
            Photosynthetic_Nitrogen_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Final_LAI)) / (Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) - self.leaf_object.Min_Spec_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Final_LAI)) / Direct_Beam_Extinction_Coefficient
            
            # Absorption of PAR by sunlit leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Final_LAI)
            
            # Calculating potential photosynthesis and dark respiration for sunlit leaves
            Potential_Photosynthesis_Sunlit, Dark_Respiration_Sunlit = Leaf.PHOTOSYN(self.leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Sunlit_Leaf_Temp, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Sunlit, self.leaf_object.Activation_Energy_Jmax, self.leaf_object.Vcmax_LeafN_Slope, self.leaf_object.Jmax_LeafN_Slope, self.leaf_object.Light_Response_Factor)
            
            # Logging potential photosynthesis for sunlit leaves
            print("Potential Photosynthesis (Sunlit):", Potential_Photosynthesis_Sunlit)
            
            # Appending the calculated photosynthesis and respiration rates to their respective hourly lists
            hourly_Photosynthesis.append(Potential_Photosynthesis_Sunlit)
            hourly_Dark_Respiration.append(Dark_Respiration_Sunlit)
            
            self.leaf_object.hourly_Photosynthesis_Sunlit = hourly_Photosynthesis
            self.leaf_object.hourly_Dark_Respiration_Sunlit = hourly_Dark_Respiration

        # Aggregate hourly data back to daily totals
        Daily_Potential_Photosynthesis_Sunlit = Leaf.aggregate_to_daily(hourly_Photosynthesis, Day_Length)
        Daily_Dark_Respiration_Sunlit = Leaf.aggregate_to_daily(hourly_Dark_Respiration, Day_Length)
        
        self.leaf_object.Daily_Potential_Photosynthesis_Sunlit = Daily_Potential_Photosynthesis_Sunlit
        self.leaf_object.Daily_Dark_Respiration_Sunlit = Daily_Dark_Respiration_Sunlit
                                    
                
        
    
    def Calculate_Potential_Transpiration(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Noon_Max, Plant_Height):
        # Utilizing updated attributes
        Total_LAI_ini = self.leaf_object.leaf_area_output['Total_LAI_ini']
        Final_LAI = self.leaf_object.leaf_area_output['Final_LAI']
        Ext_Coeff_Wind = self.leaf_object.leaf_area_output['KW']
    
        hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Noon_Max)
    
        hourly_Sunlit_Leaf_Temp = self.leaf_object.hourly_Sunlit_Leaf_Temp
        hourly_Photosynthesis_Sunlit = self.leaf_object.hourly_Photosynthesis_Sunlit
        hourly_Dark_Respiration_Sunlit = self.leaf_object.hourly_Dark_Respiration_Sunlit
    
        hourly_Transpiration_Sunlit = []
        hourly_Absorbed_Radiation_Sunlit = []
        hourly_Stomatal_Resistance_Water_Sunlit = []
        hourly_Slope_VPD_Sunlit = []
        hourly_Absorbed_PAR_Sunlit = []
    
        for (_, Hourly_Temp, Sin_Beam, Diffuse_Ratio, Wind_Speed), Sunlit_Leaf_Temp, Photosynthesis_Sunlit, Dark_Respiration_Sunlit in zip(hourly_data, hourly_Sunlit_Leaf_Temp, hourly_Photosynthesis_Sunlit, hourly_Dark_Respiration_Sunlit):
            Atmospheric_Transmissivity = Diffuse_Ratio / (0.5 * Solar_Constant * Sin_Beam)
            Diffuse_Light_Fraction = max(self.adjust_diffuse_light_fraction(Atmospheric_Transmissivity, Sin_Beam), 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
    
            Diffuse_PAR = 0.5 * Solar_Radiation * Diffuse_Light_Fraction
            Direct_PAR = 0.5 * Solar_Radiation - Diffuse_PAR
            NIR = 0.5 * Solar_Radiation
            
            # Calculation of Vapor Pressure Deficit effect
            Vapor_Pressure_Deficit_Response =  0.116214
            Sat_Vapor_Pressure, Intercellular_CO2 = Leaf.INTERNAL_CO2(Sunlit_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
    
            # Incoming Diffuse and Direct NIR (Near-Infrared Radiation)
            Diffuse_NIR = NIR * Diffuse_Light_Fraction
            Direct_NIR = NIR - Diffuse_NIR
            
            # Leaf Blade Angle and Direct Beam Extinction Coefficient
            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            # Scattering Coefficients for PAR and NIR
            Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
            Scattering_Coefficient_NIR = 0.8  # Leaf scattering coefficient for NIR
            
            # Diffuse Extinction Coefficients for PAR and NIR
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI_ini, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_LAI_ini, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
            
            # Reflection Coefficients for PAR and NIR
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            Scattered_Beam_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient)
            
            # Canopy Diffuse Reflection Coefficients for PAR and NIR
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057
            Canopy_Diffuse_Reflection_Coefficient_NIR = 0.389
            
            # Absorption of PAR and NIR by Sunlit Leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Final_LAI)
            Absorbed_NIR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_NIR, Diffuse_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR, Canopy_Diffuse_Reflection_Coefficient_NIR, Direct_NIR, Diffuse_NIR, Final_LAI)
            
            # Total Absorbed Radiation (PAR + NIR) by Sunlit Leaves
            Total_Absorbed_Radiation_Sunlit = Absorbed_PAR_Sunlit + Absorbed_NIR_Sunlit
            
            # Fraction of Sunlit and Shaded Components in Canopy
            Sunlit_Fraction = 1. / Direct_Beam_Extinction_Coefficient / Final_LAI * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Final_LAI))
            
            # Boundary Layer Resistance for Heat and Water Vapor for Sunlit Leaves
            Boundary_Layer_Resistance_Heat_Sunlit = 1. / (0.01 * np.sqrt(Wind_Speed / self.leaf_object.Leaf_Width))
            Boundary_Layer_Resistance_Heat_Sunlit_Adjusted = (1 - np.exp(-(0.5 * Ext_Coeff_Wind + Direct_Beam_Extinction_Coefficient) * Final_LAI)) / (0.5 * Ext_Coeff_Wind + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Resistance_Heat_Sunlit
            Boundary_Layer_Resistance_Water_Sunlit = 0.93 * Boundary_Layer_Resistance_Heat_Sunlit_Adjusted
            
                    
            Sat_Vapor_Pressure_Ambient, Intercellular_CO2_Ambient = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            Sat_Vapor_Pressure_Leaf, Intercellular_CO2_Leaf = Leaf.INTERNAL_CO2(Sunlit_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            
            Vapour_Pressure_Deficit = max(0, Sat_Vapor_Pressure_Ambient - Vapour_Pressure)
            
            # Dynamic slope calculation for the response of saturation vapor pressure to temperature
            Slope_VPD = (Sat_Vapor_Pressure_Leaf - Sat_Vapor_Pressure_Ambient) / self.Non_Zero(Sunlit_Leaf_Temp - Hourly_Temp)
            
            # Calculating potential CO2 conductance for sunlit leaves based on photosynthesis minus dark respiration
            Potential_CO2_Conductance_Sunlit = (Photosynthesis_Sunlit - Dark_Respiration_Sunlit) * (273.15 + Sunlit_Leaf_Temp) / 0.53717 / (self.leaf_object.Ambient_CO2 - Intercellular_CO2_Leaf)
            
            # Turbulence resistance calculation based on canopy structure and wind speed
            Turbulence_Resistance_Canopy = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            
            # Calculating sunlit leaf-specific stomatal resistance to water vapor
            Stomatal_Resistance_Water_Sunlit = max(1E-30, 1 / Potential_CO2_Conductance_Sunlit - Boundary_Layer_Resistance_Water_Sunlit * 1.3 - Turbulence_Resistance_Canopy * Sunlit_Fraction) / 1.6
            
            # Applying the Penman-Monteith equation for potential transpiration estimation
            Potential_Transpiration_Sunlit, Absorbed_Radiation_Sunlit = Leaf.Penman_Monteith(Stomatal_Resistance_Water_Sunlit, Turbulence_Resistance_Canopy * Sunlit_Fraction, Boundary_Layer_Resistance_Water_Sunlit, Boundary_Layer_Resistance_Heat_Sunlit_Adjusted, Total_Absorbed_Radiation_Sunlit, Atmospheric_Transmissivity, Sunlit_Fraction, Sunlit_Leaf_Temp, Vapour_Pressure, Slope_VPD, Vapour_Pressure_Deficit)
            
            # Updating the leaf object with calculated hourly and potential daily transpiration rates, absorbed radiation, etc.
            hourly_Transpiration_Sunlit.append(Potential_Transpiration_Sunlit)
            hourly_Absorbed_Radiation_Sunlit.append(Absorbed_Radiation_Sunlit)
            hourly_Stomatal_Resistance_Water_Sunlit.append(Stomatal_Resistance_Water_Sunlit)
            hourly_Slope_VPD_Sunlit.append(Slope_VPD)
            hourly_Absorbed_PAR_Sunlit.append(Absorbed_PAR_Sunlit)
            
        # Final aggregation to daily values might follow here
        self.leaf_object.Daily_Potential_Transpiration_Sunlit = Leaf.aggregate_to_daily(hourly_Transpiration_Sunlit, Day_Length)
        self.leaf_object.Daily_Absorbed_Radiation_Sunlit = Leaf.aggregate_to_daily(hourly_Absorbed_Radiation_Sunlit, Day_Length)
        
                
        self.leaf_object.hourly_Transpiration_Sunlit = hourly_Transpiration_Sunlit
        self.leaf_object.hourly_Absorbed_Radiation_Sunlit = hourly_Absorbed_Radiation_Sunlit
        self.leaf_object.hourly_Stomatal_Resistance_Water_Sunlit = hourly_Stomatal_Resistance_Water_Sunlit
        self.leaf_object.hourly_Slope_VPD_Sunlit = hourly_Slope_VPD_Sunlit
        self.leaf_object.hourly_Absorbed_PAR_Sunlit = hourly_Absorbed_PAR_Sunlit
        
        # Aggregate hourly data back to daily totals for a comprehensive daily overview
        Daily_Absorbed_Radiation_Sunlit = Leaf.aggregate_to_daily(hourly_Absorbed_Radiation_Sunlit, Day_Length)
        self.leaf_object.Daily_Potential_Transpiration_Sunlit = Leaf.aggregate_to_daily(hourly_Transpiration_Sunlit, Day_Length)
        # Optionally, aggregate other metrics like absorbed PAR if needed for further analysis
        # self.leaf_object.Daily_Absorbed_PAR_Sunlit = Leaf.aggregate_to_daily(hourly_Absorbed_PAR_Sunlit, Day_Length)
        
        
        
        
    def Update_LeafTemp_Photosynthesis_if_WaterStress(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Noon_Max, Plant_Height, Daily_Water_Supply, Soil_Depth_1, Root_Depth, hourly_Sunlit_Leaf_Temp, hourly_Shaded_Leaf_Temp, hourly_Soil_Evaporation):
        Gaussian_Points = 5
        Gaussian_Weights = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
    
        # Using updated attributes
        Total_LAI_ini = self.leaf_object.leaf_area_output['Total_LAI_ini']
        Final_LAI = self.leaf_object.leaf_area_output['Final_LAI']
        Leaf_Nitrogen_Extinction_Coefficient = self.leaf_object.leaf_area_output['KN']
        Wind_Extinction_Coefficient = self.leaf_object.leaf_area_output['KW']
        Spec_Leaf_N_Top = self.leaf_object.specific_leaf_n_output['Spec_Leaf_N_Top']
        Spec_Leaf_N_Bottom = self.leaf_object.specific_leaf_n_output['Spec_Leaf_N_Top_Increment']
    
        hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Noon_Max)
    
        Actual_Photosynthesis_list = []
        Actual_Transpiration_list = []
        Actual_Air_Leaf_Temperature_Difference_list = []
        Actual_Leaf_Temperature_list = []
    
        for i, (_, Hourly_Temp, Sin_Beam, Diffuse_Ratio, Wind_Speed), Transpiration_Sunlit, Transpiration_Shaded, Absorbed_Radiation_Sunlit, Stomatal_Resistance_Water_Sunlit, Slope_VPD, Soil_Evaporation in zip(range(Gaussian_Points), hourly_data, self.leaf_object.hourly_Transpiration_Sunlit, self.leaf_object.hourly_Transpiration_Shaded, self.leaf_object.hourly_Absorbed_Radiation_Sunlit, self.leaf_object.hourly_Stomatal_Resistance_Water_Sunlit, self.leaf_object.hourly_Slope_VPD_Sunlit, hourly_Soil_Evaporation):
    
            Hour = 12 - 0.5 * Day_Length + Day_Length * Gaussian_Weights[i]
            Sin_Beam_Adjusted = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (Hour - 12.) / 24.))
            
            # Calculating the hourly availability of soil water supply
            Water_Supply_Hourly = Daily_Water_Supply * (Sin_Beam_Adjusted * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
            
            # Water available for evaporation from the top soil layer
            Water_Supply_Evaporation = Water_Supply_Hourly * (Soil_Depth_1 / Root_Depth)
            
            # Total potential canopy transpiration for the hour
            Total_Potential_Canopy_Transpiration = Transpiration_Sunlit + Transpiration_Shaded
            

            # Calculation for the amount of transpiration from the top soil layer
            Transpiration_From_Top_Soil_Layer = Total_Potential_Canopy_Transpiration * (Soil_Depth_1 / Root_Depth)
            
            # Determining the maximum possible transpiration considering available water for evaporation and soil water supply
            Maximum_Possible_Transpiration = Transpiration_From_Top_Soil_Layer / (Transpiration_From_Top_Soil_Layer + Soil_Evaporation) * Water_Supply_Evaporation + (Water_Supply_Hourly - Water_Supply_Evaporation)

            # Decision-making based on the comparison between the maximum possible transpiration and the total potential canopy transpiration
            if Maximum_Possible_Transpiration < Total_Potential_Canopy_Transpiration:
                Actual_Canopy_Transpiration = Maximum_Possible_Transpiration
            else:
                Actual_Canopy_Transpiration = Total_Potential_Canopy_Transpiration
            
            # This ensures that the plant's transpiration rates are adjusted in real-time based on soil water availability, reflecting a realistic physiological response to varying water stress conditions. The model thus accounts for the dynamic interplay between environmental water availability a

            # Calculating the combined potential transpiration for the canopy
            Combined_Potential_Transpiration = Transpiration_Sunlit + Transpiration_Shaded
            
            # Adjusting actual transpiration for sunlit leaves based on canopy water availability
            Actual_Transpiration_Sunlit = Transpiration_Sunlit * (Actual_Canopy_Transpiration / Combined_Potential_Transpiration)
            
            # Solar radiation partitioning into PAR and determining the diffuse light fraction
            PAR = 0.5 * Daily_Sin_Beam_Exposure
            
            Atmospheric_Transmissivity = PAR / (0.5 * Solar_Constant * Sin_Beam)
            Diffuse_Light_Fraction = max([
                1 if Atmospheric_Transmissivity < 0.22 else
                1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2 if Atmospheric_Transmissivity <= 0.35 else
                1.47 - 1.66 * Atmospheric_Transmissivity,
                0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam))
            ])
            
            # Calculating direct and diffuse PAR based on atmospheric transmissivity
            Diffuse_PAR = PAR * Diffuse_Light_Fraction
            Direct_PAR = PAR - Diffuse_PAR
            
            # Adjusting for vapor pressure deficit's impact on intercellular CO2 concentration
            Vapor_Pressure_Deficit_Response =  0.116214
            
            # Extinction coefficients for sunlight penetration through the canopy
            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            # Scattering coefficient for PAR and adjustments for leaf and canopy level interactions
            Scattering_Coefficient_PAR = 0.2
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Final_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057
            
            # Calculating boundary layer resistance to heat and water vapor for sunlit leaves
            Boundary_Layer_Resistance_Heat_Flux = 0.01 * np.sqrt(Wind_Speed / self.leaf_object.Leaf_Width)
            Boundary_Layer_Resistance_Heat_Sunlit = (1 - np.exp(-(0.5 * Wind_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) * Final_LAI)) / (0.5 * Wind_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Resistance_Heat_Flux
            Boundary_Layer_Resistance_Water_Sunlit = 0.93 * Boundary_Layer_Resistance_Heat_Sunlit
            
            
            # Calculating the fraction of sunlit and shaded components in the canopy
            Sunlit_Fraction = 1.0 / Direct_Beam_Extinction_Coefficient / Final_LAI * (1 - np.exp(-Direct_Beam_Extinction_Coefficient * Final_LAI))
            
            # Adjusting turbulence resistance for the canopy
            Turbulence_Resistance_Canopy = 0.74 * (np.log((2 - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            
            # Calculating leaf temperature adjustments due to water stress
            Latent_Heat_Vaporization = 2.4E6  # Latent heat of vaporization (J/kg)
            Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity of air (J/m^3/°C)
            
            Adjusted_Turbulence_Resistance = Turbulence_Resistance_Canopy * Sunlit_Fraction
            Temperature_Difference = Leaf.Limit_Function(-25, 25, (Absorbed_Radiation_Sunlit - Latent_Heat_Vaporization * Actual_Transpiration_Sunlit) * (Boundary_Layer_Resistance_Heat_Sunlit + Adjusted_Turbulence_Resistance) / Volumetric_Heat_Capacity_Air)
            
            Adjusted_Leaf_Temperature = Hourly_Temp + Temperature_Difference
            
            Psychrometric_Constant = 0.067  # Psychrometric constant (kPa/°C)
            
            # Adjusting stomatal resistance to water under water stress conditions
            Adjusted_Stomatal_Resistance_Water = (Transpiration_Sunlit - Actual_Transpiration_Sunlit) * (Slope_VPD * (Boundary_Layer_Resistance_Heat_Sunlit + Adjusted_Turbulence_Resistance) + Psychrometric_Constant * (Boundary_Layer_Resistance_Water_Sunlit + Adjusted_Turbulence_Resistance)) / Actual_Transpiration_Sunlit / Psychrometric_Constant + Transpiration_Sunlit / Actual_Transpiration_Sunlit * Stomatal_Resistance_Water_Sunlit
            
            # Absorbed PAR calculation for sunlit leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Final_LAI)
            
            # Adjusting photosynthetic nitrogen for sunlit parts of the canopy
            Photosynthetic_Nitrogen_Sunlit = Spec_Leaf_N_Top * (1. - np.exp(-(Leaf_Nitrogen_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) * Final_LAI)) / (Leaf_Nitrogen_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) - self.leaf_object.Spec_Leaf_N_Min * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Final_LAI)) / Direct_Beam_Extinction_Coefficient
            
            # Calculating internal CO2 and photosynthesis adjustments
            Sat_Vapor_Pressure_Leaf, Intercellular_CO2_Leaf = Leaf.INTERNAL_CO2(Adjusted_Leaf_Temperature, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            Actual_Photosynthesis_Rate, Dark_Respiration_Rate = Leaf.PHOTOSYN(self.leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Adjusted_Leaf_Temperature, Intercellular_CO2_Leaf, Photosynthetic_Nitrogen_Sunlit, self.leaf_object.Activation_Energy_Jmax, self.leaf_object.VCMAX_LeafN_Slope, self.leaf_object.JMAX_LeafN_Slope, self.leaf_object.Photosynthetic_Light_Response_Factor)
            
            # Actual photosynthesis under water stress condition
            Actual_Photosynthesis = (1.6 * Stomatal_Resistance_Water_Sunlit + 1.3 * Boundary_Layer_Resistance_Water_Sunlit + Adjusted_Turbulence_Resistance) / (1.6 * Adjusted_Stomatal_Resistance_Water + 1.3 * Boundary_Layer_Resistance_Water_Sunlit + Adjusted_Turbulence_Resistance) * (Actual_Photosynthesis_Rate - Dark_Respiration_Rate) + Dark_Respiration_Rate
            
            # Appending results for analysis
            Actual_Photosynthesis_list.append(Actual_Photosynthesis)
            Actual_Transpiration_list.append(Actual_Transpiration_Sunlit)
            Actual_Air_Leaf_Temperature_Difference_list.append(Temperature_Difference)
            Actual_Leaf_Temperature_list.append(Adjusted_Leaf_Temperature)
        # Updating the leaf model object with the calculated hourly values
        self.leaf_object.hourly_Actual_Photosynthesis_Sunlit = Actual_Photosynthesis_list
        self.leaf_object.hourly_Actual_Transpiration_Sunlit = Actual_Transpiration_list
        self.leaf_object.hourly_Actual_AirLeafTempDiff_Sunlit = Actual_Air_Leaf_Temperature_Difference_list
        self.leaf_object.hourly_Actual_LeafTemp_Sunlit = Actual_Leaf_Temperature_list

    
    
class Leaf_Shaded(Leaf):
    def __init__(self, leaf_object):
        self.leaf_object = leaf_object
        
    def Calculate_leaf_temp(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind, Plant_Height):
        Total_LAI = self.leaf_object.leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.leaf_object.leaf_area_output['Final_LAI']
        Wind_Extinction_Coeff = self.leaf_object.leaf_area_output['Ext_Coeff_Wind']
        Leaf_N_Extinction_Coeff = self.leaf_object.leaf_area_output['Leaf_N_Ext_Coeff']
        
        Specific_Leaf_N_Top = self.leaf_object.specific_leaf_n_output['Specific_Leaf_N_Top']
        
        hourly_data, w_Gauss = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind)
        hourly_Air_Temperature_Difference_Shaded = []
        hourly_Shaded_Leaf_Temperature = []
        
        for Solar_Constant, Hourly_Temp, Sin_Beam, Solar_Radiation, Wind_Speed in hourly_data:
        
            Incoming_PAR = 0.5 * Solar_Radiation
            NIR = 0.5 * Solar_Radiation
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)
        
            Vapor_Pressure_Deficit_Response = 0.116214 # Slope for the linear effect of VPD on Ci/Ca
            Sat_Vapor_Pressure, Intercellular_CO2 = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2_Concentration, self.leaf_object.C3C4_Pathway)
        
            Vapour_Pressure_Deficit = max(0, Sat_Vapor_Pressure - Vapour_Pressure)
            Slope_SVP = 4158.6 * Sat_Vapor_Pressure / (Hourly_Temp + 239.) ** 2
        
        
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity

            
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
    
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR   
            
           
            Diffuse_NIR = NIR * Diffuse_Light_Fraction
            Direct_NIR = NIR - Diffuse_NIR


            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            
            Direct_Beam_Extinction_Coeff = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            Scattering_Coeff_PAR = 0.2  # Scattering coefficient for PAR
            Scattering_Coeff_NIR = 0.8  # Scattering coefficient for NIR
            
            Diffuse_Extinction_Coeff_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coeff_PAR)
            Diffuse_Extinction_Coeff_NIR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coeff_NIR)
            Scattered_Beam_Extinction_Coeff_PAR, Canopy_Beam_REF_Coeff_PAR = Leaf.REFLECTION_Coeff(Scattering_Coeff_PAR, Direct_Beam_Extinction_Coeff)
            Scattered_Beam_Extinction_Coeff_NIR, Canopy_Beam_REF_Coeff_NIR = Leaf.REFLECTION_Coeff(Scattering_Coeff_NIR, Direct_Beam_Extinction_Coeff)
            Canopy_Diffuse_REF_Coeff_PAR = 0.057  # Canopy reflection coefficient for diffuse PAR
            Canopy_Diffuse_REF_Coeff_NIR = 0.389  # Canopy reflection coefficient for diffuse NIR
            
            # Calculating absorbed PAR and NIR by shaded leaves
            _, Absorbed_PAR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coeff_PAR, Direct_Beam_Extinction_Coeff, Scattered_Beam_Extinction_Coeff_PAR, Diffuse_Extinction_Coeff_PAR, Canopy_Beam_REF_Coeff_PAR, Canopy_Diffuse_REF_Coeff_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            _, Absorbed_NIR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coeff_NIR, Direct_Beam_Extinction_Coeff, Scattered_Beam_Extinction_Coeff_NIR, Diffuse_Extinction_Coeff_NIR, Canopy_Beam_REF_Coeff_NIR, Canopy_Diffuse_REF_Coeff_NIR, Direct_NIR, Diffuse_NIR, Leaf_Area_Index)



            Absorbed_Total_Radiation_Shaded = Absorbed_PAR_Shaded + Absorbed_NIR_Shaded
            
            # Adjusting the slope for the vapor pressure deficit's effect on internal CO2 concentration
            Vapor_Pressure_Slope =  0.116214
            Saturated_Vapor_Pressure, Internal_CO2 = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Slope, self.leaf_object.Ambient_CO2, self.leaf_object.C3_C4)
            
            # Total photosynthetic nitrogen across the entire canopy
            Total_Photosynthetic_N_Canopy = (Specific_Leaf_N_Top * (1. - np.exp(-Leaf_N_Extinction_Coeff * Leaf_Area_Index)) / Leaf_N_Extinction_Coeff - self.leaf_object.Min_SLN * Leaf_Area_Index)
            
            # Photosynthetic nitrogen allocated to sunlit and shaded leaf portions
            N_Photosynthetic_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coeff) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coeff) - self.leaf_object.Min_SLN * (1. - np.exp(-Direct_Beam_Extinction_Coeff * Leaf_Area_Index)) / Direct_Beam_Extinction_Coeff
            N_Photosynthetic_Shaded = Total_Photosynthetic_N_Canopy - N_Photosynthetic_Sunlit
            
            Potential_Photosynthesis_Shaded, Dark_Respiration_Shaded = Leaf.PHOTOSYN(self.leaf_object.C3_C4, Absorbed_PAR_Shaded, Hourly_Temp, Internal_CO2, N_Photosynthetic_Shaded, self.leaf_object.Activation_Energy_JMAX, self.leaf_object.VCMAX_LeafN_Slope, self.leaf_object.JMAX_LeafN_Slope, self.leaf_object.Photosynthetic_Light_Response_Factor)
            
            # Calculating the fraction of sunlit and shaded leaf area
            Fraction_Sunlit = 1. / Direct_Beam_Extinction_Coeff / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coeff * Leaf_Area_Index))
            Fraction_Shaded = 1 - Fraction_Sunlit

            # Boundary layer resistance for canopy, sunlit and shaded leaves
            Boundary_Layer_Leaf_Flow = 0.01 * np.sqrt(Wind_Speed / self.leaf_object.Leaf_Width)
            Canopy_Boundary_Layer_Flow = (1 - np.exp(-0.5 * Wind_Extinction_Coeff * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff) * Boundary_Layer_Leaf_Flow
            Sunlit_Boundary_Layer_Flow = (1 - np.exp(-(0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coeff) * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coeff) * Boundary_Layer_Leaf_Flow
            Shaded_Boundary_Layer_Flow = Canopy_Boundary_Layer_Flow - Sunlit_Boundary_Layer_Flow
            Boundary_Layer_Heat_Resist_Shaded = 1. / Shaded_Boundary_Layer_Flow  # Boundary layer resistance to heat, shaded part
            Boundary_Layer_Water_Resist_Shaded = 0.93 * Boundary_Layer_Heat_Resist_Shaded  # Boundary layer resistance to water, shaded part   
            
            # Potential conductance for CO2 for shaded leaves
            Conductance_CO2_Shaded = (Potential_Photosynthesis_Shaded - Dark_Respiration_Shaded) * (273. + Hourly_Temp) / 0.53717 / (self.leaf_object.Ambient_CO2 - Internal_CO2)
            
            # Turbulence resistance for canopy
            Turbulence_Resist_Canopy = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
                            
            Stomatal_Resist_Water_Shaded = max(1E-30, 1 / Conductance_CO2_Shaded - Boundary_Layer_Water_Resist_Shaded * 1.3 - Turbulence_Resist_Canopy) / 1.6
            
            Potential_Transpiration_Shaded, Net_Radiation_Absorbed_Shaded = Leaf.Penman_Monteith(Stomatal_Resist_Water_Shaded, Turbulence_Resist_Canopy, Boundary_Layer_Water_Resist_Shaded, Boundary_Layer_Heat_Resist_Shaded, Absorbed_Total_Radiation_Shaded, Atmospheric_Transmissivity, Fraction_Shaded, Hourly_Temp, Vapour_Pressure, Slope_SVP, Vapour_Pressure_Deficit)
            # Calculate leaf temperature for shaded leaves
            Latent_Heat_Vaporization = 2.4E6  # Latent heat of water vaporization (J/kg)
            Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity (J/m^3/°C)
            
            Temperature_Difference_Shaded = (Net_Radiation_Absorbed_Shaded - Latent_Heat_Vaporization * Potential_Transpiration_Shaded) * (Boundary_Layer_Heat_Resist_Shaded + (Turbulence_Resist_Canopy * Fraction_Shaded)) / Volumetric_Heat_Capacity_Air
            hourly_Air_Temperature_Difference_Shaded.append(Temperature_Difference_Shaded)
            
            Shaded_Leaf_Temperature = Hourly_Temp + Temperature_Difference_Shaded
            hourly_Shaded_Leaf_Temperature.append(Shaded_Leaf_Temperature)

        self.leaf_object.hourly_SH_leaf_T=hourly_Shaded_Leaf_Temperature
        self.leaf_object.hourly_Air_SH_leaf_T_diff=hourly_Air_Temperature_Difference_Shaded


    def Calculate_Potential_Photosynthesis(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed):
        # Use updated attributes for the calculations
        Total_LAI = self.leaf_object.leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.leaf_object.leaf_area_output['Final_LAI']
        Leaf_N_Extinction_Coeff = self.leaf_object.leaf_area_output['Leaf_N_Ext_Coeff']
        
        Specific_Leaf_N_Top = self.leaf_object.specific_leaf_n_output['Specific_Leaf_N_Top']
        
        # Convert daily conditions to hourly data for detailed analysis
        hourly_conditions, w_Gauss = self.leaf_object.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Daily_Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)
        hourly_Shaded_Leaf_Temp = self.leaf_object.hourly_Shaded_Leaf_Temp

        # Initialize lists to store outputs for hourly calculations
        hourly_Shaded_Leaf_Photosynthesis = []
        hourly_Shaded_Dark_Respiration = []
        for (Solar_Constant, Daytime_Temperature, Sin_Beam, Solar_Radiation, Wind_Speed), Shaded_Leaf_Temp in zip(hourly_conditions, hourly_Shaded_Leaf_Temp):

            # Use Shaded_Leaf_Temp where Daytime_Temperature was used
            Incoming_PAR = 0.5 * Solar_Radiation
            
            # Calculate the diffuse light fraction based on atmospheric transmissivity
            Atmospheric_Transmissivity = Solar_Radiation / (0.5 * Solar_Constant * Sin_Beam)
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity
            
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
        
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR

            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180 
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            
            Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
        
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
        
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
        
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057  # Canopy diffuse PAR reflection coefficient
        
            # Adjusting for vapor pressure deficit influence on intercellular CO2 concentration
            Vapor_Pressure_Deficit_Response = 0.116214  # Slope for linear effect of VPD on Ci/Ca
            Sat_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Shaded_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
        
            # Calculating photosynthetic nitrogen availability for shaded canopy parts
            Photosynthetic_Nitrogen_Shaded = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) - self.leaf_object.Minimum_Specific_Leaf_Nitrogen * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
        
            # Absorption of PAR by shaded leaves
            Absorbed_PAR_Shaded, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
        
            # Calculating potential photosynthesis and dark respiration for shaded leaves
            Potential_Photosynthesis_Shaded, Dark_Respiration_Shaded = Leaf.PHOTOSYN(self.leaf_object.C3C4_Pathway, Absorbed_PAR_Shaded, Shaded_Leaf_Temp, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Shaded, self.leaf_object.Activation_Energy_for_Jmax)
            
            
            hourly_Shaded_Leaf_Photosynthesis.append(Potential_Photosynthesis_Shaded)
            hourly_Shaded_Dark_Respiration.append(Dark_Respiration_Shaded)
            
            self.leaf_object.Hourly_Photosynthesis_Shaded = hourly_Shaded_Leaf_Photosynthesis
            self.leaf_object.Hourly_Dark_Respiration_Shaded = hourly_Shaded_Dark_Respiration
            
        # Aggregate hourly data back to daily totals
        Daily_Potential_Photosynthesis_Shaded = Leaf.Aggregate_To_Daily(hourly_Shaded_Leaf_Photosynthesis, Day_Length)
        Daily_Dark_Respiration_Shaded = Leaf.Aggregate_To_Daily(hourly_Shaded_Dark_Respiration, Day_Length)
        
        self.leaf_object.Daily_Potential_Photosynthesis_Shaded = Daily_Potential_Photosynthesis_Shaded
        self.leaf_object.Daily_Dark_Respiration_Shaded = Daily_Dark_Respiration_Shaded
              
    def Calculate_Potential_Transpiration(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height):
        # Use updated attributes
        Total_LAI = self.leaf_object.leaf_area_output['Total_LAI_ini']
        Leaf_Area_Index = self.leaf_object.leaf_area_output['Final_LAI']
        Wind_Extinction_Coeff = self.leaf_object.leaf_area_output['Ext_Coeff_Wind']
    
    
        hourly_conditions, w_Gauss = self.leaf_object.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)

    
        hourly_Shaded_Leaf_Temp = self.leaf_object.hourly_Shaded_Leaf_Temp
        hourly_Photosynthesis_Shaded = self.leaf_object.hourly_photosyn_SH
        hourly_Dark_Respiration_Shaded = self.leaf_object.hourly_darkresp_SH
    
        hourly_Transpiration_Shaded = []
        hourly_Absorbed_Radiation_Shaded = []
        hourly_Stomatal_Resist_Water_Shaded = []
        hourly_SVP_Slope_Shaded = []
        hourly_Absorbed_PAR_Shaded = []         
            
        
        for (_, Hourly_Temperature, Sin_Beam, Diffuse_Ratio, Wind_Speed), Shaded_Leaf_Temp, Photosynthesis_Shaded, Dark_Respiration_Shaded in zip(hourly_conditions, hourly_Shaded_Leaf_Temp, hourly_Photosynthesis_Shaded, hourly_Dark_Respiration_Shaded):
            
            Incoming_PAR = 0.5 * Solar_Radiation
            # Calculation of diffuse light fraction
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Solar_Declination)
            
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity
            
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Solar_Declination)))
            
            # Calculating the diffuse and direct components of PAR and NIR
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR
            
            Diffuse_NIR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_NIR = Incoming_PAR - Diffuse_NIR

            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Solar_Declination, Leaf_Blade_Angle_Radians)
            
            Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
            Scattering_Coefficient_NIR = 0.8  # Leaf scattering coefficient for NIR
            
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
            
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            Scattered_Beam_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient)
            
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057  # Canopy diffuse PAR reflection coefficient
            Canopy_Diffuse_Reflection_Coefficient_NIR = 0.389  # Canopy diffuse NIR reflection coefficient
            
            # Calculating the absorption of PAR and NIR by shaded leaves
            Absorbed_PAR_Shaded, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            Absorbed_NIR_Shaded, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_NIR, Diffuse_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR, Canopy_Diffuse_Reflection_Coefficient_NIR, Direct_NIR, Diffuse_NIR, Leaf_Area_Index)

        
            Total_Absorbed_Radiation_Shaded = Absorbed_PAR_Shaded + Absorbed_NIR_Shaded
    
            Fraction_Sunlit = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            Fraction_Shaded = 1 - Fraction_Sunlit
            
            # Calculating boundary layer resistance for shaded leaves
            Boundary_Layer_Resistance_Heat_Canopy_Factor = 0.01 * np.sqrt(Wind_Speed / self.leaf_object.Leaf_Width)
            Canopy_Boundary_Layer_Resistance_Heat = (1 - np.exp(-0.5 * Wind_Extinction_Coeff * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff) * Boundary_Layer_Resistance_Heat_Canopy_Factor
            Sunlit_Leaf_Boundary_Layer_Resistance_Heat = (1 - np.exp(-(0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Resistance_Heat_Canopy_Factor
            Shaded_Leaf_Boundary_Layer_Resistance_Heat = Canopy_Boundary_Layer_Resistance_Heat - Sunlit_Leaf_Boundary_Layer_Resistance_Heat

            Shaded_Leaf_Boundary_Layer_Resistance_Heat = 1. / Shaded_Leaf_Boundary_Layer_Resistance_Heat
            Shaded_Leaf_Boundary_Layer_Resistance_Water = 0.93 * Shaded_Leaf_Boundary_Layer_Resistance_Heat
            
            Vapor_Pressure_Deficit_Response =  0.116214 # Slope for linear effect of VPD on Ci/Ca

            Saturation_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Hourly_Temperature, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            Adjusted_Saturation_Vapor_Pressure, Adjusted_Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Shaded_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            
            Vapor_Pressure_Deficit = max(0, Saturation_Vapor_Pressure - Vapour_Pressure)
            
            Slope_VPD_Response = (Adjusted_Saturation_Vapor_Pressure - Saturation_Vapor_Pressure) / self.Avoid_Zero_Division(Shaded_Leaf_Temp- Hourly_Temperature)
               
            
            
            # Adjusted CO2 conductance for transpiration
            Adjusted_CO2_Conductance = (Photosynthesis_Shaded - Dark_Respiration_Shaded) * (273.15 + Shaded_Leaf_Temp) / 0.53717 / (self.leaf_object.Ambient_CO2 - Adjusted_Intercellular_CO2_Concentration)
            
            # Turbulence resistance for canopy
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            
            # Adjusted stomatal resistance to water vapor
            Adjusted_Stomatal_Resistance = max(1E-30, 1 / Adjusted_CO2_Conductance - Shaded_Leaf_Boundary_Layer_Resistance_Water * 1.3 - Turbulence_Resistance * Fraction_Shaded) / 1.6
            
            # Calculating potential transpiration for shaded leaves and absorbed radiation
            Potential_Transpiration_Shaded, Absorbed_Radiation_Shaded = Leaf.Penman_Monteith(Adjusted_Stomatal_Resistance, Turbulence_Resistance * Fraction_Shaded, Shaded_Leaf_Boundary_Layer_Resistance_Water, Shaded_Leaf_Boundary_Layer_Resistance_Heat, Total_Absorbed_Radiation_Shaded, Atmospheric_Transmissivity, Fraction_Shaded, Shaded_Leaf_Temp, Vapour_Pressure, Slope_VPD_Response, Vapor_Pressure_Deficit)
            
            # Appending calculated values to their respective lists
            hourly_Transpiration_Shaded.append(Potential_Transpiration_Shaded)
            hourly_Absorbed_Radiation_Shaded.append(Absorbed_Radiation_Shaded)
            hourly_Stomatal_Resist_Water_Shaded.append(Adjusted_Stomatal_Resistance)
            hourly_SVP_Slope_Shaded.append(Slope_VPD_Response)
            hourly_Absorbed_PAR_Shaded.append(Absorbed_PAR_Shaded)
            
        self.leaf_object.hourly_transpiration_SH = hourly_Transpiration_Shaded
        self.leaf_object.hourly_Aradiation_SH = hourly_Absorbed_Radiation_Shaded
        self.leaf_object.hourly_rsw_SH = hourly_Stomatal_Resist_Water_Shaded
        self.leaf_object.hourly_slope_SH = hourly_SVP_Slope_Shaded
        self.leaf_object.hourly_apar_SH = hourly_Absorbed_PAR_Shaded


        
        potential_transpiration_SH_daily = Leaf.aggregate_to_daily(hourly_Transpiration_Shaded, Day_Length)
        Aradiation_SH_daily = Leaf.aggregate_to_daily(hourly_Absorbed_Radiation_Shaded, Day_Length)
        self.leaf_object.potential_transpiration_SH_daily = potential_transpiration_SH_daily
        self.leaf_object.Aradiation_SH_daily = Aradiation_SH_daily
                            
    def Update_LeafTemp_Photosynthesis_if_WaterStress(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height,
                                                      Daily_Water_Supply, Soil_Depth_1, Root_Depth,
                                                      hourly_Sunlit_Leaf_Temp, hourly_Shaded_Leaf_Temp,
                                                      hourly_Soil_Evap):
    
        Gauss_Points = 5
        Gauss_Weights = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
        # Now use updated attributes
        Total_LAI = self.leaf_object.leaf_area_output['Total_LAI_ini']
        Leaf_Area_Index = self.leaf_object.leaf_area_output['Final_LAI']
        Leaf_N_Extinction_Coeff = self.leaf_object.leaf_area_output['Leaf_N_Ext_Coeff']
        Wind_Extinction_Coeff = self.leaf_object.leaf_area_output['Ext_Coeff_Wind']
        Specific_Leaf_N_Top = self.leaf_object.specific_leaf_n_output['Specific_Leaf_N_Top']
    
        # Convert daily data to hourly
        hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)
                           
        Actual_Photosynthesis = []
        Actual_Transpiration = []
        Actual_Air_Shaded_Leaf_Temp_Difference = []
        Actual_Shaded_Leaf_Temperature = []

        for i, (Solar_Constant, Hourly_Temp, Sin_Beam, Diffuse_Ratio, Wind_Speed), Potential_Transpiration_Sunlit, Potential_Transpiration_Shaded, Radiation_Shaded, Stomatal_Resist_Water_Shaded, SVP_Slope_Shaded, Soil_Evaporation in zip(range(Gauss_Points), hourly_data,
                                                                                                            self.leaf_object.hourly_transpiration_Sunlit,
                                                                                                            self.leaf_object.hourly_transpiration_Shaded,
                                                                                                            self.leaf_object.hourly_Aradiation_Shaded,
                                                                                                            self.leaf_object.hourly_rsw_Shaded,
                                                                                                            self.leaf_object.hourly_slope_Shaded,
                                                                                                            hourly_Soil_Evap):
        
            Hour = 12 - 0.5 * Day_Length + Day_Length * Gauss_Weights[i]
            Sin_Beam_Adjusted = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (Hour - 12.) / 24.))
            
            # Diurnal availability of soil water supply 
            Water_Supply_Hourly = Daily_Water_Supply * (Sin_Beam_Adjusted * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
            
            # Available water for evaporation (from the top soil layer)
            Water_Supply_Evaporation = Water_Supply_Hourly * (Soil_Depth_1 / Root_Depth)
            
            # Sum of sunlit and shaded hourly (or sub-daily) transpiration
            Potential_Canopy_Transpiration = Potential_Transpiration_Sunlit + Potential_Transpiration_Shaded
            
            # Amount of transpiration from the top soil layer
            Transpiration_From_Top_Soil_Layer = Potential_Canopy_Transpiration * (Soil_Depth_1 / Root_Depth)
            
            # Maximum possible transpiration
            Max_Possible_Transpiration = Transpiration_From_Top_Soil_Layer / (Transpiration_From_Top_Soil_Layer + Soil_Evaporation) * Water_Supply_Evaporation + (Water_Supply_Hourly - Water_Supply_Evaporation)
            
            if Max_Possible_Transpiration < Potential_Canopy_Transpiration:
                Actual_Canopy_Transpiration = Max_Possible_Transpiration
            else:
                Actual_Canopy_Transpiration = Potential_Canopy_Transpiration


            Potential_Transpiration_Canopy = Potential_Transpiration_Sunlit + Potential_Transpiration_Shaded
            Actual_Transpiration_Shaded = Potential_Transpiration_Shaded * (Actual_Canopy_Transpiration / Potential_Transpiration_Canopy)   # Actual transpiration of shaded leaves mm s-1


            Incoming_PAR = 0.5 * Solar_Radiation
            
            # Calculation of diffuse light fraction
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity
            
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
            
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR

            Vapor_Pressure_Deficit_Response =  0.116214  # Slope for linear effect of VPD on Ci/Ca
            
            Leaf_Blade_Angle_Radians = self.leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Solar_Declination, Leaf_Blade_Angle_Radians)
            
            Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
            
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Leaf_Area_Index, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            
            Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057
            
            Fraction_Sunlit = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            Fraction_Shaded = 1 - Fraction_Sunlit
            
            
            # Calculating boundary layer resistance for shaded leaves
            Boundary_Layer_Resistance_Heat_Canopy_Factor = 0.01 * np.sqrt(Wind_Speed / self.leaf_object.Leaf_Width)
            Canopy_Boundary_Layer_Resistance_Heat = (1 - np.exp(-0.5 * Wind_Extinction_Coeff * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff) * Boundary_Layer_Resistance_Heat_Canopy_Factor
            Sunlit_Leaf_Boundary_Layer_Resistance_Heat = (1 - np.exp(-(0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Resistance_Heat_Canopy_Factor
            Shaded_Leaf_Boundary_Layer_Resistance_Heat = Canopy_Boundary_Layer_Resistance_Heat - Sunlit_Leaf_Boundary_Layer_Resistance_Heat
           
            Shaded_Leaf_Boundary_Layer_Resistance_Heat = 1. / Shaded_Leaf_Boundary_Layer_Resistance_Heat
            Shaded_Leaf_Boundary_Layer_Resistance_Water = 0.93 * Shaded_Leaf_Boundary_Layer_Resistance_Heat
            
            # Turbulence resistance for canopy
            Turbulence_Resistance_Canopy = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            
            # Leaf temperature if water stress occurs
            Latent_Heat_of_Vaporization = 2.4E6  # Latent heat of water vaporization (J/kg)
            Volumetric_Heat_Capacity = 1200  # Volumetric heat capacity (J/m^3/°C)
            
            RT_Canopy = Turbulence_Resistance_Canopy * Fraction_Shaded
            Temperature_Difference = Leaf.Limit_Function(-25., 25., (Radiation_Shaded - Latent_Heat_of_Vaporization * Actual_Transpiration_Shaded) * (Shaded_Leaf_Boundary_Layer_Resistance_Heat + RT_Canopy) / Volumetric_Heat_Capacity)
            
            Leaf_Temperature = Hourly_Temp + Temperature_Difference
            
            Psychrometric_Constant = 0.067  # psychrometric constant (kPa/oC)


            # Stomatal resistance to water if water stress occurs
            Adjusted_Stomatal_Resistance_Water_Stress = (Potential_Transpiration_Shaded - self.leaf_object.hourly_Aradiation_Shaded) * (self.leaf_object.hourly_slope_Shaded * (Shaded_Leaf_Boundary_Layer_Resistance_Heat + Turbulence_Resistance_Canopy) + Psychrometric_Constant * (Shaded_Leaf_Boundary_Layer_Resistance_Water + Turbulence_Resistance_Canopy)) / self.leaf_object.hourly_Aradiation_Shaded / Psychrometric_Constant + Potential_Transpiration_Shaded / self.leaf_object.hourly_Aradiation_Shaded * self.leaf_object.hourly_rsw_Shaded
            
            # Absorbed PAR by shaded leaves
            _, Absorbed_PAR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            
            # Total photosynthetic nitrogen in the canopy
            Photosynthetic_Nitrogen_Canopy = (Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) - self.leaf_object.Minimum_Specific_Leaf_Nitrogen * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient)
            

            # Photosynthetic nitrogen for shaded parts of the canopy
            Photosynthetic_Nitrogen_Shaded = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) - self.leaf_object.Minimum_Specific_Leaf_Nitrogen * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            Photosynthetic_Nitrogen_Sunlit = Photosynthetic_Nitrogen_Canopy - Photosynthetic_Nitrogen_Shaded
            
            # Calculating internal CO2 concentration and saturation vapor pressure
            Saturation_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(self.leaf_object.hourly_Shaded_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.leaf_object.Ambient_CO2, self.leaf_object.C3C4_Pathway)
            Absorbed_PAR_Leaf, Dark_Respiration = Leaf.PHOTOSYN(self.leaf_object.C3C4_Pathway, Absorbed_PAR_Shaded, self.leaf_object.hourly_Shaded_Leaf_Temp, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Shaded, self.leaf_object.Activation_Energy_for_Jmax, self.leaf_object.Vcmax_Nitrogen_Response, self.leaf_object.Jmax_Nitrogen_Response, self.leaf_object.Photosynthetic_Thermal_Tolerance)
            
            # Actual photosynthesis under water stress condition
            Actual_Photosynthesis_Water_Stress = (1.6 * self.leaf_object.hourly_rsw_Shaded + 1.3 * Shaded_Leaf_Boundary_Layer_Resistance_Water + Turbulence_Resistance_Canopy) / (1.6 * Adjusted_Stomatal_Resistance_Water_Stress + 1.3 * Shaded_Leaf_Boundary_Layer_Resistance_Water + Turbulence_Resistance_Canopy) * (Absorbed_PAR_Leaf - Dark_Respiration) + Dark_Respiration
            
            # Appending the calculated values to their respective lists
            Actual_Photosynthesis.append(Actual_Photosynthesis_Water_Stress)
            Actual_Transpiration.append(self.leaf_object.hourly_Aradiation_Shaded)
            Actual_Air_Shaded_Leaf_Temp_Difference.append(Temperature_Difference)
            Actual_Shaded_Leaf_Temperature.append(self.leaf_object.hourly_Shaded_Leaf_Temp)
            
        self.leaf_object.Hourly_Actual_Photosynthesis_Shaded = Actual_Photosynthesis
        self.leaf_object.Hourly_Actual_Transpiration_Shaded = Actual_Transpiration
        self.leaf_object.Hourly_Actual_Air_Shaded_Leaf_Temp_Diff = Actual_Air_Shaded_Leaf_Temp_Difference
        self.leaf_object.Hourly_Actual_Shaded_Leaf_Temp = Actual_Shaded_Leaf_Temperature
        










