import numpy as np
import math
import Canopy

def SWITCH_FUN(x, y1, y2):
    return y1 if x < 0 else y2

# Constants
Stefan_Boltzmann_Const = 5.668E-8  # Stefan-Boltzmann constant (J/m2/s/K4)
Latent_Heat_Vaporization = 2.4E6  # Latent heat of water vaporization (J/kg)
Volumetric_Heat_Capacity_Air = 1200  # Volumetric heat capacity (J/m3/°C)
Psychrometric_Constant = 0.067  # Psychrometric constant (kPa/°C)

# Activation energies and other constants related to the Farquhar model
O2_Concentration = 210  # Oxygen concentration (mmol/mol)
Activation_Energy_VCMAX = 65330  # Energy of activation for VCMAX (J/mol)
Activation_Energy_KMC = 79430  # Energy of activation for KMC (J/mol)
Activation_Energy_KMO = 36380  # Energy of activation for KMO (J/mol)
Activation_Energy_Dark_Respiration = 46390  # Energy of activation for dark respiration (J/mol)
Deactivation_Energy_Jmax = 200000  # Energy of deactivation for JMAX (J/mol)
Entropy_Term_JT_Equation = 650  # Entropy term in JT equation (J/mol/K)
Maximum_Electron_Transport_Efficiency = 0.85  # Maximum electron transport efficiency of PS II
Protons_For_ATP_Synthesis = 3  # Number of protons required to synthesize 1 ATP
Dark_Respiration_VCMAX_Ratio_25C = 0.0089  



Scattering_Coefficient_PAR = 0.2  # Leaf scattering coefficient for PAR
Scattering_Coefficient_NIR = 0.8  # Leaf scattering coefficient for NIR
Canopy_Diffuse_Reflection_Coefficient_PAR = 0.057  # Canopy diffuse PAR reflection coefficient
Canopy_Diffuse_Reflection_Coefficient_NIR = 0.389  # Canopy diffuse NIR reflection coefficient


##self.Plant_Height ... take care of Plant_Height because it's should be comming from canopy
class Leaf: 
    def __init__(self,  SLA_Const,Min_Specific_Leaf_N, Leaf_Blade_Angle, Leaf_Width, C3C4_Pathway, Ambient_CO2, Activation_Energy_Jmax,
                 Vcmax_LeafN_Slope, Jmax_LeafN_Slope, Photosynthetic_Light_Response_Factor ):
        self.Specific_Leaf_Area = SLA_Const
        self.Leaf_Blade_Angle = Leaf_Blade_Angle
        self.Leaf_Width = Leaf_Width
        self.Min_Specific_Leaf_N = Min_Specific_Leaf_N
        self.C3C4_Pathway = C3C4_Pathway
        self.Ambient_CO2 = Ambient_CO2
        self.Activation_Energy_Jmax = Activation_Energy_Jmax
        self.Vcmax_LeafN_Slope = Vcmax_LeafN_Slope
        self.Jmax_LeafN_Slope = Jmax_LeafN_Slope
        self.Photosynthetic_Light_Response_Factor = Photosynthetic_Light_Response_Factor

        # Initial conditions
        self.Carbon_Dead_Leaves=self.Carbon_Dead_Leaves_Litters_Soil = self.Dif_Air_Leaf_Temp=0         

        # Initialize attributes to store outputs
        self.Leaf_area_output = {}
        self.specific_Leaf_N_output = {}
               
        
        self.Hourly_Air_Sunlit_Leaf_Temp_diff=[]
        self.Hourly_Air_Shaded_Leaf_Temp_diff=[]

        
        self.Hourly_Sunlit_Leaf_Temp=[]
        self.Hourly_Shaded_Leaf_Temp=[]
        
        self.Hourly_Photosynthesis_Sunlit=[]
        self.Hourly_Dark_Respiration_Sunlit=[]
        self.Hourly_Transpiration_Sunlit=[]
        self.Hourly_Absorbed_Radiation_Sunlit=[]
        self.Hourly_Stomatal_Resistance_Water_Sunlit=[]
        self.Hourly_Slope_VPD_Sunlit=[]
        self.Hourly_Absorbed_PAR_Sunlit=[]


        self.Hourly_Actual_Photosynthesis_Sunlit=[]
        self.Hourly_Actual_Photosynthesis_Sunlit_DELTA=[]
        # self.Hourly_Actual_Transpiration_Sunlit=[]
        self.Hourly_Actual_Air_Sunlit_Leaf_Temp_Diff=[]
        self.Hourly_Actual_Sunlit_Leaf_Temp=[]
        
        
        self.Hourly_Photosynthesis_Shaded=[]
        self.Hourly_Actual_Photosynthesis_Shaded_DELTA=[]
        self.Hourly_Dark_Respiration_Shaded=[]
        self.Hourly_Transpiration_Shaded=[]
        self.Hourly_Absorbed_Radiation_Shaded=[]
        self.Hourly_Stomatal_Resistance_Water_Shaded=[]
        self.Hourly_Slope_VPD_Shaded=[]
        self.Hourly_Absorbed_PAR_Shaded=[]

        
        self.Hourly_Actual_Photosynthesis_Shaded=[]
        # self.Hourly_Actual_Transpiration_Shaded=[]
        self.Hourly_Actual_Air_Shaded_Leaf_Temp_Diff=[]
        self.Hourly_Actual_Shaded_Leaf_Temp=[]

        # self.potential_photosyn_Sunlit_daily=0
        # self.dark_rsp_Sunlit_daily=0
        # self.potential_photosyn_Shaded_daily=0
        # self.dark_rsp_Shaded_daily=0
        # self.potential_transpiration_Sunlit_daily=0
        # self.potential_transpiration_Shaded_daily=0
        # self.Apar_Shaded_daily=0  
        # self.Apar_Sunlit_daily=0  


    
    

    
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


    def Update_Leaf_Area(self,Tot_Leaf_N,Leaf_N,CarbonFrac_Veg,Carbon_determined_LAI):

        # Calculating delta and total Leaf area index (LAI)
        Dead_LAI = (self.Carbon_Dead_Leaves - self.Carbon_Dead_Leaves_Litters_Soil) / CarbonFrac_Veg * self.Specific_Leaf_Area
        Total_LAI = Carbon_determined_LAI + self.Carbon_Dead_Leaves / CarbonFrac_Veg * self.Specific_Leaf_Area

        # Nitrogen and light extinction coefficients
        Light_Ext_Coeff = Leaf.KDF_Coeff(Total_LAI, self.Leaf_Blade_Angle * np.pi / 180., 0.2)
        Nitro_Ext_Coeff = Light_Ext_Coeff * (Tot_Leaf_N - self.Min_Specific_Leaf_N * Total_LAI)
        intermediate_var = self.Min_Specific_Leaf_N * (1.0 - np.exp(-Light_Ext_Coeff * Total_LAI))
        
        # Assuming wind extinction coefficient is similar to light for simplicity
        Wind_Ext_Coeff = Light_Ext_Coeff  
        Leaf_Nitro_Ext_Coeff = 1.0 / Total_LAI * math.log((Nitro_Ext_Coeff + intermediate_var) / (Nitro_Ext_Coeff * math.exp(-Light_Ext_Coeff * Total_LAI) + intermediate_var))

        # Integrating LAI considering nitrogen effect
        N_determined_LAI = math.log(1. + Leaf_Nitro_Ext_Coeff * max(0., Leaf_N) / self.Min_Specific_Leaf_N) / Leaf_Nitro_Ext_Coeff
        Leaf_Area_Index = min(N_determined_LAI, Carbon_determined_LAI)
        # print(N_determined_LAI, Carbon_determined_LAI)
        # Updating Leaf area outputs
        self.Leaf_area_output = {
            'Dead_LAI': Dead_LAI,
            'Total_LAI': Total_LAI,
            'Light_Ext_Coeff': Light_Ext_Coeff,
            'Nitro_Ext_Coeff': Nitro_Ext_Coeff,
            'Wind_Ext_Coeff': Wind_Ext_Coeff,
            'Leaf_Nitro_Ext_Coeff': Leaf_Nitro_Ext_Coeff,
            'Leaf_Area_Index': Leaf_Area_Index,
            'N_determined_LAI': N_determined_LAI
        }
        # print(self.Leaf_area_output)


    def Update_Specific_Leaf_N(self,Leaf_N):
        # Extracting final LAI and Leaf nitrogen extinction coefficient
        Leaf_Area_Index = self.Leaf_area_output['Leaf_Area_Index']
        Leaf_Nitro_Ext_Coeff = self.Leaf_area_output['Leaf_Nitro_Ext_Coeff']

        # Calculating specific Leaf nitrogen
        Specific_Leaf_N = Leaf_N / Leaf_Area_Index  # Average specific Leaf nitrogen
        Specific_Leaf_N_Top = Leaf_N * Leaf_Nitro_Ext_Coeff / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Leaf_Area_Index))  # For top leaves
        # print(Tot_Leaf_N)
        
        # Specific_Leaf_N_Bottom_Exponential_with_Depth = Tot_Leaf_N * Leaf_Nitro_Ext_Coeff * np.exp(-Leaf_Nitro_Ext_Coeff * Leaf_Area_Index) / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Leaf_Area_Index))  # Exponential nitrogen profile
        Specific_Leaf_N_Bottom = Leaf_N * Leaf_Nitro_Ext_Coeff * np.exp(-Leaf_Nitro_Ext_Coeff * Leaf_Area_Index) / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Leaf_Area_Index))  # Exponential nitrogen profile

        Specific_Leaf_N_Top_Increment = (Leaf_N + 0.001 * Leaf_N) * Leaf_Nitro_Ext_Coeff / (1. - np.exp(-Leaf_Nitro_Ext_Coeff * Leaf_Area_Index))  # With small nitrogen increment

        # Updating specific Leaf nitrogen outputs
        self.specific_Leaf_n_output = {
            'Specific_Leaf_N': Specific_Leaf_N,
            'Specific_Leaf_N_Top': Specific_Leaf_N_Top,
            'Specific_Leaf_N_Top_Increment': Specific_Leaf_N_Top_Increment,
            'Specific_Leaf_N_Bottom': Specific_Leaf_N_Bottom,
        }

    def aggregate_to_daily(Hourly_results, Day_Length):
       
        wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])

        var_daily = 0
        
        for i, var_hourly in enumerate(Hourly_results):
            var_daily += var_hourly * wgauss[i]
            # print(var_daily)
        # Convert to daily totals
        var_daily *= Day_Length * 3600
    
        return var_daily
    
    def convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind):

         Gaussian_Points = 5
         x_Gauss = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
         w_Gauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
         
         Instantaneous_data = []
         for i in range(Gaussian_Points):
             hour = 12 - 0.5 * Day_Length + Day_Length * x_Gauss[i]
             Solar_Elevation_Sin = max(0., Sin_Solar_Declination + Cos_Solar_Declination * np.cos(2. * np.pi * (hour - 12.) / 24.))
             Instantaneous_radiation = Solar_Radiation * (Solar_Elevation_Sin * Solar_Constant/1367) / Daily_Sin_Beam_Exposure
             Instantaneous_Temp = Min_Temp + (Max_Temp - Min_Temp) * np.sin(np.pi * (hour + Day_Length / 2 - 12) / (Day_Length + 3))
             Instantaneous_data.append((Solar_Constant, Instantaneous_Temp, Solar_Elevation_Sin, Instantaneous_radiation, Wind))

         return Instantaneous_data, w_Gauss
    
    def INTERNAL_CO2(Leaf_Temp, VPD, VPD_Slope, Ambient_CO2, C3C4_Pathway):

        # Air-to-Leaf vapor pressure deficit
        Saturated_Vapor_Pressure_Leaf = 0.611 * np.exp(17.4 * Leaf_Temp / (Leaf_Temp + 239.))
        Vapor_Pressure_Deficit_Leaf = max(0, Saturated_Vapor_Pressure_Leaf - VPD)
        
        # Constants based on crop type
        Michaelis_Menten_CO2_25C = 650 if C3C4_Pathway == -1 else 404.9
        Michaelis_Menten_O2_25C = 450 if C3C4_Pathway == -1 else 278.4
        
        # Adjustment for temperature
        KMC = Michaelis_Menten_CO2_25C * np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 79430 / 8.314)
        KMO = Michaelis_Menten_O2_25C * np.exp((1./298. - 1./(Leaf_Temp + 273.)) * 36380 / 8.314)
        
        # CO2 compensation point without dark respiration
        CO2_compensation_point_no_resp = 0.5 * np.exp(-3.3801 + 5220./(Leaf_Temp + 273.) / 8.314) * 210 * KMC / KMO
        dark_respiration_Vcmax_ratio  = Dark_Respiration_VCMAX_Ratio_25C * np.exp((1/298 - 1/(Leaf_Temp + 273)) * (46390 - 65330) / 8.314)
        CO2_compensation_point_conditional =(CO2_compensation_point_no_resp + dark_respiration_Vcmax_ratio * KMC * (1 + 210 / KMO)) / (1 - dark_respiration_Vcmax_ratio) 
        CO2_compensation_point =  CO2_compensation_point_conditional/10 if C3C4_Pathway == -1 else CO2_compensation_point_conditional
        
        Intercellular_CO2_Ratio = 1 - (1 - CO2_compensation_point / Ambient_CO2) * (0.14 + VPD_Slope * Vapor_Pressure_Deficit_Leaf)
        
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
    
    def PHOTOSYN(C3C4_Pathway, Absorbed_PAR, Leaf_Temp, Intercellular_CO2, Photosynthetic_Active_Nitrogen, Activation_Energy_Jmax, Vcmax_LeafN_Slope
                 , Jmax_LeafN_Slope,    Photosynthetic_Light_Response_Factor):
       
        # Constants and initial calculations for both C3 and C4 plants, adapted for temperature
        Michaelis_Menten_CO2_25C = 650 if C3C4_Pathway == -1 else 404.9  # Michaelis-Menten constant for CO2 at 25°C, adjusted for C3/C4
        Michaelis_Menten_O2_25C = 450 if C3C4_Pathway == -1 else 278.4  # Michaelis-Menten constant for O2 at 25°C, adjusted for C3/C4
        

        
        # Conversion of absorbed PAR to photon flux density
        Photon_Flux_Density = 4.56 * Absorbed_PAR  # Conversion factor to umol/m^2/s
    
        # Adjusting Michaelis-Menten constants for CO2 and O2 with temperature
        CO2_Michaelis_Menten_Temp_Adjusted = Michaelis_Menten_CO2_25C * math.exp((1./298. - 1./(Leaf_Temp + 273.)) * Activation_Energy_KMC / 8.314)
        O2_Michaelis_Menten_Temp_Adjusted = Michaelis_Menten_O2_25C * math.exp((1./298. - 1./(Leaf_Temp + 273.)) * Activation_Energy_KMO / 8.314)
    
        # CO2 compensation point without dark respiration
        CO2_Compensation_No_Respiration = 0.5 * math.exp(-3.3801 + 5220. / (Leaf_Temp + 273.) / 8.314) * O2_Concentration * CO2_Michaelis_Menten_Temp_Adjusted / O2_Michaelis_Menten_Temp_Adjusted
    
        # Temperature effects on carboxylation and electron transport
        Carboxylation_Temperature_Effect = math.exp((1./298. - 1./(Leaf_Temp + 273.)) * Activation_Energy_VCMAX / 8.314)
        
        Electron_Transport_Temperature_Effect = math.exp((1./298. - 1./(Leaf_Temp + 273.)) * Activation_Energy_Jmax / 8.314) * \
            (1. + math.exp(Entropy_Term_JT_Equation / 8.314 - Deactivation_Energy_Jmax / 298. / 8.314)) / \
            (1. + math.exp(Entropy_Term_JT_Equation / 8.314 - 1. / (Leaf_Temp + 273.) * Deactivation_Energy_Jmax / 8.314))
    
        # Adjusted VCMAX and JMAX based on Leaf nitrogen content
        Adjusted_VCMAX = Vcmax_LeafN_Slope * Carboxylation_Temperature_Effect * Photosynthetic_Active_Nitrogen
        Adjusted_JMAX = Jmax_LeafN_Slope * Electron_Transport_Temperature_Effect * Photosynthetic_Active_Nitrogen

        

         # Assumption for electron transport
        Pseudocyclic_Electron_Transport = 0  # Assuming no pseudocyclic electron transport

        if C3C4_Pathway == -1:
            CO2_Leakage_Factor = 0.2  # CO2 leakage from bundle-sheath to mesophyll in C4 plants
            Concentrated_CO2 = 10 * Intercellular_CO2  # Mimicking C4 CO2 concentrating mechanism
            Scarcity_Factor = 2 * (Concentrated_CO2 - CO2_Compensation_No_Respiration) / (1. - CO2_Leakage_Factor)
            Quantum_Yield_Factor = 1 - Pseudocyclic_Electron_Transport - 2 * (4 * Concentrated_CO2 + 8 * CO2_Compensation_No_Respiration) / Protons_For_ATP_Synthesis / (Scarcity_Factor + 3 * Concentrated_CO2 + 7 * CO2_Compensation_No_Respiration)
            Cyclic_Electron_Flow = Quantum_Yield_Factor
        else:
            Concentrated_CO2 = Intercellular_CO2  # Directly using intercellular CO2 concentration for C3 plants
            Scarcity_Factor = 0
            Quantum_Yield_Factor = 0
            Cyclic_Electron_Flow = 1 - (Pseudocyclic_Electron_Transport * Protons_For_ATP_Synthesis * (Scarcity_Factor + 3 * Concentrated_CO2 + 7 * CO2_Compensation_No_Respiration) / (4 * Concentrated_CO2 + 8 * CO2_Compensation_No_Respiration) + 1) / \
                                   (Protons_For_ATP_Synthesis * (Scarcity_Factor + 3 * Concentrated_CO2 + 7 * CO2_Compensation_No_Respiration) / (4 * Concentrated_CO2 + 8 * CO2_Compensation_No_Respiration) - 1)

        # Electron transport rate in response to absorbed PAR photon flux
        Quantum_Efficiency_Adjustment = (1 - Cyclic_Electron_Flow) / (1 + (1 - Cyclic_Electron_Flow) / Maximum_Electron_Transport_Efficiency)
        Electron_Transport_Ratio = Quantum_Efficiency_Adjustment * Photon_Flux_Density / max(1E-10, Adjusted_JMAX)
        Adjusted_Electron_Transport_Rate = Adjusted_JMAX * (1 + Electron_Transport_Ratio - ((1 + Electron_Transport_Ratio)**2 - 4 * Electron_Transport_Ratio * Photosynthetic_Light_Response_Factor)**0.5) / 2 / Photosynthetic_Light_Response_Factor
    
        # Carboxylation rates limited by Rubisco activity and electron transport
        Carboxylation_Rate_Rubisco_Limited = Adjusted_VCMAX * Concentrated_CO2 / (Concentrated_CO2 + CO2_Michaelis_Menten_Temp_Adjusted * (O2_Concentration / O2_Michaelis_Menten_Temp_Adjusted + 1.))
        Carboxylation_Rate_Electron_Transport_Limited = Adjusted_Electron_Transport_Rate * Concentrated_CO2 * (2 + Quantum_Yield_Factor - Cyclic_Electron_Flow) / Protons_For_ATP_Synthesis / (Scarcity_Factor + 3 * Concentrated_CO2 + 7 * CO2_Compensation_No_Respiration) / (1 - Cyclic_Electron_Flow)

        # Gross rate of Leaf photosynthesis
        Photosynthesis_Efficiency = (1 - CO2_Compensation_No_Respiration / Concentrated_CO2) * min(Carboxylation_Rate_Rubisco_Limited, Carboxylation_Rate_Electron_Transport_Limited)
        Gross_Leaf_Photosynthesis = max(1E-10, (1E-6) * 44 * Photosynthesis_Efficiency)
    
        # Rate of Leaf dark respiration
        
        Temporary_var = math.exp((1/298 - 1/(Leaf_Temp + 273)) * Activation_Energy_Dark_Respiration / 8.314)
        Leaf_Dark_Respiration = (1E-6) * 44 * Dark_Respiration_VCMAX_Ratio_25C * (Vcmax_LeafN_Slope * Photosynthetic_Active_Nitrogen) * Temporary_var
        
        
        
        return Gross_Leaf_Photosynthesis, Leaf_Dark_Respiration

    def Penman_Monteith(Stomatal_Resist_Water, Turbulence_Resist, Boundary_Layer_Resist_Water, 
                        Boundary_Layer_Resist_Heat, Absorbed_Global_Radiation, Atmospheric_Transmissivity, 
                        Fraction_Leaf_Classes, Leaf_Temp, Vapour_Pressure, Saturated_Vapour_Pressure_Slope, VPD):


    
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
    
     # Potential Leaf transpiration calculation
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




class Leaf_sunlit(Leaf):
    def __init__(self, Leaf_object):
        self.Leaf_object = Leaf_object
        
    def Calculate_Leaf_temp(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, 
                            Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height,C3C4_Pathway):
        # Using updated attributes
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Wind_Ext_Coeff = self.Leaf_object.Leaf_area_output['Wind_Ext_Coeff']
        Leaf_Nitro_Ext_Coeff = self.Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff']
        Specific_Leaf_N_Top = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top']
        #print(Leaf_Area_Index)
        Hourly_data, w_Gauss = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)        
        Hourly_Sunlit_Leaf_Temp = []
        Hourly_Air_Leaf_Temp_Diff = []
        # print(Hourly_data)

        for Solar_Constant, Hourly_Temp, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed in Hourly_data:
            
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
            Incoming_NIR = 0.5 * Hourly_Solar_Radiation

            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)

            # Calculation of Vapor Pressure Deficit effect
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            Sat_Vapor_Pressure, Intercellular_CO2 = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
    
            
            # Dynamic slope calculation for the response of saturation vapor pressure to temperature
            Vapour_Pressure_Deficit = max(0, Sat_Vapor_Pressure - Vapour_Pressure)
            Slope_SVP = 4158.6 * Sat_Vapor_Pressure / (Hourly_Temp + 239.)**2





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
            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            

            
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
            
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            Scattered_Beam_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient)
            

            
                        
            # Calculating photosynthetic nitrogen for sunlit canopy parts
            Photosynthetic_Nitrogen_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            # print(Specific_Leaf_N_Sunlit)

            # Absorption of PAR and NIR by sunlit leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR,
                                                       Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR,
                                                       Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            Absorbed_NIR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_NIR, Diffuse_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR, Canopy_Diffuse_Reflection_Coefficient_NIR, Direct_NIR, Diffuse_NIR, Leaf_Area_Index)
            # print("Absorbed_PAR_Sunlit")
            # print(Absorbed_PAR_Sunlit)
            
            
            # Calculating potential photosynthesis and dark respiration for sunlit leaves
            Potential_Photosynthesis_Sunlit, Dark_Respiration_Sunlit = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Hourly_Temp, Intercellular_CO2, Photosynthetic_Nitrogen_Sunlit, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)
            # print(Potential_Photosynthesis_Sunlit)

            # Total absorbed radiation (PAR + NIR) by sunlit leaves
            Total_Absorbed_Radiation_Sunlit = Absorbed_PAR_Sunlit + Absorbed_NIR_Sunlit

            # Calculating the fraction of sunlit and shaded canopy components
            Sunlit_Fraction = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            # print(Sunlit_Fraction)

            # Boundary layer resistance for heat and water vapor for sunlit leaves
            Boundary_Layer_Conductance_Heat=(0.01 * np.sqrt(Wind_Speed / self.Leaf_object.Leaf_Width))
            Boundary_Layer_Conductance_Heat_Sunlit = Boundary_Layer_Conductance_Heat*(1 - np.exp(-(0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient)
            Boundary_Layer_Resistance_Heat_Sunlit = 1. / Boundary_Layer_Conductance_Heat_Sunlit
            Boundary_Layer_Resistance_Water_Sunlit = 0.93 * Boundary_Layer_Resistance_Heat_Sunlit

            
            # print(Boundary_Layer_Resistance_Heat_Sunlit)
            # Calculating potential CO2 conductance for sunlit leaves
            Potential_CO2_Conductance_Sunlit = (Potential_Photosynthesis_Sunlit - Dark_Respiration_Sunlit) * (273.15 + Hourly_Temp) / 0.53717 / (self.Leaf_object.Ambient_CO2 - Intercellular_CO2)
            # print(Potential_CO2_Conductance_Sunlit)
            # Calculating Turbulence Resistance for the canopy
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            Turbulence_Resistance_Sunlit=Turbulence_Resistance * Sunlit_Fraction
            Stomatal_Resistance_Sunlit = max(1E-30, 1 / Potential_CO2_Conductance_Sunlit - Boundary_Layer_Resistance_Water_Sunlit * 1.3 - Turbulence_Resistance_Sunlit) / 1.6
            # print(Stomatal_Resistance_Sunlit)

            # Applying Penman-Monteith equation to calculate transpiration and net radiation absorbed
            Transpiration_Sunlit, Net_Radiation_Absorbed_Sunlit = Leaf.Penman_Monteith(Stomatal_Resistance_Sunlit,Turbulence_Resistance_Sunlit ,
                                                                                       Boundary_Layer_Resistance_Water_Sunlit, Boundary_Layer_Resistance_Heat_Sunlit,
                                                                                       Total_Absorbed_Radiation_Sunlit, Atmospheric_Transmissivity,
                                                                                       Sunlit_Fraction, Hourly_Temp, Vapour_Pressure, 
                                                                                       Slope_SVP, Vapour_Pressure_Deficit)
            # print(Net_Radiation_Absorbed_Sunlit)

            # Calculating Leaf temperature difference
            Air_Leaf_Temp_Diff_Sunlit = (Net_Radiation_Absorbed_Sunlit - Latent_Heat_Vaporization * Transpiration_Sunlit) * (Boundary_Layer_Resistance_Heat_Sunlit + Turbulence_Resistance_Sunlit) / Volumetric_Heat_Capacity_Air
            Air_Leaf_Temp_Diff_Sunlit=max(-25, min(Air_Leaf_Temp_Diff_Sunlit, 25)) 
            # print("leaf temp")
            # print(Air_Leaf_Temp_Diff_Sunlit)
            # Appending calculated temperature difference and adjusted Leaf temperature to lists
            Hourly_Air_Leaf_Temp_Diff.append(Air_Leaf_Temp_Diff_Sunlit)
            Adjusted_Leaf_Temp_Sunlit = Hourly_Temp + Air_Leaf_Temp_Diff_Sunlit
            # print(Boundary_Layer_Resistance_Heat_Sunlit + Turbulence_Resistance_Sunlit)

            Hourly_Sunlit_Leaf_Temp.append(Adjusted_Leaf_Temp_Sunlit)
        self.Leaf_object.Hourly_Sunlit_Leaf_Temp=Hourly_Sunlit_Leaf_Temp
        self.Leaf_object.Hourly_Air_Sunlit_Leaf_Temp_diff=Hourly_Air_Leaf_Temp_Diff             
        
        
    def Calculate_Potential_Photosynthesis(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, 
                                           Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed,C3C4_Pathway):
        # Use updated attributes
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Leaf_Nitro_Ext_Coeff = self.Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff']
        Specific_Leaf_N_Top = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top']

        # Convert daily data to hourly
        Hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)
        
        Hourly_Sunlit_Leaf_Temp = self.Leaf_object.Hourly_Sunlit_Leaf_Temp
        
        Hourly_Photosynthesis = []
        Hourly_Dark_Respiration = []
    
        for (_, Hourly_Temp, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed), Sunlit_Leaf_Temp in zip(Hourly_data, Hourly_Sunlit_Leaf_Temp):
            # Use sunlit Leaf temperature if it differs from ambient

    
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
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
            
            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
                        
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            
            
            # Adjusting for vapor pressure deficit influence on intercellular CO2 concentration
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            Sat_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Sunlit_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            
            # Calculating photosynthetic nitrogen availability for sunlit canopy parts
            Photosynthetic_Nitrogen_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_Nitro_Ext_Coeff + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient

            # Absorption of PAR by sunlit leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            
            
            # Calculating potential photosynthesis and dark respiration for sunlit leaves
            Potential_Photosynthesis_Sunlit, Dark_Respiration_Sunlit = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Sunlit_Leaf_Temp, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Sunlit, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)
            # print( Absorbed_PAR_Sunlit, Sunlit_Leaf_Temp, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Sunlit,)
                   
            # Appending the calculated photosynthesis and respiration rates to their respective hourly lists
            Hourly_Photosynthesis.append(Potential_Photosynthesis_Sunlit)
            Hourly_Dark_Respiration.append(Dark_Respiration_Sunlit)
            # print(Specific_Leaf_N_Top)
        # print('Next day')

        self.Leaf_object.Hourly_Photosynthesis_Sunlit = Hourly_Photosynthesis
        self.Leaf_object.Hourly_Dark_Respiration_Sunlit = Hourly_Dark_Respiration
        # print(self.Leaf_object.Hourly_Photosynthesis_Sunlit)

        # # Aggregate hourly data back to daily totals
        # Daily_Potential_Photosynthesis_Sunlit = Leaf.aggregate_to_daily(Hourly_Photosynthesis, Day_Length)
        # Daily_Dark_Respiration_Sunlit = Leaf.aggregate_to_daily(Hourly_Dark_Respiration, Day_Length)
        # #print(Daily_Potential_Photosynthesis_Sunlit)

        # self.Leaf_object.Daily_Potential_Photosynthesis_Sunlit = Daily_Potential_Photosynthesis_Sunlit
        # self.Leaf_object.Daily_Dark_Respiration_Sunlit = Daily_Dark_Respiration_Sunlit
                                    
                
        
    
    def Calculate_Potential_Transpiration(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, 
                                          Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height,C3C4_Pathway):
        # Utilizing updated attributes
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Wind_Ext_Coeff = self.Leaf_object.Leaf_area_output['Wind_Ext_Coeff']
    
        Hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)
    
        Hourly_Sunlit_Leaf_Temp = self.Leaf_object.Hourly_Sunlit_Leaf_Temp
        Hourly_Photosynthesis_Sunlit = self.Leaf_object.Hourly_Photosynthesis_Sunlit
        Hourly_Dark_Respiration_Sunlit = self.Leaf_object.Hourly_Dark_Respiration_Sunlit
    
        Hourly_Transpiration_Sunlit = []
        Hourly_Absorbed_Radiation_Sunlit = []
        Hourly_Stomatal_Resistance_Water_Sunlit = []
        Hourly_Slope_VPD_Sunlit = []
        Hourly_Absorbed_PAR_Sunlit = []
    
        for (_, Hourly_Temp, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed), Sunlit_Leaf_Temp, Photosynthesis_Sunlit, Dark_Respiration_Sunlit in zip(Hourly_data, Hourly_Sunlit_Leaf_Temp, Hourly_Photosynthesis_Sunlit, Hourly_Dark_Respiration_Sunlit):

            
            
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
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
            NIR = 0.5 * Hourly_Solar_Radiation
            
            # Calculation of Vapor Pressure Deficit effect
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            Sat_Vapor_Pressure, Intercellular_CO2 = Leaf.INTERNAL_CO2(Sunlit_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
    
            # Incoming Diffuse and Direct NIR (Near-Infrared Radiation)
            Diffuse_NIR = NIR * Diffuse_Light_Fraction
            Direct_NIR = NIR - Diffuse_NIR
            
            # Leaf Blade Angle and Direct Beam Extinction Coefficient
            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            

            
            # Diffuse Extinction Coefficients for PAR and NIR
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
            
            # Reflection Coefficients for PAR and NIR
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            Scattered_Beam_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient)

            
            # Absorption of PAR and NIR by Sunlit Leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            Absorbed_NIR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_NIR, Diffuse_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR, Canopy_Diffuse_Reflection_Coefficient_NIR, Direct_NIR, Diffuse_NIR, Leaf_Area_Index)
            
            # Total Absorbed Radiation (PAR + NIR) by Sunlit Leaves
            Total_Absorbed_Radiation_Sunlit = Absorbed_PAR_Sunlit + Absorbed_NIR_Sunlit
            
            # Fraction of Sunlit and Shaded Components in Canopy
            Sunlit_Fraction = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            
            # Boundary Layer Resistance for Heat and Water Vapor for Sunlit Leaves
            Boundary_Layer_Conductance_Heat=(0.01 * np.sqrt(Wind_Speed / self.Leaf_object.Leaf_Width))
            Boundary_Layer_Conductance_Heat_Sunlit = Boundary_Layer_Conductance_Heat*(1 - np.exp(-(0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient)
            Boundary_Layer_Resistance_Heat_Sunlit = 1. / Boundary_Layer_Conductance_Heat_Sunlit
            Boundary_Layer_Resistance_Water_Sunlit = 0.93 * Boundary_Layer_Resistance_Heat_Sunlit

                    
            Sat_Vapor_Pressure_Ambient, _ = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            Sat_Vapor_Pressure_Leaf, Intercellular_CO2_Leaf = Leaf.INTERNAL_CO2(Sunlit_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            
            Vapour_Pressure_Deficit = max(0, Sat_Vapor_Pressure_Ambient - Vapour_Pressure)
            
            # Dynamic slope calculation for the response of saturation vapor pressure to temperature
            Slope_VPD = (Sat_Vapor_Pressure_Leaf - Sat_Vapor_Pressure_Ambient) / self.Avoid_Zero_Division(Sunlit_Leaf_Temp - Hourly_Temp)
            
            # Calculating potential CO2 conductance for sunlit leaves based on photosynthesis minus dark respiration
            Potential_CO2_Conductance_Sunlit = (Photosynthesis_Sunlit - Dark_Respiration_Sunlit) * (273.15 + Sunlit_Leaf_Temp) / 0.53717 / (self.Leaf_object.Ambient_CO2 - Intercellular_CO2_Leaf)
            
            # Turbulence resistance calculation based on canopy structure and wind speed
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            Turbulence_Resistance_Sunlit=Turbulence_Resistance * Sunlit_Fraction
            # Calculating sunlit Leaf-specific stomatal resistance to water vapor
            Stomatal_Resistance_Water_Sunlit = max(1E-30, 1 / Potential_CO2_Conductance_Sunlit - Boundary_Layer_Resistance_Water_Sunlit * 1.3 - Turbulence_Resistance_Sunlit) / 1.6
            
            # Applying the Penman-Monteith equation for potential transpiration estimation
            Potential_Transpiration_Sunlit, Absorbed_Radiation_Sunlit = Leaf.Penman_Monteith(Stomatal_Resistance_Water_Sunlit, Turbulence_Resistance_Sunlit, Boundary_Layer_Resistance_Water_Sunlit, Boundary_Layer_Resistance_Heat_Sunlit, Total_Absorbed_Radiation_Sunlit, Atmospheric_Transmissivity, Sunlit_Fraction, Sunlit_Leaf_Temp, Vapour_Pressure, Slope_VPD, Vapour_Pressure_Deficit)
            #print(Potential_Transpiration_Sunlit)
            # Updating the Leaf object with calculated hourly and potential daily transpiration rates, absorbed radiation, etc.
            Hourly_Transpiration_Sunlit.append(Potential_Transpiration_Sunlit)
            Hourly_Absorbed_Radiation_Sunlit.append(Absorbed_Radiation_Sunlit)
            Hourly_Stomatal_Resistance_Water_Sunlit.append(Stomatal_Resistance_Water_Sunlit)
            Hourly_Slope_VPD_Sunlit.append(Slope_VPD)
            Hourly_Absorbed_PAR_Sunlit.append(Absorbed_PAR_Sunlit)
            
        # Final aggregation to daily values might follow here
        # self.Leaf_object.Daily_Potential_Transpiration_Sunlit = Leaf.aggregate_to_daily(Hourly_Transpiration_Sunlit, Day_Length)
        # self.Leaf_object.Daily_Absorbed_Radiation_Sunlit = Leaf.aggregate_to_daily(Hourly_Absorbed_Radiation_Sunlit, Day_Length)
        
                
        self.Leaf_object.Hourly_Transpiration_Sunlit = Hourly_Transpiration_Sunlit
        self.Leaf_object.Hourly_Absorbed_Radiation_Sunlit = Hourly_Absorbed_Radiation_Sunlit
        self.Leaf_object.Hourly_Stomatal_Resistance_Water_Sunlit = Hourly_Stomatal_Resistance_Water_Sunlit
        self.Leaf_object.Hourly_Slope_VPD_Sunlit = Hourly_Slope_VPD_Sunlit
        self.Leaf_object.Hourly_Absorbed_PAR_Sunlit = Hourly_Absorbed_PAR_Sunlit
        
        # Aggregate hourly data back to daily totals for a comprehensive daily overview
        # Daily_Absorbed_Radiation_Sunlit = Leaf.aggregate_to_daily(Hourly_Absorbed_Radiation_Sunlit, Day_Length)
        # self.Leaf_object.Daily_Potential_Transpiration_Sunlit = Leaf.aggregate_to_daily(Hourly_Transpiration_Sunlit, Day_Length)
        # Optionally, aggregate other metrics like absorbed PAR if needed for further analysis
        # self.Leaf_object.Daily_Absorbed_PAR_Sunlit = Leaf.aggregate_to_daily(Hourly_Absorbed_PAR_Sunlit, Day_Length)
        
        
        
        
    def Update_LeafTemp_Photosynthesis_if_WaterStress(self,water_supply_for_Transpiration,Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,
                                                      Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height, average_root_zone_water_content,
                                                      water_supply_for_evaporation, Root_Depth,
                                                      Hourly_Sunlit_Leaf_Temp, Hourly_Shaded_Leaf_Temp,
                                                      Potential_Canopy_Transpiration,Actual_Soil_Evaporation, Hourly_Soil_Evaporation,C3C4_Pathway):
        Gaussian_Points = 5
        Gaussian_Weights = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
    
        # Using updated attributes
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Leaf_Nitrogen_Extinction_Coefficient = self.Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff']
        Wind_Ext_Coeff = self.Leaf_object.Leaf_area_output['Wind_Ext_Coeff']
        Specific_Leaf_N_Top = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top']
        Specific_Leaf_N_Top_Increment = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top_Increment']
    
        Hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)
    
        Actual_Photosynthesis_list = []
        Actual_Photosynthesis_DELTA_list = []
        Actual_Transpiration_list = []
        Actual_Air_Leaf_Temperature_Difference_list = []
        Actual_Leaf_Temperature_list = []
        # print(water_supply_for_evaporation,Actual_Soil_Evaporation)

        # Evaporation_layer_leftover=max(0,(water_supply_for_evaporation-Actual_Soil_Evaporation))
        # Maximum_Possible_Transpiration = max(1e-32, 1000*average_root_zone_water_content -water_supply_for_evaporation+Evaporation_layer_leftover)
        Maximum_Possible_Transpiration=water_supply_for_Transpiration
        
        if Potential_Canopy_Transpiration > Maximum_Possible_Transpiration :
            Water_Stress_Fraction=Maximum_Possible_Transpiration/Potential_Canopy_Transpiration
            # print('L 776: ' ,Water_Stress_Fraction,Potential_Canopy_Transpiration, Maximum_Possible_Transpiration,self.Leaf_object.Hourly_Transpiration_Sunlit)

            Actual_Hourly_Transpiration_Sunlit=np.array(self.Leaf_object.Hourly_Transpiration_Sunlit)*Water_Stress_Fraction
            Actual_Hourly_Transpiration_Shaded=np.array(self.Leaf_object.Hourly_Transpiration_Shaded)*Water_Stress_Fraction

        else:
            Water_Stress_Fraction=1
            Actual_Hourly_Transpiration_Sunlit=np.array(self.Leaf_object.Hourly_Transpiration_Sunlit)
            Actual_Hourly_Transpiration_Shaded=np.array(self.Leaf_object.Hourly_Transpiration_Shaded)
        # print('L 785: ' ,Water_Stress_Fraction, self.Leaf_object.Hourly_Transpiration_Sunlit)

        Hourly_Transpiration_Sunlit=[]
        for i, (_, Hourly_Temp, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed), Transpiration_Sunlit, Transpiration_Shaded,Actual_Transpiration_Sunlit,Actual_Transpiration_Shaded, Absorbed_Radiation_Sunlit, Stomatal_Resistance_Water_Sunlit, Slope_VPD, Soil_Evaporation in zip(range(Gaussian_Points), Hourly_data,
                                                                                                                                                                                                                  self.Leaf_object.Hourly_Transpiration_Sunlit,
                                                                                                                                                                                                                  self.Leaf_object.Hourly_Transpiration_Shaded,
                                                                                                                                                                                                                  Actual_Hourly_Transpiration_Sunlit,
                                                                                                                                                                                                                  Actual_Hourly_Transpiration_Shaded,
                                                                                                                                                                                                                  self.Leaf_object.Hourly_Absorbed_Radiation_Sunlit,
                                                                                                                                                                                                                  self.Leaf_object.Hourly_Stomatal_Resistance_Water_Sunlit,
                                                                                                                                                                                                                  self.Leaf_object.Hourly_Slope_VPD_Sunlit,
                                                                                                                                                                                                                  Hourly_Soil_Evaporation):
            Hourly_Transpiration_Sunlit.append(Actual_Transpiration_Sunlit)
            Hour = 12 - 0.5 * Day_Length + Day_Length * Gaussian_Weights[i]
            Sin_Beam_Adjusted = max(0., Sin_Beam + Cos_Solar_Declination * np.cos(2. * np.pi * (Hour - 12.) / 24.))
            
            # Root_Zone_Water_Supply_Hourly = 1000* average_root_zone_water_content/(3600*Day_Length) 
            # water_supply_for_evaporation_Hourly = water_supply_for_evaporation/(3600*Day_Length)

            # Evaporation_layer_leftover=max(0,(water_supply_for_evaporation_Hourly-Soil_Evaporation))
            # Maximum_Possible_Transpiration = max(1e-32,Root_Zone_Water_Supply_Hourly -water_supply_for_evaporation_Hourly+Evaporation_layer_leftover)

            # # Total potential canopy transpiration for the hour
            # Total_Potential_Canopy_Transpiration = Transpiration_Sunlit + Transpiration_Shaded
            # print(Maximum_Possible_Transpiration)


            # # Decision-making based on the comparison between the maximum possible transpiration and the total potential canopy transpiration
            # if Maximum_Possible_Transpiration < Total_Potential_Canopy_Transpiration:
            #     Actual_Canopy_Transpiration = Maximum_Possible_Transpiration
            # else:
            #     Actual_Canopy_Transpiration = Total_Potential_Canopy_Transpiration
            
            # This ensures that the plant's transpiration rates are adjusted in real-time based on soil water availability, reflecting a realistic physiological response to varying water stress conditions. The model thus accounts for the dynamic interplay between environmental water availability a
            # Actual_Canopy_Transpiration = Transpiration_Sunlit + Transpiration_Shaded
            
            # Adjusting actual transpiration for sunlit leaves based on canopy water availability
            # Actual_Transpiration_Sunlit = Transpiration_Sunlit * (Actual_Canopy_Transpiration / Total_Potential_Canopy_Transpiration)
            # Actual_Transpiration_Sunlit = Transpiration_Sunlit

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
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            
            # Extinction coefficients for sunlight penetration through the canopy
            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            # Scattering coefficient for PAR and adjustments for Leaf and canopy level interactions

            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            

            
            # Calculating boundary layer resistance to heat and water vapor for sunlit leaves
            Boundary_Layer_Conductance_Heat=(0.01 * np.sqrt(Wind_Speed / self.Leaf_object.Leaf_Width))
            Boundary_Layer_Conductance_Heat_Sunlit = Boundary_Layer_Conductance_Heat*(1 - np.exp(-(0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Ext_Coeff + Direct_Beam_Extinction_Coefficient)
            Boundary_Layer_Resistance_Heat_Sunlit = 1. / Boundary_Layer_Conductance_Heat_Sunlit
            Boundary_Layer_Resistance_Water_Sunlit = 0.93 * Boundary_Layer_Resistance_Heat_Sunlit

            
            # Calculating the fraction of sunlit and shaded components in the canopy
            Sunlit_Fraction = 1.0 / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1 - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            
            # Adjusting turbulence resistance for the canopy
            Turbulence_Resistance = 0.74 * (np.log((2 - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            
            # Calculating Leaf temperature adjustments due to water stress
            Turbulence_Resistance_Sunlit = Turbulence_Resistance * Sunlit_Fraction
            Temperature_Difference =  (Absorbed_Radiation_Sunlit - Latent_Heat_Vaporization * Actual_Transpiration_Sunlit) * (Boundary_Layer_Resistance_Heat_Sunlit + Turbulence_Resistance_Sunlit) / Volumetric_Heat_Capacity_Air
            Temperature_Difference=max(-25, min(Temperature_Difference, 25)) 

            Adjusted_Leaf_Temperature = Hourly_Temp + Temperature_Difference
                        
            # Adjusting stomatal resistance to water under water stress conditions
            Adjusted_Stomatal_Resistance_Water = (Transpiration_Sunlit - Actual_Transpiration_Sunlit) * (Slope_VPD * (Boundary_Layer_Resistance_Heat_Sunlit + Turbulence_Resistance_Sunlit) + Psychrometric_Constant * (Boundary_Layer_Resistance_Water_Sunlit + Turbulence_Resistance_Sunlit)) / Actual_Transpiration_Sunlit / Psychrometric_Constant + Transpiration_Sunlit / Actual_Transpiration_Sunlit * Stomatal_Resistance_Water_Sunlit
            
            # Absorbed PAR calculation for sunlit leaves
            Absorbed_PAR_Sunlit, _ = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            
            # Adjusting photosynthetic nitrogen for sunlit parts of the canopy
            Photosynthetic_Nitrogen_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_Nitrogen_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_Nitrogen_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            Photosynthetic_Nitrogen_Sunlit_DELTA = Specific_Leaf_N_Top_Increment * (1. - np.exp(-(Leaf_Nitrogen_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_Nitrogen_Extinction_Coefficient + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient

            # Calculating internal CO2 and photosynthesis adjustments
            Sat_Vapor_Pressure_Leaf, Intercellular_CO2_Leaf = Leaf.INTERNAL_CO2(Adjusted_Leaf_Temperature, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            
            Actual_Photosynthesis_Rate, Dark_Respiration_Rate = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Adjusted_Leaf_Temperature, Intercellular_CO2_Leaf, Photosynthetic_Nitrogen_Sunlit, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)
            Actual_Photosynthesis_Rate_DELTA, Dark_Respiration_Rate_DELTA = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Sunlit, Adjusted_Leaf_Temperature, Intercellular_CO2_Leaf, Photosynthetic_Nitrogen_Sunlit_DELTA, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)

            # Actual photosynthesis under water stress condition
            Actual_Photosynthesis = (1.6 * Stomatal_Resistance_Water_Sunlit + 1.3 * Boundary_Layer_Resistance_Water_Sunlit + Turbulence_Resistance_Sunlit) / (1.6 * Adjusted_Stomatal_Resistance_Water + 1.3 * Boundary_Layer_Resistance_Water_Sunlit + Turbulence_Resistance_Sunlit) * (Actual_Photosynthesis_Rate - Dark_Respiration_Rate) + Dark_Respiration_Rate
            Actual_Photosynthesis_DELTA = (1.6 * Stomatal_Resistance_Water_Sunlit + 1.3 * Boundary_Layer_Resistance_Water_Sunlit + Turbulence_Resistance_Sunlit) / (1.6 * Adjusted_Stomatal_Resistance_Water + 1.3 * Boundary_Layer_Resistance_Water_Sunlit + Turbulence_Resistance_Sunlit) * (Actual_Photosynthesis_Rate_DELTA - Dark_Respiration_Rate_DELTA) + Dark_Respiration_Rate_DELTA
            # print(Actual_Photosynthesis_DELTA,Actual_Photosynthesis)
            # Appending results for analysis
            Actual_Photosynthesis_list.append(Actual_Photosynthesis)
            Actual_Photosynthesis_DELTA_list.append(Actual_Photosynthesis_DELTA)
            #print(Actual_Photosynthesis_DELTA_list)
            Actual_Transpiration_list.append(Actual_Transpiration_Sunlit)
            Actual_Air_Leaf_Temperature_Difference_list.append(Temperature_Difference)
            Actual_Leaf_Temperature_list.append(Adjusted_Leaf_Temperature)
            #print(Actual_Leaf_Temperature_list)

        # Updating the Leaf model object with the calculated hourly values
        self.Leaf_object.Hourly_Transpiration_Sunlit = Hourly_Transpiration_Sunlit
        # print('L 905: ', self.Leaf_object.Hourly_Transpiration_Sunlit)

        self.Leaf_object.Hourly_Actual_Photosynthesis_Sunlit = Actual_Photosynthesis_list
        self.Leaf_object.Hourly_Actual_Photosynthesis_Sunlit_DELTA = Actual_Photosynthesis_DELTA_list
        # self.Leaf_object.Hourly_Actual_Transpiration_Sunlit = Actual_Transpiration_list
        self.Leaf_object.Hourly_Actual_Air_Sunlit_Leaf_Temp_Diff = Actual_Air_Leaf_Temperature_Difference_list
        self.Leaf_object.Hourly_Actual_Sunlit_Leaf_Temp = Actual_Leaf_Temperature_list
        #print(Actual_Leaf_Temperature_list)

    
    
class Leaf_Shaded(Leaf):
    def __init__(self, Leaf_object):
        self.Leaf_object = Leaf_object
        
    def Calculate_Leaf_temp(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation,
                            Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height,C3C4_Pathway):
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Wind_Extinction_Coeff = self.Leaf_object.Leaf_area_output['Wind_Ext_Coeff']
        Leaf_N_Extinction_Coeff = self.Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff']
        
        Specific_Leaf_N_Top = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top']
        
        Hourly_data, w_Gauss = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)

        Hourly_Air_Temperature_Difference_Shaded = []
        Hourly_Shaded_Leaf_Temperature = []
        
        for Solar_Constant, Hourly_Temp, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed in Hourly_data:
        
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
            NIR = 0.5 * Hourly_Solar_Radiation
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)
        
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            Sat_Vapor_Pressure, Intercellular_CO2 = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
        
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


            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            
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
            Saturated_Vapor_Pressure, Internal_CO2 = Leaf.INTERNAL_CO2(Hourly_Temp, Vapour_Pressure, Vapor_Pressure_Slope, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            
            # Total photosynthetic nitrogen across the entire canopy
            Total_Photosynthetic_N_Canopy = (Specific_Leaf_N_Top * (1. - np.exp(-Leaf_N_Extinction_Coeff * Leaf_Area_Index)) / Leaf_N_Extinction_Coeff - self.Leaf_object.Min_Specific_Leaf_N * Leaf_Area_Index)
            
            # Photosynthetic nitrogen allocated to sunlit and shaded Leaf portions
            Photosynthetic_N_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coeff) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coeff) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coeff * Leaf_Area_Index)) / Direct_Beam_Extinction_Coeff
            Photosynthetic_N_Shaded = Total_Photosynthetic_N_Canopy - Photosynthetic_N_Sunlit
            
            Potential_Photosynthesis_Shaded, Dark_Respiration_Shaded = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Shaded, Hourly_Temp, Internal_CO2, Photosynthetic_N_Shaded, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)
            
            # Calculating the fraction of sunlit and shaded Leaf area
            Fraction_Sunlit = 1. / Direct_Beam_Extinction_Coeff / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coeff * Leaf_Area_Index))
            Fraction_Shaded = 1 - Fraction_Sunlit

            # Boundary layer resistance for canopy, sunlit and shaded leaves
            Boundary_Layer_Leaf_Flow = 0.01 * np.sqrt(Wind_Speed / self.Leaf_object.Leaf_Width)
            Canopy_Boundary_Layer_Flow = (1 - np.exp(-0.5 * Wind_Extinction_Coeff * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff) * Boundary_Layer_Leaf_Flow
            Sunlit_Boundary_Layer_Flow = (1 - np.exp(-(0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coeff) * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coeff) * Boundary_Layer_Leaf_Flow
            Shaded_Boundary_Layer_Flow = Canopy_Boundary_Layer_Flow - Sunlit_Boundary_Layer_Flow
            Boundary_Layer_Heat_Resist_Shaded = 1. / Shaded_Boundary_Layer_Flow  # Boundary layer resistance to heat, shaded part
            Boundary_Layer_Water_Resist_Shaded = 0.93 * Boundary_Layer_Heat_Resist_Shaded  # Boundary layer resistance to water, shaded part   
                    
            
            # Potential conductance for CO2 for shaded leaves
            Conductance_CO2_Shaded = (Potential_Photosynthesis_Shaded - Dark_Respiration_Shaded) * (273. + Hourly_Temp) / 0.53717 / (self.Leaf_object.Ambient_CO2 - Internal_CO2)
            
            # Turbulence resistance for canopy
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            Turbulence_Resistance_Shaded=Turbulence_Resistance * Fraction_Shaded        
            Stomatal_Resist_Water_Shaded = max(1E-30, 1 / Conductance_CO2_Shaded - Boundary_Layer_Water_Resist_Shaded * 1.3 - Turbulence_Resistance_Shaded) / 1.6
            
            Potential_Transpiration_Shaded, Net_Radiation_Absorbed_Shaded = Leaf.Penman_Monteith(Stomatal_Resist_Water_Shaded, Turbulence_Resistance_Shaded, Boundary_Layer_Water_Resist_Shaded, Boundary_Layer_Heat_Resist_Shaded, Absorbed_Total_Radiation_Shaded, Atmospheric_Transmissivity, Fraction_Shaded, Hourly_Temp, Vapour_Pressure, Slope_SVP, Vapour_Pressure_Deficit)
            # Calculate Leaf temperature for shaded leaves
            Temperature_Difference_Shaded = (Net_Radiation_Absorbed_Shaded - Latent_Heat_Vaporization * Potential_Transpiration_Shaded) * (Boundary_Layer_Heat_Resist_Shaded + Turbulence_Resistance_Shaded) / Volumetric_Heat_Capacity_Air
            Temperature_Difference_Shaded=max(-25, min(Temperature_Difference_Shaded, 25)) 

            
            Hourly_Air_Temperature_Difference_Shaded.append(Temperature_Difference_Shaded)
            
            Shaded_Leaf_Temperature = Hourly_Temp + Temperature_Difference_Shaded
            Hourly_Shaded_Leaf_Temperature.append(Shaded_Leaf_Temperature)
            #print(Temperature_Difference_Shaded)
        self.Leaf_object.Hourly_Shaded_Leaf_Temp=Hourly_Shaded_Leaf_Temperature
        self.Leaf_object.Hourly_Air_Shaded_Leaf_Temp_diff=Hourly_Air_Temperature_Difference_Shaded


    def Calculate_Potential_Photosynthesis(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,
                                           Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed,C3C4_Pathway):
        # Use updated attributes for the calculations
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Leaf_N_Extinction_Coeff = self.Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff']
        
        Specific_Leaf_N_Top = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top']
        
        # Convert daily conditions to hourly data for detailed analysis
        Hourly_conditions, w_Gauss = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)
        
        Hourly_Shaded_Leaf_Temp = self.Leaf_object.Hourly_Shaded_Leaf_Temp

        # Initialize lists to store outputs for hourly calculations
        Hourly_Shaded_Leaf_Photosynthesis = []
        Hourly_Shaded_Dark_Respiration = []
        for (Solar_Constant, Daytime_Temperature, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed), Shaded_Leaf_Temp in zip(Hourly_conditions, Hourly_Shaded_Leaf_Temp):

            # Use Shaded_Leaf_Temp where Daytime_Temperature was used
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
            
            # Calculate the diffuse light fraction based on atmospheric transmissivity
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity
            
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
            #print(Solar_Constant , Sin_Solar_Declination)
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR

            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180 
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            
        
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
        
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
        
        
            # Adjusting for vapor pressure deficit influence on intercellular CO2 concentration
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway ==-1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            Sat_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Shaded_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
        
            # Total photosynthetic nitrogen across the entire canopy
            Total_Photosynthetic_N_Canopy = (Specific_Leaf_N_Top * (1. - np.exp(-Leaf_N_Extinction_Coeff * Leaf_Area_Index)) / Leaf_N_Extinction_Coeff - self.Leaf_object.Min_Specific_Leaf_N * Leaf_Area_Index)
            
            # Calculating photosynthetic nitrogen availability for shaded canopy parts
            Photosynthetic_N_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            Photosynthetic_N_Shaded = Total_Photosynthetic_N_Canopy - Photosynthetic_N_Sunlit
        
        
        
            # Absorption of PAR by shaded leaves
            _, Absorbed_PAR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            # print(Absorbed_PAR_Shaded)
            # Calculating potential photosynthesis and dark respiration for shaded leaves
            Potential_Photosynthesis_Shaded, Dark_Respiration_Shaded = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Shaded, Shaded_Leaf_Temp, Intercellular_CO2_Concentration, Photosynthetic_N_Shaded, self.Leaf_object.Activation_Energy_Jmax,self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)
            #print( Photosynthetic_N_Shaded,)
           
            Hourly_Shaded_Leaf_Photosynthesis.append(Potential_Photosynthesis_Shaded)
            Hourly_Shaded_Dark_Respiration.append(Dark_Respiration_Shaded)
            #print(Hourly_Shaded_Leaf_Photosynthesis)

        self.Leaf_object.Hourly_Photosynthesis_Shaded = Hourly_Shaded_Leaf_Photosynthesis
        self.Leaf_object.Hourly_Dark_Respiration_Shaded = Hourly_Shaded_Dark_Respiration
        
        # # Aggregate hourly data back to daily totals
        # Daily_Potential_Photosynthesis_Shaded = Leaf.aggregate_to_daily(Hourly_Shaded_Leaf_Photosynthesis, Day_Length)
        # Daily_Dark_Respiration_Shaded = Leaf.aggregate_to_daily(Hourly_Shaded_Dark_Respiration, Day_Length)
        # #print(Daily_Potential_Photosynthesis_Shaded)

        # self.Leaf_object.Daily_Potential_Photosynthesis_Shaded = Daily_Potential_Photosynthesis_Shaded
        # self.Leaf_object.Daily_Dark_Respiration_Shaded = Daily_Dark_Respiration_Shaded
              
    def Calculate_Potential_Transpiration(self, Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure,
                                          Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height,C3C4_Pathway):
        # Use updated attributes
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Wind_Extinction_Coeff = self.Leaf_object.Leaf_area_output['Wind_Ext_Coeff']
    
    
        Hourly_conditions, w_Gauss = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)

    
        Hourly_Shaded_Leaf_Temp = self.Leaf_object.Hourly_Shaded_Leaf_Temp
        Hourly_Photosynthesis_Shaded = self.Leaf_object.Hourly_Photosynthesis_Shaded
        Hourly_Dark_Respiration_Shaded = self.Leaf_object.Hourly_Dark_Respiration_Shaded
    
        Hourly_Transpiration_Shaded = []
        Hourly_Absorbed_Radiation_Shaded = []
        Hourly_Stomatal_Resistance_Water_Shaded = []
        Hourly_Slope_VPD_Shaded = []
        Hourly_Absorbed_PAR_Shaded = []         
            
        
        for (Solar_Constant, Hourly_Temperature, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed), Shaded_Leaf_Temp, Photosynthesis_Shaded, Dark_Respiration_Shaded in zip(Hourly_conditions, Hourly_Shaded_Leaf_Temp, Hourly_Photosynthesis_Shaded, Hourly_Dark_Respiration_Shaded):
            
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
            # Calculation of diffuse light fraction
            Atmospheric_Transmissivity = Incoming_PAR / (0.5 * Solar_Constant * Sin_Beam)
            
            if Atmospheric_Transmissivity < 0.22:
                Diffuse_Light_Fraction = 1
            elif 0.22 < Atmospheric_Transmissivity <= 0.35:
                Diffuse_Light_Fraction = 1 - 6.4 * (Atmospheric_Transmissivity - 0.22) ** 2
            else:
                Diffuse_Light_Fraction = 1.47 - 1.66 * Atmospheric_Transmissivity
            
            Diffuse_Light_Fraction = max(Diffuse_Light_Fraction, 0.15 + 0.85 * (1 - np.exp(-0.1 / Sin_Beam)))
            
            # Calculating the diffuse and direct components of PAR and NIR
            Diffuse_PAR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_PAR = Incoming_PAR - Diffuse_PAR
            
            Diffuse_NIR = Incoming_PAR * Diffuse_Light_Fraction
            Direct_NIR = Incoming_PAR - Diffuse_NIR

            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)

            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Diffuse_Extinction_Coefficient_NIR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_NIR)
            
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            Scattered_Beam_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient)

            # Calculating the absorption of PAR and NIR by shaded leaves
            _, Absorbed_PAR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            _, Absorbed_NIR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coefficient_NIR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_NIR, Diffuse_Extinction_Coefficient_NIR, Canopy_Beam_Reflection_Coefficient_NIR, Canopy_Diffuse_Reflection_Coefficient_NIR, Direct_NIR, Diffuse_NIR, Leaf_Area_Index)

        
            Total_Absorbed_Radiation_Shaded = Absorbed_PAR_Shaded + Absorbed_NIR_Shaded
    
            Fraction_Sunlit = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            Fraction_Shaded = 1 - Fraction_Sunlit
            
            # Calculating boundary layer resistance for shaded leaves
            Boundary_Layer_Leaf_Flow = 0.01 * np.sqrt(Wind_Speed / self.Leaf_object.Leaf_Width)
            Canopy_Boundary_Layer_Flow = (1 - np.exp(-0.5 * Wind_Extinction_Coeff * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff) * Boundary_Layer_Leaf_Flow
            Sunlit_Boundary_Layer_Flow = (1 - np.exp(-(0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Leaf_Flow
            Shaded_Boundary_Layer_Flow = Canopy_Boundary_Layer_Flow - Sunlit_Boundary_Layer_Flow
            Boundary_Layer_Heat_Resist_Shaded = 1. / Shaded_Boundary_Layer_Flow  # Boundary layer resistance to heat, shaded part
            Boundary_Layer_Water_Resist_Shaded = 0.93 * Boundary_Layer_Heat_Resist_Shaded  # Boundary layer resistance to water, shaded part   
   
            
            
            
            
            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)

            Saturation_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Hourly_Temperature, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            Adjusted_Saturation_Vapor_Pressure, Adjusted_Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Shaded_Leaf_Temp, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            
            Vapor_Pressure_Deficit = max(0, Saturation_Vapor_Pressure - Vapour_Pressure)
            
            Slope_VPD_Response = (Adjusted_Saturation_Vapor_Pressure - Saturation_Vapor_Pressure) / self.Avoid_Zero_Division(Shaded_Leaf_Temp- Hourly_Temperature)
               
            
            
            # Adjusted CO2 conductance for transpiration
            Adjusted_CO2_Conductance = (Photosynthesis_Shaded - Dark_Respiration_Shaded) * (273.15 + Shaded_Leaf_Temp) / 0.53717 / (self.Leaf_object.Ambient_CO2 - Adjusted_Intercellular_CO2_Concentration)
            
            # Turbulence resistance for canopy
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            Turbulence_Resistance_Shaded =Turbulence_Resistance*Fraction_Shaded
            # Adjusted stomatal resistance to water vapor
            Adjusted_Stomatal_Resistance = max(1E-30, 1 / Adjusted_CO2_Conductance - Boundary_Layer_Water_Resist_Shaded * 1.3 - Turbulence_Resistance_Shaded) / 1.6
            
            # Calculating potential transpiration for shaded leaves and absorbed radiation
            Potential_Transpiration_Shaded, Absorbed_Radiation_Shaded = Leaf.Penman_Monteith(Adjusted_Stomatal_Resistance, Turbulence_Resistance_Shaded, Boundary_Layer_Water_Resist_Shaded, Boundary_Layer_Heat_Resist_Shaded, Total_Absorbed_Radiation_Shaded, Atmospheric_Transmissivity, Fraction_Shaded, Shaded_Leaf_Temp, Vapour_Pressure, Slope_VPD_Response, Vapor_Pressure_Deficit)
            #print(Potential_Transpiration_Shaded)
            
            # Appending calculated values to their respective lists
            Hourly_Transpiration_Shaded.append(Potential_Transpiration_Shaded)
            Hourly_Absorbed_Radiation_Shaded.append(Absorbed_Radiation_Shaded)
            Hourly_Stomatal_Resistance_Water_Shaded.append(Adjusted_Stomatal_Resistance)
            Hourly_Slope_VPD_Shaded.append(Slope_VPD_Response)
            Hourly_Absorbed_PAR_Shaded.append(Absorbed_PAR_Shaded)
            
        self.Leaf_object.Hourly_Transpiration_Shaded = Hourly_Transpiration_Shaded
        self.Leaf_object.Hourly_Absorbed_Radiation_Shaded = Hourly_Absorbed_Radiation_Shaded
        self.Leaf_object.Hourly_Stomatal_Resistance_Water_Shaded = Hourly_Stomatal_Resistance_Water_Shaded
        self.Leaf_object.Hourly_Slope_VPD_Shaded = Hourly_Slope_VPD_Shaded
        self.Leaf_object.Hourly_Absorbed_PAR_Shaded = Hourly_Absorbed_PAR_Shaded


        
        # potential_transpiration_Shaded_daily = Leaf.aggregate_to_daily(Hourly_Transpiration_Shaded, Day_Length)
        # Daily_Absorbed_Radiation_Shaded = Leaf.aggregate_to_daily(Hourly_Absorbed_Radiation_Shaded, Day_Length)
        # self.Leaf_object.potential_transpiration_Shaded_daily = potential_transpiration_Shaded_daily
        # self.Leaf_object.Daily_Absorbed_Radiation_Shaded = Daily_Absorbed_Radiation_Shaded
        #print(potential_transpiration_Shaded_daily)                    
    def Update_LeafTemp_Photosynthesis_if_WaterStress(self, water_supply_for_Transpiration,Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Vapour_Pressure, Wind_Speed, Plant_Height,
                                                      average_root_zone_water_content,water_supply_for_evaporation, Root_Depth,
                                                      Hourly_Sunlit_Leaf_Temp, Hourly_Shaded_Leaf_Temp,
                                                      Potential_Canopy_Transpiration,Actual_Soil_Evaporation,Hourly_Soil_Evap,C3C4_Pathway):
    
        Gauss_Points = 5
        Gauss_Weights = np.array([0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899])
        # Now use updated attributes
        Total_LAI = self.Leaf_object.Leaf_area_output['Total_LAI']
        Leaf_Area_Index = self.Leaf_object.Leaf_area_output['Leaf_Area_Index']
        Leaf_N_Extinction_Coeff = self.Leaf_object.Leaf_area_output['Leaf_Nitro_Ext_Coeff']
        Wind_Extinction_Coeff = self.Leaf_object.Leaf_area_output['Wind_Ext_Coeff']
        Specific_Leaf_N_Top = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top']
        Specific_Leaf_N_Top_Increment = self.Leaf_object.specific_Leaf_n_output['Specific_Leaf_N_Top_Increment']

        # Convert daily data to hourly
        Hourly_data, _ = Leaf.convert_daily_to_hourly(Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Daily_Sin_Beam_Exposure, Solar_Radiation, Max_Temp, Min_Temp, Wind_Speed)

        
        Actual_Photosynthesis = []
        Actual_Photosynthesis_DELTA = []
        Actual_Transpiration = []
        Actual_Air_Shaded_Leaf_Temp_Difference = []
        Actual_Shaded_Leaf_Temperature = []

        # Evaporation_layer_leftover=max(0,(irrigation_amount+water_supply_for_evaporation-Actual_Soil_Evaporation))
        # if Evaporation_layer_leftover < -0.00001:
        #     print(water_supply_for_evaporation,Actual_Soil_Evaporation,Evaporation_layer_leftover)
        #     raise ValueError('Evaporation_layer_leftover does not match up!')
        # Maximum_Possible_Transpiration = max(1e-32, 1000*average_root_zone_water_content -water_supply_for_evaporation+Evaporation_layer_leftover)
        Maximum_Possible_Transpiration=water_supply_for_Transpiration
        if Potential_Canopy_Transpiration > Maximum_Possible_Transpiration :
            Water_Stress_Fraction=Maximum_Possible_Transpiration/Potential_Canopy_Transpiration
            Actual_Hourly_Transpiration_Sunlit=np.array(self.Leaf_object.Hourly_Transpiration_Sunlit)*Water_Stress_Fraction
            Actual_Hourly_Transpiration_Shaded=np.array(self.Leaf_object.Hourly_Transpiration_Shaded)*Water_Stress_Fraction
        else:
            Actual_Hourly_Transpiration_Sunlit=np.array(self.Leaf_object.Hourly_Transpiration_Sunlit)
            Actual_Hourly_Transpiration_Shaded=np.array(self.Leaf_object.Hourly_Transpiration_Shaded)
        
        Hourly_Transpiration_Shaded=[]
        for i, (Solar_Constant, Hourly_Temp, Sin_Beam, Hourly_Solar_Radiation, Wind_Speed), Transpiration_Sunlit, Transpiration_Shaded,Actual_Transpiration_Sunlit,Actual_Transpiration_Shaded, Radiation_Shaded, Stomatal_Resist_Water_Shaded, Slope_VPD, Soil_Evaporation in zip(range(Gauss_Points), Hourly_data,  
                                                                                                                                                                                                                                             self.Leaf_object.Hourly_Transpiration_Sunlit,
                                                                                                                                                                                                                                             self.Leaf_object.Hourly_Transpiration_Shaded, 
                                                                                                                                                                                                                                             Actual_Hourly_Transpiration_Sunlit,
                                                                                                                                                                                                                                             Actual_Hourly_Transpiration_Shaded,
                                                                                                                                                                                                                                             self.Leaf_object.Hourly_Absorbed_Radiation_Shaded,
                                                                                                                                                                                                                                             self.Leaf_object.Hourly_Stomatal_Resistance_Water_Shaded,
                                                                                                                                                                                                                                             self.Leaf_object.Hourly_Slope_VPD_Shaded,
                                                                                                                                                                                                                                             Hourly_Soil_Evap):
            Hourly_Transpiration_Shaded.append(Actual_Transpiration_Shaded)
            Hour = 12 - 0.5 * Day_Length + Day_Length * Gauss_Weights[i]
            Sin_Beam_Adjusted = max(0., Sin_Beam + Cos_Solar_Declination * np.cos(2. * np.pi * (Hour - 12.) / 24.))
            
            # Diurnal availability of soil water supply 
            
            # Root_Zone_Water_Supply_Hourly = 1000* average_root_zone_water_content* (Sin_Beam_Adjusted * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
            # water_supply_for_evaporation_Hourly = water_supply_for_evaporation * (Sin_Beam_Adjusted * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure

            # Root_Zone_Water_Supply_Hourly = 1000* average_root_zone_water_content/(3600*Day_Length) 
            # water_supply_for_evaporation_Hourly = water_supply_for_evaporation/(3600*Day_Length)
            
            # Evaporation_layer_leftover=max(0,(water_supply_for_evaporation_Hourly-Soil_Evaporation))
            # Maximum_Possible_Transpiration = max(1e-32,Root_Zone_Water_Supply_Hourly -water_supply_for_evaporation_Hourly+Evaporation_layer_leftover)

            # # Total potential canopy transpiration for the hour
            # Potential_Canopy_Transpiration = Transpiration_Sunlit + Transpiration_Shaded
                   
 
            # if Maximum_Possible_Transpiration < Potential_Canopy_Transpiration:
            #     Actual_Canopy_Transpiration = Maximum_Possible_Transpiration
            # else:
            #     Actual_Canopy_Transpiration = Potential_Canopy_Transpiration


            # Actual_Transpiration_Shaded = Transpiration_Shaded * (Actual_Canopy_Transpiration / Potential_Canopy_Transpiration)   # Actual transpiration of shaded leaves mm s-1
            # print('Actual_Canopy_Transpiration:',Actual_Canopy_Transpiration)
            # Actual_Transpiration_Shaded = Actual_Transpiration_Sunlit
            Incoming_PAR = 0.5 * Hourly_Solar_Radiation
            
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

            Vapor_Pressure_Deficit_Response = 0.195127 if C3C4_Pathway == -1 else 0.116214  # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-Leaf vapour pressure deficit)
            
            Leaf_Blade_Angle_Radians = self.Leaf_object.Leaf_Blade_Angle * np.pi / 180
            Direct_Beam_Extinction_Coefficient = Leaf.KDR_Coeff(Sin_Beam, Leaf_Blade_Angle_Radians)
            
            
            Diffuse_Extinction_Coefficient_PAR = Leaf.KDF_Coeff(Total_LAI, Leaf_Blade_Angle_Radians, Scattering_Coefficient_PAR)
            Scattered_Beam_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR = Leaf.REFLECTION_Coeff(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient)
            
            
            Fraction_Sunlit = 1. / Direct_Beam_Extinction_Coefficient / Leaf_Area_Index * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index))
            Fraction_Shaded = 1 - Fraction_Sunlit
            
            
            # Calculating boundary layer resistance for shaded leaves
            Boundary_Layer_Leaf_Flow = 0.01 * np.sqrt(Wind_Speed / self.Leaf_object.Leaf_Width)
            Canopy_Boundary_Layer_Flow = (1 - np.exp(-0.5 * Wind_Extinction_Coeff * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff) * Boundary_Layer_Leaf_Flow
            Sunlit_Boundary_Layer_Flow = (1 - np.exp(-(0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (0.5 * Wind_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Boundary_Layer_Leaf_Flow
            Shaded_Boundary_Layer_Flow = Canopy_Boundary_Layer_Flow - Sunlit_Boundary_Layer_Flow
            Boundary_Layer_Heat_Resist_Shaded = 1. / Shaded_Boundary_Layer_Flow  # Boundary layer resistance to heat, shaded part
            Boundary_Layer_Water_Resist_Shaded = 0.93 * Boundary_Layer_Heat_Resist_Shaded  # Boundary layer resistance to water, shaded part   
   
            # Turbulence resistance for canopy
            Turbulence_Resistance = 0.74 * (np.log((2. - 0.7 * Plant_Height) / (0.1 * Plant_Height))) ** 2 / (0.4 ** 2 * Wind_Speed)
            
            #print(Turbulence_Resistance)
            
            # Leaf temperature if water stress occurs
            Turbulence_Resistance_Shaded = Turbulence_Resistance * Fraction_Shaded
            Temperature_Difference =  (Radiation_Shaded - Latent_Heat_Vaporization * Actual_Transpiration_Shaded) * (Boundary_Layer_Heat_Resist_Shaded + Turbulence_Resistance_Shaded) / Volumetric_Heat_Capacity_Air
            Temperature_Difference=max(-25, min(Temperature_Difference, 25)) 

            Adjusted_Leaf_Temperature = Hourly_Temp + Temperature_Difference
            


            # Stomatal resistance to water if water stress occurs
            Adjusted_Stomatal_Resistance_Water_Stress = (Transpiration_Shaded - Actual_Transpiration_Shaded) * (Slope_VPD * (Boundary_Layer_Heat_Resist_Shaded + Turbulence_Resistance_Shaded) + Psychrometric_Constant * (Boundary_Layer_Water_Resist_Shaded + Turbulence_Resistance_Shaded)) / Actual_Transpiration_Shaded / Psychrometric_Constant + Transpiration_Shaded / Actual_Transpiration_Shaded * Stomatal_Resist_Water_Shaded


            # Absorbed PAR by shaded leaves
            _, Absorbed_PAR_Shaded = Leaf.LIGHT_ABSORB(Scattering_Coefficient_PAR, Direct_Beam_Extinction_Coefficient, Scattered_Beam_Extinction_Coefficient_PAR, Diffuse_Extinction_Coefficient_PAR, Canopy_Beam_Reflection_Coefficient_PAR, Canopy_Diffuse_Reflection_Coefficient_PAR, Direct_PAR, Diffuse_PAR, Leaf_Area_Index)
            
            # Total photosynthetic nitrogen in the canopy
            Photosynthetic_Nitrogen_Canopy =(Specific_Leaf_N_Top * (1. - np.exp(-Leaf_N_Extinction_Coeff * Leaf_Area_Index)) / Leaf_N_Extinction_Coeff - self.Leaf_object.Min_Specific_Leaf_N * Leaf_Area_Index)
            # Photosynthetic nitrogen for shaded parts of the canopy
            Photosynthetic_Nitrogen_Sunlit = Specific_Leaf_N_Top * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            Photosynthetic_Nitrogen_Shaded = Photosynthetic_Nitrogen_Canopy - Photosynthetic_Nitrogen_Sunlit
            
            Photosynthetic_Nitrogen_Canopy_DELTA =(Specific_Leaf_N_Top_Increment * (1. - np.exp(-Leaf_N_Extinction_Coeff * Leaf_Area_Index)) / Leaf_N_Extinction_Coeff - self.Leaf_object.Min_Specific_Leaf_N * Leaf_Area_Index)
            Photosynthetic_Nitrogen_Sunlit_DELTA = Specific_Leaf_N_Top_Increment * (1. - np.exp(-(Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) * Leaf_Area_Index)) / (Leaf_N_Extinction_Coeff + Direct_Beam_Extinction_Coefficient) - self.Leaf_object.Min_Specific_Leaf_N * (1. - np.exp(-Direct_Beam_Extinction_Coefficient * Leaf_Area_Index)) / Direct_Beam_Extinction_Coefficient
            Photosynthetic_Nitrogen_Shaded_DELTA = Photosynthetic_Nitrogen_Canopy_DELTA - Photosynthetic_Nitrogen_Sunlit_DELTA

            
            # Calculating internal CO2 concentration and saturation vapor pressure
            Saturation_Vapor_Pressure, Intercellular_CO2_Concentration = Leaf.INTERNAL_CO2(Adjusted_Leaf_Temperature, Vapour_Pressure, Vapor_Pressure_Deficit_Response, self.Leaf_object.Ambient_CO2, self.Leaf_object.C3C4_Pathway)
            
            Actual_Photosynthesis_Rate, Dark_Respiration = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Shaded, Adjusted_Leaf_Temperature, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Shaded, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)
            Actual_Photosynthesis_Rate_DELTA, Dark_Respiration_DELTA = Leaf.PHOTOSYN(self.Leaf_object.C3C4_Pathway, Absorbed_PAR_Shaded, Adjusted_Leaf_Temperature, Intercellular_CO2_Concentration, Photosynthetic_Nitrogen_Shaded_DELTA, self.Leaf_object.Activation_Energy_Jmax, self.Leaf_object.Vcmax_LeafN_Slope, self.Leaf_object.Jmax_LeafN_Slope, self.Leaf_object.Photosynthetic_Light_Response_Factor)

            
            
            # print(Absorbed_PAR_Leaf, Dark_Respiration)
            # Actual photosynthesis under water stress condition
            Actual_Photosynthesis_Water_Stress = (1.6 * Stomatal_Resist_Water_Shaded + 1.3 * Boundary_Layer_Water_Resist_Shaded + Turbulence_Resistance_Shaded) / (1.6 * Adjusted_Stomatal_Resistance_Water_Stress + 1.3 * Boundary_Layer_Water_Resist_Shaded + Turbulence_Resistance_Shaded) * (Actual_Photosynthesis_Rate - Dark_Respiration) + Dark_Respiration
            Actual_Photosynthesis_Water_Stress_DELTA = (1.6 * Stomatal_Resist_Water_Shaded + 1.3 * Boundary_Layer_Water_Resist_Shaded + Turbulence_Resistance_Shaded) / (1.6 * Adjusted_Stomatal_Resistance_Water_Stress + 1.3 * Boundary_Layer_Water_Resist_Shaded + Turbulence_Resistance_Shaded) * (Actual_Photosynthesis_Rate_DELTA - Dark_Respiration_DELTA) + Dark_Respiration_DELTA


            # Appending the calculated values to their respective lists
            Actual_Photosynthesis.append(Actual_Photosynthesis_Water_Stress)
            Actual_Photosynthesis_DELTA.append(Actual_Photosynthesis_Water_Stress_DELTA)
            Actual_Transpiration.append(Actual_Transpiration_Shaded)
            Actual_Air_Shaded_Leaf_Temp_Difference.append(Temperature_Difference)
            Actual_Shaded_Leaf_Temperature.append(Adjusted_Leaf_Temperature)
            # Updating the Leaf model object with the calculated hourly values

        self.Leaf_object.Hourly_Transpiration_Shaded = Hourly_Transpiration_Shaded
        self.Leaf_object.Hourly_Actual_Photosynthesis_Shaded = Actual_Photosynthesis
        self.Leaf_object.Hourly_Actual_Photosynthesis_Shaded_DELTA = Actual_Photosynthesis_DELTA
        # self.Leaf_object.Hourly_Actual_Transpiration_Shaded = Actual_Transpiration
        self.Leaf_object.Hourly_Actual_Air_Shaded_Leaf_Temp_Diff = Actual_Air_Shaded_Leaf_Temp_Difference
        self.Leaf_object.Hourly_Actual_Shaded_Leaf_Temp = Actual_Shaded_Leaf_Temperature
        #print("******************************************************************************************")
        #print(Actual_Transpiration)






 # # Root_Zone_Water_Supply_Hourly = 1000* average_root_zone_water_content* (Sin_Beam_Adjusted * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
 # # water_supply_for_evaporation_Hourly = water_supply_for_evaporation * (Sin_Beam_Adjusted * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure

 # instant_canopy_potential_transpiration=(Transpiration_Shaded + Transpiration_Sunlit)
 # Tr_Fraction=(instant_canopy_potential_transpiration/total_instant_transpiration)

 # instant_canopy_potential_transpiration=instant_canopy_potential_transpiration*3600*gaussian_weight*Day_Length #convert it to total mm transpiration during the day
 # instant_soil_evaporation=Soil_Evaporation*3600*gaussian_weight*Day_Length
 # # # Diurnal availability of soil water supply
 # # water_supply_for_evaporation_Hourly = self.water_supply_for_evaporation * (Hourly_sin_beam * Solar_Constant / 1367) / Daily_Sin_Beam_Exposure
 # # water_supply_for_evaporation_Hourly = water_supply_Hourly * self.evaporation_depth / Root_Depth
 # water_supply_for_evaporation_Hourly = self.water_supply_for_evaporation *Tr_Fraction
 # Root_Zone_Water_Supply_Hourly = 1000* average_root_zone_water_content*Tr_Fraction


 # # Water available for evaporation from the top soil layer
 # # Water_Supply_Evaporation = Water_Supply_Hourly * (evaporation_depth / Root_Depth)
 # Maximum_Possible_Transpiration = Root_Zone_Water_Supply_Hourly -water_supply_for_evaporation_Hourly+(water_supply_for_evaporation_Hourly-instant_soil_evaporation)

 # # # Total potential canopy transpiration for the hour
 # # Potential_Canopy_Transpiration = Transpiration_Sunlit + Transpiration_Shaded
        

 # if Maximum_Possible_Transpiration < instant_canopy_potential_transpiration:
 #     Actual_Canopy_Transpiration = Maximum_Possible_Transpiration
 # else:
 #     Actual_Canopy_Transpiration = instant_canopy_potential_transpiration


 # Actual_Transpiration_Shaded = Transpiration_Shaded * (Actual_Canopy_Transpiration / instant_canopy_potential_transpiration)   # Actual transpiration of shaded leaves mm s-1









