import math
import Leaf
import numpy as np

#initializing intermediate variables
Germination_Efficiency = 0.25  # Growth efficiency
CarbonFrac_Veg = 0.48  # Carbon fraction in vegetative biomass
FractionProtein_StorageOrgans = 0.13  # Fraction of protein in storage organs
FractionCarbs_StorageOrgans = 0.71  # Fraction of carbohydrates in storage organs
CarbonFraction_Organs = 0.47  # Carbon fraction in the storage organs
Growth_Efficiency_Veg = 0.8  # Growth efficiency for storage organs




class Canopy:
    def __init__(self, BTemp_Phen, OTemp_Phen, CTemp_Phen, TempCurve_Res, CropType_Photoperiod, StartPhotoperiod_Phase, EndPhotoperiod_Phase, Photoperiod_Sensitivity,
                 MinThermal_Day_Veg, MinThermal_Day_Rep,
                 CarbonAlloc_Shoot, NitrogenAlloc_Shoot, Plant_Density, Seed_Weight, 
                 Crop_TypeDet, EndSeedNum_DetPeriod, SeedN_RemobFract, SLA_Const, Min_Specific_Leaf_N, MinRootN_Conc,
                 CarbonCost_NFix, MaxN_Uptake, Stem_Nitrogen,
                 IniLeafN_Conc, MaxPlant_Height,
                 legume,
                 MaxStemGrowth_DS, MaxSeedGrowth_DS, StemDW_Height, Model_TimeStep):
        # Thermal Units
        self.BTemp_Phen = BTemp_Phen
        self.OTemp_Phen = OTemp_Phen
        self.CTemp_Phen = CTemp_Phen
        self.TempCurve_Res = TempCurve_Res

        # Development rate
        self.CropType_Photoperiod = CropType_Photoperiod
        self.StartPhotoperiod_Phase = StartPhotoperiod_Phase
        self.EndPhotoperiod_Phase = EndPhotoperiod_Phase
        self.Photoperiod_Sensitivity = Photoperiod_Sensitivity
        self.MinThermal_Day_Veg = MinThermal_Day_Veg
        self.MinThermal_Day_Rep = MinThermal_Day_Rep

        # Calculate biomass
        self.CarbonAlloc_Shoot = CarbonAlloc_Shoot 
        self.NitrogenAlloc_Shoot = NitrogenAlloc_Shoot  
        self.Plant_Density = Plant_Density
        self.Seed_Weight = Seed_Weight
        self.MaxPlant_Height = MaxPlant_Height

        # Seed properties
        self.Crop_TypeDet = Crop_TypeDet
        self.EndSeedNum_DetPeriod = EndSeedNum_DetPeriod
        self.SeedN_RemobFract = SeedN_RemobFract
        self.SLA_Const = SLA_Const
        self.Min_Specific_Leaf_N = Min_Specific_Leaf_N
        self.MinRootN_Conc = MinRootN_Conc
        self.MaxSeedGrowth_DS = MaxSeedGrowth_DS
        self.MaxStemGrowth_DS = MaxStemGrowth_DS
        self.StemDW_Height = StemDW_Height
        self.CarbonCost_NFix = CarbonCost_NFix
        self.MaxN_Uptake = MaxN_Uptake
        self.legume = legume
        self.Stem_Nitrogen = Stem_Nitrogen
        self.IniLeafN_Conc = IniLeafN_Conc
        
        #initializing intermediate variables
        self.Germination_Efficiency = Germination_Efficiency  # Growth efficiency
        self.CarbonFrac_Veg = CarbonFrac_Veg  # Carbon fraction in vegetative biomass
        self.FractionProtein_StorageOrgans = FractionProtein_StorageOrgans  # Fraction of protein in storage organs
        self.FractionCarbs_StorageOrgans = FractionCarbs_StorageOrgans  # Fraction of carbohydrates in storage organs
        self.CarbonFraction_Organs = CarbonFraction_Organs  # Carbon fraction in the storage organs
        self.Growth_Efficiency_Veg = Growth_Efficiency_Veg  # Growth efficiency for storage organs
        
        self.Initial_Leaf_Carbon = self.Plant_Density * self.Seed_Weight * self.CarbonFraction_Organs * self.Germination_Efficiency * self.CarbonAlloc_Shoot
        self.Initial_Root_Carbon = self.Plant_Density * self.Seed_Weight * self.CarbonFraction_Organs * self.Germination_Efficiency * (1.0 - self.CarbonAlloc_Shoot)
        self.Initial_Leaf_N = self.IniLeafN_Conc * self.Initial_Leaf_Carbon / self.CarbonFrac_Veg
        self.Initial_Root_N = (self.Plant_Density * self.Seed_Weight * self.Germination_Efficiency * self.IniLeafN_Conc * self.CarbonAlloc_Shoot / self.NitrogenAlloc_Shoot) - self.Initial_Leaf_N
        self.Initial_Plant_Height = self.MaxPlant_Height / 1000.0  # Converting to a different unit if necessary
        self.MinLeafN_Conc = self.SLA_Const * self.Min_Specific_Leaf_N
        self.Initial_LAI = self.Initial_Leaf_Carbon / self.CarbonFrac_Veg * self.SLA_Const
        self.Initial_SLN_Bottom = self.Initial_Leaf_N / self.Initial_LAI


        # Actual canopy conditions
        self.Actual_AirCanopy_TemperatureDifference = 0
        self.Actual_Canopy_Temperature = 0
        self.Actual_Canopy_Photosynthesis = 0
        self.Actual_Canopy_Transpiration = 0
        self.Actual_Hourly_Canopy_Transpiration = []
    
        # Potential canopy conditions
        self.Potential_AirCanopy_TemperatureDifference = 0
        self.Potential_Canopy_Temperature = 0
        self.Potential_Canopy_Photosynthesis = 0
        self.Potential_Canopy_Transpiration = 0
        self.Canopy_PhotosyntheticallyActiveRadiation = 0
    
        self.Potential_Hourly_Canopy_Transpiration = []
    


        self.DevelopmentRate = 0
        self.NitrogenDemand_Allocation = 0
        self.Nitrogen_Demand = 0
        self.Shoot_Harvestable_Sugar_Allocation = 0
        self.End_SeedFill_Dynamic = 0
        self.Fraction_Seed_Carbon = 0
        self.Fraction_Stem_Carbon = 0
        self.Fraction_Leaf_Carbon = 0
        self.Fraction_Stem_CarbonReserve = 0
        self.Fraction_Root_CarbonReserve = 0
        self.Root_Carbon_ChangeRate = 0
        self.Stem_Carbon_ChangeRate = 0
        self.Seed_Carbon_ChangeRate = 0
        self.Leaf_Carbon_ChangeRate = 0
        self.Carbon_determined_LAI = self.Initial_LAI

    
        self.NitrogenFixation_Reserve_Pool_ChangeRate = 0
        self.Specific_Leaf_N_Bottom_ChangeRate = 0
    
        self.Nitrogen_DeadRoot = 0
    
        self.Maintenance_Respiration = 0
        self.RootWeight_Dead = 0
    
        self.Respired_Carbon = 0
        self.Carbon_Shoot = 0
        self.Root_Carbon = 0
        self.Nitrogen_Shoot = 0

        self.Model_TimeStep= Model_TimeStep


        self.DailyCarbonDemand_SeedFill_Cumulative = 0
        self.DailyCarbonDemand_Seed = 0
        self.DailyCarbon_Supply_Stem = 0
    
        self.FlowLimit_Seed_Carbon = 0
        self.FlowLimit_Stem_Carbon = 0
        self.DailyCarbonDemand_StemGrowth = 0
        self.DailyCarbon_Supply_Total = 0
        self.Plant_Height_ChangeRate = 0
        self.CarbonDemand_Stem_ChangeRate = 0
        self.Fraction_Carbon_Shoot = 0
    
        self.Remobilized_Carbon_Stem_to_Seed = 0
        self.Remobilized_Carbon_Root_to_Seed = 0
        self.Stem_CarbonReserve_ChangeRate = 0
        self.Root_CarbonReserve_ChangeRate = 0
    
        self.TotalPhotosynthesis_Canopy = 0
        self.LAI_ChangeRate = 0
    
        # New additions
        self.Carbon_DeadRoot_Total = 0
    
        self.DailyCarbonDemand_SeedFill_Current = 0
        self.DailyCarbonDemand_Seed_Current = 0
        self.DailyCarbon_Supply_Stem_Current = 0
    
        self.FlowLimit_CarbonSeed_Current = 0
        self.FlowLimit_CarbonStem_Current = 0
        self.DailyCarbonDemand_StemGrowth_Current = 0
        self.DailyCarbon_Supply_Total_Current = 0
        self.Plant_Height_ChangeRate_Current = 0

            
        
        
        # State Variables
        # Development stage and thermal units
        self.Development_Stage = 0
        self.Cumulative_Thermal_Unit = 0
    
        # Carbon in various plant components
        self.Leaf_Carbon = self.Initial_Leaf_Carbon
        self.DeadLeaf_Carbon = 0
        self.Stem_Carbon = 0
        self.Seed_Carbon = 0
        self.Root_Carbon = self.Initial_Root_Carbon
        self.DeadRoot_Carbon = 0
        self.DeadLeaf_Carbon_Litter = 0

        # Nitrogen in various plant components
        self.Nitrogen_Root = self.Initial_Root_N
        self.Nitrogen_Stem = 0
        self.Nitrogen_Leaf = self.Initial_Leaf_N
        self.Total_Leaf_Nitrogen = self.Initial_Leaf_N
        self.Nitrogen_Seed = 0
        self.Nitrogen_DeadLeaf = 0
    
        # Carbon reserves
        # self.Shoot_CarbonReserve = 0
        self.Root_CarbonReserve = 0
        self.Stem_CarbonReserve=0
    
        # Nitrogen dynamics
        self.Remobilizable_N_seedgrowth_after_seedfilling = 0
        self.Remobilizable_N_seedgrowth_before_seedfilling = 0 
        
        # Carbon dynamics
        self.CarbonDemand_Deficit_ForSeedFill_PreviousTimeSteps = 0
        self.CarbonDemand_Deficit_Stem_PreviousTimeSteps = 0
        self.Canopy_Nitrogen_Supply_PreviousTimeStep = 0
        # Additional state variables

        self.LAI_ChangeRate = 0

    
        # Resetting rates of change and flows
        self.Maintenance_Respiration = 0
        self.Respired_Carbon = 0
        self.FlowLimit_Seed_Carbon = 0
        self.FlowLimit_Stem_Carbon = 0

        self.Leaf_Nitrogen_Loss_ChangeRate = 0
        self.Leaf_Carbon_Loss_ChangeRate = 0
        self.NitrogenFixed_Reserve_pool = 0

        # Specific Leaf Nitrogen content and Leaf Area Index
        self.Specific_Leaf_N_Bottom = self.Initial_SLN_Bottom

    
        # Respiration MineralUptakes and Nitrogen dynamics
        self.Respiration_MineralUptakes = 0
        self.Canopy_Nitrogen_Demand_PreviousTimeStep = 0
        self.Nitrogen_Supply = 0
        self.Total_NitrogenFixed = 0
        self.NitrogenFixed_Reserve_Pool = 0
    
        # Carbon dynamics for stem growth and seed filling
        self.CarbonDemand_Stem_PreviousTimeSteps = 0
        self.Plant_Height = self.Initial_Plant_Height
        self.TotalRespiration = 0
        self.Cumulative_Respiration = 0
        self.Cumulative_LitterNitrogen = 0
    
        # Photosynthesis, Transpiration, and Canopy PAR
        self.Canopy_PhotosyntheticallyActiveRadiation = 0
        self.Actual_Canopy_Photosynthesis = 0
        self.Actual_Canopy_Transpiration = 0
        self.Cumulative_Canopy_PhotosyntheticallyActiveRadiation = 0
        self.Cumulative_Actual_Canopy_Photosynthesis = 0
        self.Cumulative_Actual_Canopy_Transpiration = 0
            
         
      
        
        
        
        

    
    def Leaf_to_Canopy_Integration(self,hourly_data_list1,hourly_data_list2):
        wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
        hourly_data_Canopy = np.array([sum(x) for x in zip(hourly_data_list1, hourly_data_list2)])
        weighted_data = hourly_data_Canopy * wgauss
        return weighted_data.sum()


    def Update_Canopy_PAR(self,hourly_apar_SH, hourly_apar_SU,dayl):
        daily_average_canopy_PAR =  self.Leaf_to_Canopy_Integration(hourly_apar_SH, hourly_apar_SU)
        self.Canopy_PhotosyntheticallyActiveRadiation=daily_average_canopy_PAR * dayl*  3600




    def Update_Canopy_Temp(self,State,hourly_Air_SU_leaf_T_diff,hourly_Air_SH_leaf_T_diff,
                           hourly_SU_leaf_T, hourly_SH_leaf_T):
        daily_average_temperature_dif=0.5* self.Leaf_to_Canopy_Integration(hourly_Air_SU_leaf_T_diff, hourly_Air_SH_leaf_T_diff)
        daily_average_temperature=0.5*self.Leaf_to_Canopy_Integration(hourly_SU_leaf_T, hourly_SH_leaf_T)

        State= State.lower()
        if State=='p':
            self.Potential_Air_canpopy_Tdif=daily_average_temperature_dif
            self.Potential_Canopy_Temp=daily_average_temperature
        elif State=='a':
            self.Actual_Air_canpopy_Tdif=daily_average_temperature_dif
            self.Actual_Canopy_Temp=daily_average_temperature
        else:
            raise ValueError("Invalid State provided. State must be 'A' for Actual or 'P' for Potential.")
            
            
    def Update_Canopy_Photosyn(self,State,photosyn_SU, photosyn_SH,dayl):
        print(photosyn_SU, photosyn_SH)
        daily_average_canopy_photosyn =  self.Leaf_to_Canopy_Integration(photosyn_SU, photosyn_SH)
        State= State.lower()
        if State=='p':
            self.Potential_Canopy_Photosynthesis=daily_average_canopy_photosyn * dayl*  3600
        elif State=='a':
            self.Actual_Canopy_Photosynthesis=daily_average_canopy_photosyn * dayl*  3600
        else:
            raise ValueError("Invalid State provided. State must be 'A' for Actual or 'P' for Potential.")


    def Update_Canopy_Transpiration(self,State,transpiration_SU, transpiration_SH,dayl):
        daily_average_canopy_transpiration =  self.Leaf_to_Canopy_Integration(transpiration_SU, transpiration_SH)
        State= State.lower()

        if State=='p':
            self.Potential_Canopy_Transpiration=daily_average_canopy_transpiration * dayl*  3600
        elif State=='a':
            self.Actual_Canopy_Transpiration=daily_average_canopy_transpiration * dayl*  3600
        else:
            raise ValueError("Invalid State provided. State must be 'A' for Actual or 'P' for Potential.")

        
   

    def calculate_thermal_units(self, Development_Stage, Tmax, Tmin, DayLength, BTemp_Phen, OTemp_Phen, CTemp_Phen, TempCurve_Res):
        
        # Timing for sunrise and sunset
        SunriseTime = 12. - 0.5 * DayLength
        SunsetTime = 12. + 0.5 * DayLength
    
        # Mean daily temperature
        Tmean = (Tmax + Tmin) / 2
        TotalThermalUnits = 0  # Total thermal units
    
        # Diurnal course of temperature
        for Hour in range(1, 25):
            if SunriseTime <= Hour <= SunsetTime:
                TempHourly = Tmean + 0.5 * abs(Tmax - Tmin) * math.cos(0.2618 * (Hour - 14))
            else:
                TempHourly = Tmean + 0.5 * abs(Tmax - Tmin) * math.cos(0.2618 * (Hour - 14))
    
            # Assuming development rate at supra-optimum temperatures during the reproductive
            # phase equals that at the optimum temperature
            if Development_Stage > 1:
                TempHourly = min(TempHourly, OTemp_Phen)
    
            # Instantaneous thermal unit based on bell-shaped temperature response
            if TempHourly < BTemp_Phen or TempHourly > CTemp_Phen:
                ThermalUnit = 0
            else:
                ThermalUnit = (((CTemp_Phen - TempHourly) / (CTemp_Phen - OTemp_Phen)) * ((TempHourly - BTemp_Phen) / (OTemp_Phen - BTemp_Phen)) ** ((OTemp_Phen - BTemp_Phen) / (CTemp_Phen - OTemp_Phen))) ** TempCurve_Res
            TotalThermalUnits += ThermalUnit / 24
            
        # Daily thermal unit
        self.Daily_Thermal_Unit = TotalThermalUnits
    
        
    

        
    def calculate_developement_rate(self, Development_Stage, CropType_Photoperiod, DayLength_PhotoPeriod, StartPhotoperiod_Phase, EndPhotoperiod_Phase, Photoperiod_Sensitivity, MinThermal_Day_Veg, MinThermal_Day_Rep, Daily_Thermal_Unit):
        # Determining if it is for short-day or long-day crop
        if CropType_Photoperiod < 0:
            OptimumPhotoperiod = 18  # Minimum optimum photoperiod for long-day crop
            DayLengthResponse = min(OptimumPhotoperiod, DayLength_PhotoPeriod)
        else:
            OptimumPhotoperiod = 11  # Maximum optimum photoperiod for short-day crop
            DayLengthResponse = max(OptimumPhotoperiod, DayLength_PhotoPeriod)
    
        # Effect of photoperiod on development rate
        if Development_Stage < StartPhotoperiod_Phase or Development_Stage > EndPhotoperiod_Phase:
            EffectPhotoperiod = 1
        else:
            EffectPhotoperiod = max(0, 1 - Photoperiod_Sensitivity * (DayLengthResponse - OptimumPhotoperiod))
    
        # Development rate of vegetative and reproductive phases
        if 0 <= Development_Stage < 1.0:
            DevelopmentRate = 1 / MinThermal_Day_Veg * Daily_Thermal_Unit * EffectPhotoperiod
        else:
            DevelopmentRate = 1 / MinThermal_Day_Rep * Daily_Thermal_Unit
    
        self.DevelopmentRate = DevelopmentRate





    def Initialize_Biomass_Formation(self):
        Seed_Dry_Weight = self.Seed_Carbon / self.CarbonFraction_Organs # Dry weight of seed
        LiveLeaf_Dry_Weight = self.Leaf_Carbon / self.CarbonFrac_Veg # Dry weight of live leaves
        Stem_Dry_Weight = self.Stem_Carbon / self.CarbonFrac_Veg + self.Stem_CarbonReserve / 0.444 # Dry weight of stems
        LiveRoot_Dry_Weight = self.Root_Carbon / self.CarbonFrac_Veg + self.Root_CarbonReserve / 0.444 # Dry weight of live roots
    
        Shoot_Dry_Weight = LiveLeaf_Dry_Weight + Stem_Dry_Weight + Seed_Dry_Weight # Dry weight of live shoot (above-ground) organs
    
        TotalLiveOrgan_Dry_Weight = Shoot_Dry_Weight + LiveRoot_Dry_Weight # Dry weight of total live organs
    
        DeadLeaf_Dry_Weight = self.DeadLeaf_Carbon / self.CarbonFrac_Veg # Dry weight of dead leaves
    
        Shoot_Dry_Weight_ExcShedLeaves = Shoot_Dry_Weight + (DeadLeaf_Dry_Weight - self.DeadLeaf_Carbon_Litter / self.CarbonFrac_Veg) # Dry weight of shoot organs (excluding shedded leaves)
    
        Harvest_Index = Seed_Dry_Weight / Shoot_Dry_Weight_ExcShedLeaves if Shoot_Dry_Weight_ExcShedLeaves else 0  # Avoid division by zero
        DeadRoot_Dry_Weight = self.DeadRoot_Carbon / self.CarbonFrac_Veg # Dry weight of dead roots
        self.TotalLiveOrgan_Dry_Weight = TotalLiveOrgan_Dry_Weight
        self.Harvest_Index = Harvest_Index
        self.LiveLeaf_Dry_Weight = LiveLeaf_Dry_Weight
        self.Shoot_Dry_Weight = Shoot_Dry_Weight
        self.LiveRoot_Dry_Weight = LiveRoot_Dry_Weight
        self.Seed_Dry_Weight = Seed_Dry_Weight  
        self.DeadRoot_Dry_Weight = DeadRoot_Dry_Weight



    def Initialize_Nitrogen_Accumulation(self):
        Nitrogen_Shoot = self.Nitrogen_Stem + self.Nitrogen_Leaf + self.Nitrogen_Seed
        LiveLeaf_Dry_Weight = self.DeadLeaf_Carbon / self.CarbonFrac_Veg
        Nitrogen_Shoot_ExcShedLeaves = Nitrogen_Shoot + (LiveLeaf_Dry_Weight - self.DeadLeaf_Carbon_Litter / self.CarbonFrac_Veg) * self.MinLeafN_Conc
        Nitrogen_Total = Nitrogen_Shoot + self.Nitrogen_Root
        Leaf_Nitrogen_Conc = self.Nitrogen_Leaf / self.LiveLeaf_Dry_Weight if self.LiveLeaf_Dry_Weight else 0  # Nitrogen concentration in living leaves
    
        Shoot_Nitrogen_Conc = (self.Nitrogen_Stem + self.Nitrogen_Leaf + self.Nitrogen_Seed) / self.Shoot_Dry_Weight if self.Shoot_Dry_Weight else 0 # Nitrogen concentration in living shoot
    
        Root_Nitrogen_Conc = self.Nitrogen_Root / self.LiveRoot_Dry_Weight if self.LiveRoot_Dry_Weight else 0 # Nitrogen concentration in living roots
    
        Seed_Nitrogen_Conc = self.Switch_Function(self.Seed_Dry_Weight * -1, self.Nitrogen_Seed / max(self.Avoid_Zero_Division(self.Seed_Dry_Weight), 0), 0) # Nitrogen concentration in the storage organs 
        Plant_Nitrogen_Conc = Nitrogen_Total / self.TotalLiveOrgan_Dry_Weight if self.TotalLiveOrgan_Dry_Weight else 0 # Nitrogen concentration in living plant material
    
        
        self.Nitrogen_Shoot = Nitrogen_Shoot
        self.Leaf_Nitrogen_Conc = Leaf_Nitrogen_Conc
        self.Shoot_Nitrogen_Conc = Shoot_Nitrogen_Conc
        self.Root_Nitrogen_Conc = Root_Nitrogen_Conc
        self.Seed_Nitrogen_Conc = Seed_Nitrogen_Conc
        self.Plant_Nitrogen_Conc = Plant_Nitrogen_Conc
        self.Nitrogen_Total = Nitrogen_Total



    def Initialize_Carbon_Accumulation(self):
     Carbon_Shoot = self.Leaf_Carbon + self.Stem_Carbon + self.Stem_CarbonReserve + self.Seed_Carbon
     Root_Carbon = self.Root_Carbon + self.Root_CarbonReserve
     Carbon_Total = Carbon_Shoot + Root_Carbon
     self.Carbon_Total = Carbon_Total
     self.Carbon_Shoot = Carbon_Shoot
     self.Root_Carbon = Root_Carbon



    def Calculate_Respiration(self):
        Maintenance_Respiration_Extra = max(min(44./12.*0.218*(self.Nitrogen_Total-self.Shoot_Dry_Weight*self.MinLeafN_Conc-self.LiveRoot_Dry_Weight*self.MinRootN_Conc), self.Actual_Canopy_Photosynthesis-1.E-5-self.Maintenance_Respiration), 0.)
        Maintenance_Respiration = max(0., min(self.Actual_Canopy_Photosynthesis-1.E-5, self.Respiration_MineralUptakes) + Maintenance_Respiration_Extra)
        self.Maintenance_Respiration = Maintenance_Respiration

    def Calculate_Nitrogen_Fixation(self):
        Nitrogen_Fixed_Demand = max(0., self.Canopy_Nitrogen_Demand_PreviousTimeStep - self.Canopy_Nitrogen_Supply_PreviousTimeStep)
        Nitrogen_Fixed_Efficiency = max(0., self.Actual_Canopy_Photosynthesis - 1.E-5 - self.Maintenance_Respiration) / self.CarbonCost_NFix * 12. / 44.
        Nitrogen_Fixed = self.Switch_Function(self.legume, 0., min(Nitrogen_Fixed_Efficiency, Nitrogen_Fixed_Demand))
        self.Nitrogen_Fixed = Nitrogen_Fixed

    def Calculate_Photo_Assimilates(self):
        Respired_Carbon = 44. / 12. * (self.CarbonCost_NFix * self.Nitrogen_Fixed)
        Photo_Assimilate = self.Actual_Canopy_Photosynthesis - self.Maintenance_Respiration - Respired_Carbon
        self.Photo_Assimilate = Photo_Assimilate
        self.Respired_Carbon = Respired_Carbon

    


    def Calculate_Crop_Nitrogen_Demand(self, Specific_Leaf_Nitrogen):
        Shoot_Harvestable_Sugar_Allocation = 12.0 / 44.0 * self.Growth_Efficiency_Veg * (self.Actual_Canopy_Photosynthesis - self.Maintenance_Respiration - self.Respired_Carbon) / self.Carbon_Shoot
        Nitrogen_Demand_Allocation = self.Root_Carbon * Shoot_Harvestable_Sugar_Allocation**2 
        High_Nitrogen_Conc_Crit = self.IniLeafN_Conc * np.exp(-0.4 * self.Development_Stage)
        Nitrogen_Demand_Dynamic = self.Switch_Function(self.Development_Stage - 1.0, self.Shoot_Dry_Weight * (High_Nitrogen_Conc_Crit - self.Shoot_Nitrogen_Conc) * (1.0 + self.Nitrogen_Root / self.Nitrogen_Shoot) / self.Model_TimeStep, 0.0)
        Nitrogen_Demand_Adjusted = self.Switch_Function(self.Leaf_Nitrogen_Conc - 1.5 * self.IniLeafN_Conc, max(Nitrogen_Demand_Allocation, Nitrogen_Demand_Dynamic), 0)
        Nitrogen_Demand = self.Switch_Function(self.Min_Specific_Leaf_N - Specific_Leaf_Nitrogen + 1.0e-5, min(self.MaxN_Uptake, Nitrogen_Demand_Adjusted), 0)  # Crop nitrogen demand
    
        print(self.Root_Carbon, Shoot_Harvestable_Sugar_Allocation)
        self.Nitrogen_Demand_Allocation = Nitrogen_Demand_Allocation
        self.Nitrogen_Demand = Nitrogen_Demand
    



    def Calculate_Nitrogen_Partitioning(self, Specific_Leaf_N_Top_Increment):
        print("*******************************************************************")
        print((self.Actual_Canopy_Photosynthesis))
        Nitrogen_Carbon_Ratio = self.Switch_Function(Specific_Leaf_N_Top_Increment - self.Min_Specific_Leaf_N, 0, min(self.MaxN_Uptake, self.Nitrogen_Demand_Allocation)) / (self.Growth_Efficiency_Veg * (self.Actual_Canopy_Photosynthesis - self.Maintenance_Respiration - self.Respired_Carbon) * 12 / 44)
        Shoot_Harvestable_Sugar_Allocation = 12.0 / 44.0 * self.Growth_Efficiency_Veg * (self.Actual_Canopy_Photosynthesis - self.Maintenance_Respiration - self.Respired_Carbon) / self.Carbon_Shoot
        Fraction_Shoot_Nitrogen = 1 / (1 + Nitrogen_Carbon_Ratio * 1 / Shoot_Harvestable_Sugar_Allocation * self.Carbon_Shoot / self.Root_Carbon * self.Nitrogen_Root / self.Nitrogen_Shoot)
        
        print(Fraction_Shoot_Nitrogen)
        self.Shoot_Harvestable_Sugar_Allocation = Shoot_Harvestable_Sugar_Allocation
        self.Fraction_Shoot_Nitrogen = Fraction_Shoot_Nitrogen
        self.Nitrogen_Carbon_Ratio = Nitrogen_Carbon_Ratio


    def Calculate_Seed_Properties(self, Crop_TypeDet, EndSeedNum_DetPeriod, Remobilizable_N_seedgrowth_before_seedfilling, Remobilizable_N_seedgrowth_after_seedfilling, Fraction_SeedNitrogen_Remobilizable, Standard_SeedNitrogen_Conc, Seed_Weight):
        End_SeedFill_Dynamic = self.Switch_Function(Crop_TypeDet, EndSeedNum_DetPeriod, 1.)
        Nitrogen_Remobilization = Remobilizable_N_seedgrowth_before_seedfilling + (Remobilizable_N_seedgrowth_after_seedfilling - Remobilizable_N_seedgrowth_before_seedfilling) * (End_SeedFill_Dynamic - 1.0) / self.Avoid_Zero_Division(min(self.Development_Stage, End_SeedFill_Dynamic) - 1)
        Total_Seed_Nitrogen = Nitrogen_Remobilization / Fraction_SeedNitrogen_Remobilizable / Standard_SeedNitrogen_Conc / Seed_Weight
        Total_Seed_Weight = self.Seed_Carbon / self.CarbonFraction_Organs / self.Avoid_Zero_Division(Total_Seed_Nitrogen) * 1000
        
        self.End_SeedFill_Dynamic = End_SeedFill_Dynamic
        self.Total_Seed_Nitrogen = Total_Seed_Nitrogen
        self.Total_Seed_Weight = Total_Seed_Weight

    
    def Calculate_Senescence(self, Carbon_determined_LAI, Nitrogen_determined_LAI, Time_Step):
        LeafWeight_Loss_Max = (Carbon_determined_LAI - min(Carbon_determined_LAI, Nitrogen_determined_LAI)) / self.SLA_Const / Time_Step
        LeafWeight_Loss = min(self.LiveLeaf_Dry_Weight - 1.e-5, LeafWeight_Loss_Max + self.REANOR_Function(self.End_SeedFill_Dynamic - self.Development_Stage, LeafWeight_Loss_Max) * 0.03 * self.LiveLeaf_Dry_Weight)
        Leaf_Nitrogen_Loss_ChangeRate = min(LeafWeight_Loss, LeafWeight_Loss_Max) * self.MinLeafN_Conc + (LeafWeight_Loss - min(LeafWeight_Loss, LeafWeight_Loss_Max)) * (self.Nitrogen_Leaf / self.LiveLeaf_Dry_Weight)
        Leaf_Carbon_Loss_ChangeRate = LeafWeight_Loss_Max * self.CarbonFrac_Veg
        
        print(Leaf_Nitrogen_Loss_ChangeRate)
        self.Leaf_Nitrogen_Loss_ChangeRate = Leaf_Nitrogen_Loss_ChangeRate
        self.Leaf_Carbon_Loss_ChangeRate = Leaf_Carbon_Loss_ChangeRate




    def Calculate_Carbon_Flow(self):
        Fraction_Shoot = 1 / (1 + self.Nitrogen_Carbon_Ratio * 1 / self.Shoot_Harvestable_Sugar_Allocation)  
        DailyCarbon_Supply_Shoot = 12.0 / 44.0 * Fraction_Shoot * self.Photo_Assimilate  # Daily carbon supply for shoot growth from assimilates
        DailyCarbon_Supply_Root = 12 / 44 * (1 - Fraction_Shoot) * self.Photo_Assimilate  # Daily carbon supply for root growth from assimilates
        
        self.DailyCarbon_Supply_Shoot = DailyCarbon_Supply_Shoot
        self.DailyCarbon_Supply_Root = DailyCarbon_Supply_Root
        self.Fraction_Shoot = Fraction_Shoot
    



    def Calculate_Carbon_Flow_Seed_Filling(self):
        SeedFill_Transition = self.MaxSeedGrowth_DS * 1 # Transition period
        SeedFill_Phase = 1 # Seed filling phase length
        SeedFill_Start = 1 # Start of seed filling growth
        SeedFill_Day = self.Limit_Function(1, 2, self.Development_Stage) - 1 # Day of seed filling
    
        DailySeed_Fill_Rate = self.DevelopmentRate * (2 * SeedFill_Phase - SeedFill_Transition) * (SeedFill_Phase - SeedFill_Day) / SeedFill_Phase / (SeedFill_Phase - SeedFill_Transition) ** 2 * (SeedFill_Day / SeedFill_Phase) ** (SeedFill_Transition / (SeedFill_Phase - SeedFill_Transition))
    
        TotalSeed_Carbons = self.Total_Seed_Nitrogen * self.Seed_Weight * self.CarbonFraction_Organs
        DailyCarbonDemand_SeedFill = self.Switch_Function(self.Development_Stage-SeedFill_Start, 0., TotalSeed_Carbons/self.Growth_Efficiency_Veg*DailySeed_Fill_Rate) 
        DailySeed_Carbon_ChangeRate = DailyCarbonDemand_SeedFill + max(0, self.CarbonDemand_Deficit_ForSeedFill_PreviousTimeSteps) / self.Model_TimeStep
        FlowLimit_Seed_Carbon = min(DailySeed_Carbon_ChangeRate, self.DailyCarbon_Supply_Shoot)
        
        self.DailyCarbonDemand_SeedFill = DailyCarbonDemand_SeedFill
        self.DailySeed_Carbon_ChangeRate = DailySeed_Carbon_ChangeRate
        self.FlowLimit_Seed_Carbon = FlowLimit_Seed_Carbon



    def Calculate_Carbon_Flow_Stem_Growth(self):
        DailyCarbon_Supply_Stem = self.DailyCarbon_Supply_Shoot - self.FlowLimit_Seed_Carbon  # Daily carbon supply from current photosynthesis for structural stem growth
    
        IntegralFactor_Stress_Height = self.Limit_Function(0, 1, DailyCarbon_Supply_Stem / self.Avoid_Zero_Division(self.CarbonDemand_Stem_PreviousTimeSteps))  # Integral factor of stresses on plant height growth
        Transition_Height = self.MaxStemGrowth_DS * (1. + self.End_SeedFill_Dynamic) / 2. # Transition period for height growth
        Phase_Height = (1. + self.End_SeedFill_Dynamic) / 2.  # Phase length for height growth
        SeedFill_Day = min((1. + self.End_SeedFill_Dynamic) / 2., self.Development_Stage)  # Day of height growth phase
    
        DailyHeight_Growth_Rate = self.DevelopmentRate * (2 * Phase_Height - Transition_Height) * (Phase_Height - SeedFill_Day) / Phase_Height / (Phase_Height - Transition_Height) ** 2 * (SeedFill_Day / Phase_Height) ** (Transition_Height / (Phase_Height - Transition_Height))
    
        Start_Growth = 0 # Start of structural stem growth
        TotalStem_Carbon = self.StemDW_Height * self.MaxPlant_Height * self.CarbonFrac_Veg
        DailyCarbonDemand_StemGrowth = self.Switch_Function(self.Development_Stage - Start_Growth, 0., TotalStem_Carbon / self.Growth_Efficiency_Veg * DailyHeight_Growth_Rate * IntegralFactor_Stress_Height)
        DailyCarbon_Supply_Total = DailyCarbonDemand_StemGrowth + max(0, self.CarbonDemand_Deficit_Stem_PreviousTimeSteps) / self.Model_TimeStep
        FlowLimit_Stem_Carbon = min(DailyCarbon_Supply_Total, DailyCarbon_Supply_Stem)
    
        # Daily carbon flow for structural stem growth
        Plant_Height_ChangeRate = min(self.MaxPlant_Height - self.Plant_Height, DailyHeight_Growth_Rate * self.MaxPlant_Height * IntegralFactor_Stress_Height)  # Rate of plant height growth
        CarbonDemand_Stem_ChangeRate = (DailyCarbonDemand_StemGrowth - self.CarbonDemand_Stem_PreviousTimeSteps) / self.Model_TimeStep  # Carbon demand for structural stem growth at the previous time step
        
        self.DailyCarbon_Supply_Stem = DailyCarbon_Supply_Stem
        self.FlowLimit_Stem_Carbon = FlowLimit_Stem_Carbon
        self.DailyCarbonDemand_StemGrowth = DailyCarbonDemand_StemGrowth
        self.DailyCarbon_Supply_Total = DailyCarbon_Supply_Total
        self.Plant_Height_ChangeRate = Plant_Height_ChangeRate
        self.CarbonDemand_Stem_ChangeRate = CarbonDemand_Stem_ChangeRate
    




    
    def Calculate_Carbon_Partitioning(self, Nitrogen_determined_LAI, Carbon_determined_LAI, Root_Carbon_ChangeRateNew):
        Fraction_Seed_Carbon = self.FlowLimit_Seed_Carbon / self.DailyCarbon_Supply_Shoot
        Fraction_Stem_Carbon = self.Switch_Function(self.Development_Stage - (self.End_SeedFill_Dynamic + 0.2), self.FlowLimit_Stem_Carbon / self.DailyCarbon_Supply_Shoot, 0)
        Fraction_Leaf_Carbon = self.Indicator_Function(Nitrogen_determined_LAI - Carbon_determined_LAI, self.End_SeedFill_Dynamic - self.Development_Stage) * (1.0 - Fraction_Seed_Carbon - Fraction_Stem_Carbon)
        Fraction_Stem_CarbonReserve = 1.0 - Fraction_Leaf_Carbon - Fraction_Seed_Carbon - Fraction_Stem_Carbon  # Fraction of new shoot carbon to stem reserves
        Fraction_Root_CarbonReserve = self.Switch_Function(Root_Carbon_ChangeRateNew - self.Root_Carbon, 1.0, 0.0)  # Fraction of new root carbon to root reserves
    
        self.Fraction_Seed_Carbon = Fraction_Seed_Carbon
        self.Fraction_Stem_Carbon = Fraction_Stem_Carbon
        self.Fraction_Leaf_Carbon = Fraction_Leaf_Carbon
        self.Fraction_Stem_CarbonReserve = Fraction_Stem_CarbonReserve
        self.Fraction_Root_CarbonReserve = Fraction_Root_CarbonReserve

    def Update_Carbon_Reserve_Pools(self):
        Gap_Carbon = max(0, self.DailyCarbonDemand_Seed - self.DailyCarbon_Supply_Shoot)
        Intermediate_variable_Stem = min(0.94 * self.Stem_CarbonReserve, self.Stem_CarbonReserve / self.Avoid_Zero_Division(self.Stem_CarbonReserve + self.Root_CarbonReserve) * Gap_Carbon) / 0.94
        Intermediate_variable_Root = min(0.94 * self.Root_CarbonReserve, self.Root_CarbonReserve / self.Avoid_Zero_Division(self.Stem_CarbonReserve + self.Root_CarbonReserve) * Gap_Carbon) / 0.94
        Remobilized_Carbon_Stem_to_Seed = self.Switch_Function(self.DailyCarbonDemand_Seed - self.DailyCarbon_Supply_Shoot, 0, Intermediate_variable_Stem)
        
        Remobilized_Carbon_Root_to_Seed = self.Switch_Function(self.DailyCarbonDemand_Seed - self.DailyCarbon_Supply_Shoot, 0, Intermediate_variable_Root)
        
        Stem_CarbonReserve_ChangeRate = self.Fraction_Stem_CarbonReserve * self.DailyCarbon_Supply_Shoot - Remobilized_Carbon_Stem_to_Seed
        Root_CarbonReserve_ChangeRate = self.Fraction_Root_CarbonReserve * self.DailyCarbon_Supply_Root - Remobilized_Carbon_Root_to_Seed
    
        self.Remobilized_Carbon_Stem_to_Seed = Remobilized_Carbon_Stem_to_Seed
        self.Remobilized_Carbon_Root_to_Seed = Remobilized_Carbon_Root_to_Seed
        self.Stem_CarbonReserve_ChangeRate = Stem_CarbonReserve_ChangeRate
        self.Root_CarbonReserve_ChangeRate = Root_CarbonReserve_ChangeRate

        
    def Calculate_Carbon_Production_Rate(self, Root_carbon_loss_rate_senescence):
        Root_Carbon_ChangeRate = 12 / 44 * self.Photo_Assimilate * (1 - self.Fraction_Shoot) * (1 - self.Fraction_Root_CarbonReserve) * self.Growth_Efficiency_Veg - Root_carbon_loss_rate_senescence
        Stem_Carbon_ChangeRate = 12.0 / 44.0 * self.Photo_Assimilate * self.Fraction_Shoot * self.Fraction_Stem_Carbon * self.Growth_Efficiency_Veg
        Seed_Carbon_ChangeRate = 12.0 / 44.0 * self.Photo_Assimilate * self.Fraction_Shoot * self.Fraction_Seed_Carbon * self.Growth_Efficiency_Veg + 0.94 * (self.Remobilized_Carbon_Stem_to_Seed + self.Remobilized_Carbon_Root_to_Seed) * self.Growth_Efficiency_Veg
        Leaf_Carbon_ChangeRate = 12.0 / 44.0 * self.Photo_Assimilate * self.Fraction_Shoot * self.Fraction_Leaf_Carbon * self.Growth_Efficiency_Veg - self.Leaf_Carbon_Loss_ChangeRate
    
        self.Root_Carbon_ChangeRate = Root_Carbon_ChangeRate
        self.Stem_Carbon_ChangeRate = Stem_Carbon_ChangeRate
        self.Seed_Carbon_ChangeRate = Seed_Carbon_ChangeRate
        self.Leaf_Carbon_ChangeRate = Leaf_Carbon_ChangeRate


    def Calculate_Biomass_Change_Rate(self):
        RootWeight_Rate = self.Root_Carbon_ChangeRate / self.CarbonFrac_Veg + self.Root_CarbonReserve_ChangeRate / 0.444
        SeedWeight_Rate = self.Seed_Carbon_ChangeRate / self.CarbonFraction_Organs
        LeafWeight_Rate = self.Leaf_Carbon_ChangeRate / self.CarbonFrac_Veg
        StemWeight_Rate = self.Stem_Carbon_ChangeRate / self.CarbonFrac_Veg + self.Stem_CarbonReserve_ChangeRate / 0.444
    
        self.RootWeight_Rate = RootWeight_Rate
        self.SeedWeight_Rate = SeedWeight_Rate
        self.LeafWeight_Rate = LeafWeight_Rate
        self.StemWeight_Rate = StemWeight_Rate

    
    def Calculate_Carbon_Rates(self):
        CarbonDemand_Deficit_ForSeedFill_ChangeRate = max(0.0, (self.DailyCarbonDemand_SeedFill - self.Seed_Carbon_ChangeRate / self.Growth_Efficiency_Veg)) - (self.FlowLimit_Seed_Carbon - min(self.DailyCarbonDemand_SeedFill, self.DailyCarbon_Supply_Shoot))
        CarbonDemand_Deficit_Stem_ChangeRate = max(0., (self.DailyCarbonDemand_StemGrowth - self.Stem_Carbon_ChangeRate / self.Growth_Efficiency_Veg)) - (self.FlowLimit_Stem_Carbon - min(self.DailyCarbonDemand_StemGrowth, self.DailyCarbon_Supply_Stem))
    
        self.CarbonDemand_Deficit_ForSeedFill_ChangeRate = CarbonDemand_Deficit_ForSeedFill_ChangeRate
        self.CarbonDemand_Deficit_Stem_ChangeRate = CarbonDemand_Deficit_Stem_ChangeRate
    


    def Calculate_Nitrogen_Dynamics(self, NitrogenUptake):
        NitrogenFixation_Reserve_Pool_ChangeRate = self.Nitrogen_Fixed - min(self.Nitrogen_Demand, self.NitrogenFixed_Reserve_pool )
        NitrogenDemand_ChangeRate = (self.Nitrogen_Demand - self.Canopy_Nitrogen_Demand_PreviousTimeStep) / self.Model_TimeStep
        NitrogenSupply_ChangeRate = (NitrogenUptake - self.Canopy_Nitrogen_Supply_PreviousTimeStep) / self.Model_TimeStep
    
        self.NitrogenFixation_Reserve_Pool_ChangeRate = NitrogenFixation_Reserve_Pool_ChangeRate
        self.NitrogenDemand_ChangeRate = NitrogenDemand_ChangeRate
        self.NitrogenSupply_ChangeRate = NitrogenSupply_ChangeRate



    def Calculate_Seed_Number_Rate(self, NitrogenUptake, Root_carbon_loss_rate_senescence):
        NitrogenReserve_Rate = NitrogenUptake - (self.MinLeafN_Conc * (self.Leaf_Carbon_ChangeRate + self.Leaf_Carbon_Loss_ChangeRate) + self.MinRootN_Conc * (self.Root_Carbon_ChangeRate + Root_carbon_loss_rate_senescence) + self.Stem_Nitrogen * self.Stem_Carbon_ChangeRate) / self.CarbonFrac_Veg
        Remobilizable_N_seedgrowth_after_seedfilling_Rate = self.Switch_Function(self.Development_Stage - self.End_SeedFill_Dynamic, NitrogenReserve_Rate, 0.0)
        Remobilizable_N_seedgrowth_before_seedfilling_Rate = self.Switch_Function(self.Development_Stage - 1.0, NitrogenReserve_Rate, 0.0)
    
        self.Remobilizable_N_seedgrowth_after_seedfilling_Rate = Remobilizable_N_seedgrowth_after_seedfilling_Rate
        self.Remobilizable_N_seedgrowth_before_seedfilling_Rate = Remobilizable_N_seedgrowth_before_seedfilling_Rate

    
    def Calculate_Nitrogen_Accumulation_Rate(self, NitrogenUptake, Root_nitrogen_loss_rate_senescence,Standard_SeedNitrogen_Conc):
        # Amount of N partitioned to shoot
        Nitrogen_Shoot_New = self.Fraction_Shoot_Nitrogen * NitrogenUptake
    
        # Leaf N or root N available for remobilization
        Nitrogen_Leaf_Available = self.Switch_Function(self.MinLeafN_Conc - self.Leaf_Nitrogen_Conc, self.Nitrogen_Leaf - self.LiveLeaf_Dry_Weight * self.MinLeafN_Conc, 0.) / self.Model_TimeStep
        Nitrogen_Root_Available = self.Switch_Function(self.MinRootN_Conc - self.Root_Nitrogen_Conc, self.Nitrogen_Root - self.LiveRoot_Dry_Weight * self.MinRootN_Conc, 0.) / self.Model_TimeStep
        Nitrogen_Total_Available = Nitrogen_Leaf_Available + Nitrogen_Root_Available
    
        # Rate of N accumulation in stem
        Nitrogen_Stem_ChangeRate = self.StemWeight_Rate * self.Switch_Function(-Nitrogen_Total_Available, self.Stem_Nitrogen, 0.)
    
        # Expected N dynamics during seed filling
        Expected_SeedNitrogen_Conc =  Standard_SeedNitrogen_Conc
    
        # Rate of N accumulation in seed
        Nitrogen_Seed_Growth = Nitrogen_Shoot_New - Nitrogen_Stem_ChangeRate - Expected_SeedNitrogen_Conc * self.SeedWeight_Rate
        Nitrogen_Organic_Compound_Conc = max(0., self.Switch_Function(Nitrogen_Total_Available + Nitrogen_Seed_Growth, (Nitrogen_Total_Available + Nitrogen_Shoot_New - Nitrogen_Stem_ChangeRate) / self.Avoid_Zero_Division(self.SeedWeight_Rate), Expected_SeedNitrogen_Conc))
        Nitrogen_Seed_ChangeRate = self.SeedWeight_Rate * Nitrogen_Organic_Compound_Conc
    
        # Rate of N accumulation in leaf
        Nitrogen_Leaf_Accumulation = self.Switch_Function(Nitrogen_Total_Available + Nitrogen_Seed_Growth, -Nitrogen_Leaf_Available - self.Leaf_Nitrogen_Loss_ChangeRate, -Nitrogen_Leaf_Available / self.Avoid_Zero_Division(Nitrogen_Total_Available) * (-Nitrogen_Seed_Growth) - self.Leaf_Nitrogen_Loss_ChangeRate)
        LeafNitrogen_Growth = self.Switch_Function(Nitrogen_Seed_Growth, Nitrogen_Leaf_Accumulation, Nitrogen_Shoot_New - Nitrogen_Stem_ChangeRate - Nitrogen_Seed_ChangeRate - self.Leaf_Nitrogen_Loss_ChangeRate)
        Nitrogen_Leaf_ChangeRate = max(-self.Nitrogen_Leaf + 1E-7, LeafNitrogen_Growth)
        Nitrogen_Leaf_ChangeRate_Positive = max(0, Nitrogen_Leaf_ChangeRate)
    
        # Rate of N accumulation in root
        Nitrogen_Root_Accumulation = self.Switch_Function(Nitrogen_Total_Available + Nitrogen_Seed_Growth, NitrogenUptake - Nitrogen_Shoot_New - Nitrogen_Root_Available - Root_nitrogen_loss_rate_senescence, NitrogenUptake - Nitrogen_Shoot_New - Nitrogen_Root_Available / self.Avoid_Zero_Division(Nitrogen_Total_Available) * (-Nitrogen_Seed_Growth) - Root_nitrogen_loss_rate_senescence)
        RootNitrogen_Growth = self.Switch_Function(Nitrogen_Seed_Growth, Nitrogen_Root_Accumulation, NitrogenUptake - Nitrogen_Shoot_New - Root_nitrogen_loss_rate_senescence)
        Nitrogen_Root_ChangeRate = max(-self.Nitrogen_Root + 5E-8, RootNitrogen_Growth)
    
        self.Nitrogen_Root_ChangeRate = Nitrogen_Root_ChangeRate
        self.Nitrogen_Stem_ChangeRate = Nitrogen_Stem_ChangeRate
        self.Nitrogen_Leaf_ChangeRate = Nitrogen_Leaf_ChangeRate
        self.Nitrogen_Leaf_ChangeRate_Positive = Nitrogen_Leaf_ChangeRate_Positive
        self.Nitrogen_Seed_ChangeRate = Nitrogen_Seed_ChangeRate

    def Calculate_Leaf_Area_ChangeRate(self, Carbon_determined_LAI, Specific_Leaf_N_Bottom_Exponential_with_Depth, LeafNitrogen_ExtinctionCoefficient):
        Specific_Leaf_N_Bottom_ChangeRate = (Specific_Leaf_N_Bottom_Exponential_with_Depth - self.Specific_Leaf_N_Bottom) / self.Model_TimeStep
    
        # Rate of LAI driven by carbon supply
        LAI_ChangeRate = self.Switch_Function(self.LeafWeight_Rate, max(-Carbon_determined_LAI + 1E-5, self.SLA_Const * self.LeafWeight_Rate), self.SLA_Const * self.LeafWeight_Rate)
    
        # Adjusting LAI rate based on nitrogen during the juvenile phase
        if Carbon_determined_LAI < 1 and self.Development_Stage < 0.5:
            LAI_ChangeRate = (self.Specific_Leaf_N_Bottom * self.Nitrogen_Leaf_ChangeRate - self.Nitrogen_Leaf * Specific_Leaf_N_Bottom_ChangeRate) / self.Specific_Leaf_N_Bottom / (self.Specific_Leaf_N_Bottom + LeafNitrogen_ExtinctionCoefficient * self.Nitrogen_Leaf)
    
        self.LAI_ChangeRate = LAI_ChangeRate
        self.Specific_Leaf_N_Bottom_ChangeRate = Specific_Leaf_N_Bottom_ChangeRate


    
    
    def Calculate_Respiration_Rates(self, Nitrogen_uptake, Root_carbon_loss_rate_senescence):
        MaintenanceLeaf_DarkRespiration = 0.06 * (1 - self.Fraction_Shoot) * self.Photo_Assimilate
        MaintenanceStem_UndergroundStorage = 0.06 * 0.05 / 0.454 * self.Growth_Efficiency_Veg * self.Photo_Assimilate
        MaintenanceNitrogen_UptakeNitrate = 44 / 12 * 2.05 * Nitrogen_uptake
        # MaintenanceNitrogen_UptakeAmmonium = 44.0 / 12.0 * 0.17 * NitrogenUptake_Ammonium
        MaintenanceNitrogen_UptakeAmmonium=0
        Rate_Respiration_MineralUptakes = (MaintenanceNitrogen_UptakeNitrate + MaintenanceNitrogen_UptakeAmmonium + MaintenanceStem_UndergroundStorage + MaintenanceLeaf_DarkRespiration - self.Maintenance_Respiration) / self.Model_TimeStep
        GrowthRespiration = 44.0 / 12.0 * (((1.0 - self.Growth_Efficiency_Veg) / self.Growth_Efficiency_Veg * (self.Leaf_Carbon_ChangeRate + self.Stem_Carbon_ChangeRate + self.Root_Carbon_ChangeRate + self.Leaf_Carbon_Loss_ChangeRate + Root_carbon_loss_rate_senescence)) + ((1.0 - self.Growth_Efficiency_Veg) / self.Growth_Efficiency_Veg * self.Seed_Carbon_ChangeRate))
        TotalRespiration = self.Maintenance_Respiration + self.Respired_Carbon + GrowthRespiration + 44.0 / 12.0 * 0.06 * (self.Remobilized_Carbon_Stem_to_Seed + self.Remobilized_Carbon_Root_to_Seed)
    
        self.MaintenanceLeaf_DarkRespiration = MaintenanceLeaf_DarkRespiration
        self.MaintenanceStem_UndergroundStorage = MaintenanceStem_UndergroundStorage
        self.MaintenanceNitrogen_UptakeNitrate = MaintenanceNitrogen_UptakeNitrate
        self.MaintenanceNitrogen_UptakeAmmonium = MaintenanceNitrogen_UptakeAmmonium
        self.Rate_Respiration_MineralUptakes = Rate_Respiration_MineralUptakes
        self.GrowthRespiration = GrowthRespiration
        self.TotalRespiration = TotalRespiration

    def Calculate_Carbon_Nitrogen_Returns(self, Root_carbon_loss_rate_senescence, Root_nitrogen_loss_rate_senescence, average_soil_temperature):
        LeafDetritus_Carbon = (self.DeadLeaf_Carbon - self.DeadLeaf_Carbon_Litter) / 10. * (average_soil_temperature - self.BTemp_Phen) / (self.OTemp_Phen - self.BTemp_Phen)
        LitterCarbon_Total = Root_carbon_loss_rate_senescence + LeafDetritus_Carbon
        LitterNitrogen_Total = Root_nitrogen_loss_rate_senescence + LeafDetritus_Carbon / self.CarbonFrac_Veg * self.MinLeafN_Conc 
        NitrogenReturns_Total = self.Cumulative_LitterNitrogen + self.Switch_Function(self.Development_Stage - 2, 0, self.Nitrogen_Leaf + self.Nitrogen_Stem + self.Nitrogen_Root + self.NitrogenFixed_Reserve_pool + (self.DeadLeaf_Carbon - self.DeadLeaf_Carbon_Litter) / self.CarbonFrac_Veg * self.MinLeafN_Conc * (1) / 2)
    
        self.LeafDetritus_Carbon = LeafDetritus_Carbon
        self.LitterCarbon_Total = LitterCarbon_Total
        self.LitterNitrogen_Total = LitterNitrogen_Total
        self.NitrogenReturns_Total = NitrogenReturns_Total


    # def Check_Carbon_Nitrogen_Balance(self, TotalNitrogenUptake):
    #     CarbonCheck_Input = self.Carbon_Total + self.DeadLeaf_Carbon + self.DeadRoot_Dry_Weight - self.Initial_Leaf_Carbon - self.Initial_Root_Carbon
    #     CarbonCheck_Balance = (CarbonCheck_Input - (self.TotalPhotosynthesis_Canopy - self.Cumulative_Respiration) * 12.0 / 44.0) / self.Avoid_Zero_Division(CarbonCheck_Input) * 100
    #     NitrogenCheck_Input = self.Nitrogen_Total + self.Nitrogen_DeadLeaf + self.Nitrogen_DeadRoot - self.Initial_Leaf_N - self.Initial_Root_N
    #     NitrogenCheck_Balance = (NitrogenCheck_Input - TotalNitrogenUptake) / self.Avoid_Zero_Division(TotalNitrogenUptake) * 100
    
    #     self.CarbonCheck_Balance = CarbonCheck_Balance
    #     self.NitrogenCheck_Balance = NitrogenCheck_Balance




        
    def Switch_Function(self,x, y1, y2):

        if x < 0:
            out = y1
        else:
            out = y2

        return out
    
    def Avoid_Zero_Division(self, x):

        if x != 0:
            out = x
        else:
            out = 1.0
        return out
    
    def Indicator_Function(self, x1, x2):
        if x1 > 0 and x2 > 0:
            out = 1.0
        else:
            out = 0.0
        return out
    
    def Limit_Function(self,xl, xh, x):
        if xl <= x <= xh:
            out = x
        elif x > xh:
            out = xh
        else:
            out = xl
        return out
    
    def REANOR_Function(self,x1, x2):
        if x1 <= 0 and x2 <= 0:
            out = 1.0
        else:
            out = 0.0
        return out



    def Update_State_Variables(self, Root_nitrogen_loss_rate_senescence):
        self.Development_Stage += self.DevelopmentRate
    
        self.Cumulative_Thermal_Unit += self.Daily_Thermal_Unit
        
        self.Leaf_Carbon += self.Leaf_Carbon_ChangeRate
        self.DeadLeaf_Carbon += self.Leaf_Carbon_Loss_ChangeRate
        self.Stem_Carbon += self.Stem_Carbon_ChangeRate
        self.Seed_Carbon += self.Seed_Carbon_ChangeRate
        self.Root_Carbon += self.Root_Carbon_ChangeRate
        
        self.DeadLeaf_Carbon_Litter += self.LeafDetritus_Carbon #Rate of transfer of carbon from dead leaves to litter

        
        self.Nitrogen_Root += self.Nitrogen_Root_ChangeRate
        self.Nitrogen_Stem += self.Nitrogen_Stem_ChangeRate
        self.Nitrogen_Leaf += self.Nitrogen_Leaf_ChangeRate
        self.Nitrogen_Seed += self.Nitrogen_Seed_ChangeRate
        self.Total_Leaf_Nitrogen += self.Nitrogen_Leaf_ChangeRate_Positive
        self.Nitrogen_DeadLeaf += self.Leaf_Nitrogen_Loss_ChangeRate
    
        self.Stem_CarbonReserve += self.Stem_CarbonReserve_ChangeRate
        self.Root_CarbonReserve += self.Root_CarbonReserve_ChangeRate
        
        
        self.Remobilizable_N_seedgrowth_after_seedfilling += self.Remobilizable_N_seedgrowth_after_seedfilling_Rate
        self.Remobilizable_N_seedgrowth_before_seedfilling += self.Remobilizable_N_seedgrowth_before_seedfilling_Rate
        
        self.CarbonDemand_Deficit_ForSeedFill_PreviousTimeSteps += self.CarbonDemand_Deficit_ForSeedFill_ChangeRate
        self.CarbonDemand_Deficit_Stem_PreviousTimeSteps += self.CarbonDemand_Deficit_Stem_ChangeRate
        self.CarbonDemand_Stem_PreviousTimeSteps += self.CarbonDemand_Stem_ChangeRate

        self.Specific_Leaf_N_Bottom += self.Specific_Leaf_N_Bottom_ChangeRate
        self.Carbon_determined_LAI += self.LAI_ChangeRate
        
        
        self.Respiration_MineralUptakes += self.Rate_Respiration_MineralUptakes
        self.Canopy_Nitrogen_Demand_PreviousTimeStep += self.NitrogenDemand_ChangeRate
        self.Canopy_Nitrogen_Supply_PreviousTimeStep += self.NitrogenSupply_ChangeRate
        self.Total_NitrogenFixed += self.Nitrogen_Fixed
        self.NitrogenFixed_Reserve_Pool += self.NitrogenFixation_Reserve_Pool_ChangeRate
        
        self.Plant_Height += self.Plant_Height_ChangeRate
        
        self.Cumulative_Respiration += self.TotalRespiration
        self.Cumulative_LitterNitrogen += self.LitterNitrogen_Total
        self.Cumulative_Canopy_PhotosyntheticallyActiveRadiation += self.Canopy_PhotosyntheticallyActiveRadiation
        
        self.Cumulative_Actual_Canopy_Photosynthesis += self.Actual_Canopy_Photosynthesis
        self.Cumulative_Actual_Canopy_Transpiration += self.Actual_Canopy_Transpiration
        
        self.Nitrogen_DeadRoot += Root_nitrogen_loss_rate_senescence

        
 
 
 













