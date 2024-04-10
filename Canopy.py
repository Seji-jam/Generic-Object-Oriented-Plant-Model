import math
import Leaf
import numpy as np

class Canopy:
    def __init__(self, BTemp_Phen, OTemp_Phen, CTemp_Phen, TempCurve_Res, CropType_Photoperiod, StartPhotoperiod_Phase, EndPhotoperiod_Phase, Photoperiod_Sensitivity,
                 MinThermal_Day_Veg, MinThermal_Day_Rep,
                 CarbonAlloc_Shoot, NitrogenAlloc_Shoot, Plant_Density, Seed_Weight, 
                 Crop_TypeDet, EndSeedNum_DetPeriod, SeedN_RemobFract, SLA_Const, MinSLN_Photosyn, MinRootN_Conc,
                 CarbonCost_NFix, MaxN_Uptake, MinStemN_Conc,
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
        self.MinSLN_Photosyn = MinSLN_Photosyn
        self.MinRootN_Conc = MinRootN_Conc
        self.MaxSeedGrowth_DS = MaxSeedGrowth_DS
        self.MaxStemGrowth_DS = MaxStemGrowth_DS
        self.StemDW_Height = StemDW_Height
        self.CarbonCost_NFix = CarbonCost_NFix
        self.MaxN_Uptake = MaxN_Uptake
        self.legume = legume
        self.MinStemN_Conc = MinStemN_Conc
        self.IniLeafN_Conc = IniLeafN_Conc
        
        #initializing intermediate variables
        self.GrowthEfficiency = 0.25  # Growth efficiency
        self.CarbonFrac_Veg = 0.48  # Carbon fraction in vegetative biomass
        self.FractionProtein_StorageOrgans = 0.13  # Fraction of protein in storage organs
        self.FractionCarbs_StorageOrgans = 0.71  # Fraction of carbohydrates in storage organs
        self.CarbonFraction_Organs = 0.47  # Carbon fraction in the storage organs
        self.YieldGrowthOutput = 0.8  # Growth efficiency for storage organs
        
        self.InitCarbon_LivingVegBiomass = self.Plant_Density * self.Seed_Weight * self.CarbonFraction_Organs * self.GrowthEfficiency * self.CarbonAlloc_Shoot
        self.InitCarbon_ReproductiveTissues = self.Plant_Density * self.Seed_Weight * self.CarbonFraction_Organs * self.GrowthEfficiency * (1.0 - self.CarbonAlloc_Shoot)
        self.InitLeafAreaIndex = self.InitCarbon_LivingVegBiomass / self.CarbonFrac_Veg * self.SLA_Const
        self.InitNitrogen_LivingVegBiomass = self.IniLeafN_Conc * self.InitCarbon_LivingVegBiomass / self.CarbonFrac_Veg
        self.InitPlant_Height = self.MaxPlant_Height / 1000.0  # Converting to a different unit if necessary
        self.MinLeafN_Conc = self.SLA_Const * self.MinSLN_Photosyn
        self.InitNitrogen_ReproductiveTissues = (self.Plant_Density * self.Seed_Weight * self.GrowthEfficiency * self.IniLeafN_Conc * self.CarbonAlloc_Shoot / self.NitrogenAlloc_Shoot) - self.InitNitrogen_LivingVegBiomass
        self.InitSpecificLeafN_Biomass = self.InitNitrogen_LivingVegBiomass / self.InitLeafAreaIndex


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
    
        # Totals for actual conditions
        self.Total_ActualCanopy_Photosynthesis = 0
        self.Total_ActualCanopy_Transpiration = 0
        self.Total_CanopyPhotosyntheticallyActiveRadiation = 0

        self.DevelopmentRate = 0
        self.NitrogenDemand_Allocation = 0
        self.Nitrogen_Demand = 0
        self.Shoot_Harvestable_Sugar_Allocation = 0
        self.End_SeedFill_Dynamic = 0
        self.Fraction_Carbon_Seed = 0
        self.Fraction_Carbon_Stem = 0
        self.Fraction_Carbon_Leaf = 0
        self.Fraction_Carbon_StemReserve = 0
        self.Fraction_Carbon_RootReserve = 0
        self.Carbon_Supply_RootTotal = 0
        self.Carbon_Supply_StemTotal = 0
        self.Carbon_Supply_Seed = 0
        self.Carbon_Supply_Leaf = 0
        self.LAI_Current = self.Initial_LAI
        self.LAI_Next = self.Initial_LAI
    
        self.NitrogenFixation_Rate = 0
        self.Rate_SpecificLeafNitrogen_Bottom = 0
    
        self.Nitrogen_DeadRoot = 0
    
        self.Respiration_Main = 0
        self.RootWeight_Dead = 0
    
        self.Respired_Carbon = 0
        self.Carbon_Shoot = 0
        self.Carbon_Root = 0
        self.Nitrogen_Shoot = 0

        self.delt= Model_TimeStep

        self.Carbon_DeadRoot = 0

        self.DailyCarbonDemand_SeedFill_Cumulative = 0
        self.DailyCarbonDemand_Seed = 0
        self.DailyCarbon_Supply_Stem = 0
    
        self.FlowLimit_Carbon_Seed = 0
        self.FlowLimit_Carbon_Stem = 0
        self.DailyCarbonDemand_StemGrowth = 0
        self.DailyCarbon_Supply_Total = 0
        self.Rate_Height = 0
        self.RateDaily_CarbonDemand_TotalPrev = 0
        self.Fraction_Carbon_Shoot = 0
    
        self.ReserveCarbon_Shoot_Update = 0
        self.ReserveCarbon_Root_Update = 0
        self.CarbonReserve_Shoot_Increase = 0
        self.CarbonReserve_Root_Increase = 0
    
        self.TotalPhotosynthesis_Canopy = 0
        self.Rate_LeafAreaIndex = 0
    
        # New additions
        self.Carbon_DeadRoot_Total = 0
    
        self.DailyCarbonDemand_SeedFill_Current = 0
        self.DailyCarbonDemand_Seed_Current = 0
        self.DailyCarbon_Supply_Stem_Current = 0
    
        self.FlowLimit_CarbonSeed_Current = 0
        self.FlowLimit_CarbonStem_Current = 0
        self.DailyCarbonDemand_StemGrowth_Current = 0
        self.DailyCarbon_Supply_Total_Current = 0
        self.Rate_Height_Current = 0
        self.RateDaily_CarbonDemand_TotalPrev_Current = 0
    
        self.ReserveCarbon_Shoot_Current = 0
        self.ReserveCarbon_Root_Current = 0
            
        
        
        # State Variables
        # Development stage and thermal units
        self.DevelopmentStage = 0
        self.CumulativeThermalUnits = 0
    
        # Carbon in various plant components
        self.Carbon_Leaf = self.Carbon_InitialLeaf
        self.Carbon_DeadLeaf = 0
        self.Carbon_Stem = 0
        self.Carbon_Seed = 0
        self.Carbon_Root = self.Carbon_InitialRoot
        self.Carbon_DeadLeafShed = 0
    
        # Nitrogen in various plant components
        self.Nitrogen_Root = self.Nitrogen_InitialRoot
        self.Nitrogen_Stem = 0
        self.Nitrogen_Leaf = self.Nitrogen_InitialLeaf
        self.TotalNitrogen_Leaf = self.Nitrogen_InitialLeaf
        self.Nitrogen_Seed = 0
        self.Nitrogen_DeadLeaf = 0
    
        # Carbon reserves
        self.CarbonReserve_Shoot = 0
        self.CarbonReserve_Root = 0
    
        # Nitrogen dynamics
        self.NitrogenRemoval_EndSeedFill = 0
        self.NitrogenRemoval_EarlyFlowering = 0
    
        # Carbon dynamics
        self.DailyCarbonDemand_SeedReserve = 0
        self.DailyCarbonDemand_StemReserve = 0
    
        # Additional state variables
        self.SpecificLeafNitrogen_Bottom = 0
        self.Rate_LeafAreaIndex = 0
        self.RespirationMultiplier = 0
        self.NitrogenDemand_Previous = 0
        self.Nitrogen_Supply = 0
        self.NitrogenFixed_Total = 0
        self.NitrogenFixed_Rate = 0
        self.DailyCarbonDemand_TotalPrev = 0
    
        # Resetting rates of change and flows
        self.Respiration_Main = 0
        self.Respired_Carbon = 0
        self.FlowLimit_Carbon_Seed = 0
        self.FlowLimit_Carbon_Stem = 0

            # Specific Leaf Nitrogen content and Leaf Area Index
        self.SpecificLeafNitrogen_Bottom = self.SpecificLeafNitrogen_BottomInitial
        self.LAI_Current = self.Initial_LAI  # Initial Leaf Area Index
    
        # Respiration Multiplier and Nitrogen dynamics
        self.RespirationMultiplier = 0
        self.NitrogenDemand_Previous = 0
        self.Nitrogen_Supply = 0
        self.NitrogenFixed_Total = 0
        self.NitrogenFixed_Rate = 0
    
        # Carbon dynamics for stem growth and seed filling
        self.DailyCarbonDemand_TotalPrev = 0
        self.Plant_Height = self.Initial_Plant_Height
        self.TotalRespiration = 0
        self.LitterNitrogen_TotalPrevious = 0
    
        # Photosynthesis, Transpiration, and Canopy PAR
        self.Total_CanopyPhotosyntheticallyActiveRadiation = 0
        self.Total_ActualCanopy_Photosynthesis = 0
        self.Total_ActualCanopy_Transpiration = 0
            
         
      
        
        
        
        

    
    def Leaf_to_Canopy_Integration(self,hourly_data_list1,hourly_data_list2):
        wgauss = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
        hourly_data_Canopy = np.array([sum(x) for x in zip(hourly_data_list1, hourly_data_list2)])
        weighted_data = hourly_data_Canopy * wgauss
        return weighted_data.sum()


    def Update_Canopy_PAR(self,hourly_apar_SH, hourly_apar_SU,dayl):
        daily_average_canopy_PAR =  self.Leaf_to_Canopy_Integration(hourly_apar_SH, hourly_apar_SU)
        self.Canopy_PAR=daily_average_canopy_PAR * dayl*  3600




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
        daily_average_canopy_photosyn =  self.Leaf_to_Canopy_Integration(photosyn_SU, photosyn_SH)
        State= State.lower()
        if State=='p':
            self.Potential_canopy_photosyn=daily_average_canopy_photosyn * dayl*  3600
        elif State=='a':
            self.Actual_canopy_photosyn=daily_average_canopy_photosyn * dayl*  3600
        else:
            raise ValueError("Invalid State provided. State must be 'A' for Actual or 'P' for Potential.")


    def Update_Canopy_Transpiration(self,State,transpiration_SU, transpiration_SH,dayl):
        daily_average_canopy_transpiration =  self.Leaf_to_Canopy_Integration(transpiration_SU, transpiration_SH)
        State= State.lower()

        if State=='p':
            self.Potential_canopy_transpiration=daily_average_canopy_transpiration * dayl*  3600
        elif State=='a':
            self.Actual_canopy_transpiration=daily_average_canopy_transpiration * dayl*  3600
        else:
            raise ValueError("Invalid State provided. State must be 'A' for Actual or 'P' for Potential.")

        
   

    def calculate_thermal_units(self, ds, Tmax, Tmin, DayLength, BTemp_Phen, OTemp_Phen, CTemp_Phen, TempCurve_Res):
        
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
            if ds > 1:
                TempHourly = min(TempHourly, OTemp_Phen)
    
            # Instantaneous thermal unit based on bell-shaped temperature response
            if TempHourly < BTemp_Phen or TempHourly > CTemp_Phen:
                ThermalUnit = 0
            else:
                ThermalUnit = (((CTemp_Phen - TempHourly) / (CTemp_Phen - OTemp_Phen)) * ((TempHourly - BTemp_Phen) / (OTemp_Phen - BTemp_Phen)) ** ((OTemp_Phen - BTemp_Phen) / (CTemp_Phen - OTemp_Phen))) ** TempCurve_Res
            TotalThermalUnits += ThermalUnit / 24
            
        # Daily thermal unit
        self.DailyThermalUnit = TotalThermalUnits
    
        
    

        
    def calculate_development_rate(self, DevelopmentStage, CropType_Photoperiod, DayLength_PhotoPeriod, StartPhotoperiod_Phase, EndPhotoperiod_Phase, Photoperiod_Sensitivity, MinThermal_Day_Veg, MinThermal_Day_Rep, DailyThermalUnit):
        # Determining if it is for short-day or long-day crop
        if CropType_Photoperiod < 0:
            OptimumPhotoperiod = 18  # Minimum optimum photoperiod for long-day crop
            DayLengthResponse = min(OptimumPhotoperiod, DayLength_PhotoPeriod)
        else:
            OptimumPhotoperiod = 11  # Maximum optimum photoperiod for short-day crop
            DayLengthResponse = max(OptimumPhotoperiod, DayLength_PhotoPeriod)
    
        # Effect of photoperiod on development rate
        if DevelopmentStage < StartPhotoperiod_Phase or DevelopmentStage > EndPhotoperiod_Phase:
            EffectPhotoperiod = 1
        else:
            EffectPhotoperiod = max(0, 1 - Photoperiod_Sensitivity * (DayLengthResponse - OptimumPhotoperiod))
    
        # Development rate of vegetative and reproductive phases
        if 0 <= DevelopmentStage < 1.0:
            DevelopmentRate = 1 / MinThermal_Day_Veg * DailyThermalUnit * EffectPhotoperiod
        else:
            DevelopmentRate = 1 / MinThermal_Day_Rep * DailyThermalUnit
    
        self.DevelopmentRate = DevelopmentRate





    def Initialize_Biomass_Formation(self):
        Seed_Dry_Weight = self.Seed_Carbon / self.CarbonFraction_Organs # Dry weight of seed
        LiveLeaf_Dry_Weight = self.LiveLeaf_Carbon / self.CarbonFrac_Veg # Dry weight of live leaves
        Stem_Dry_Weight = self.Stem_Carbon / self.CarbonFrac_Veg + self.ReserveStem_Carbon / 0.444 # Dry weight of stems
        LiveRoot_Dry_Weight = self.LiveRoot_Carbon / self.CarbonFrac_Veg + self.ReserveRoot_Carbon / 0.444 # Dry weight of live roots
    
        Shoot_Dry_Weight = LiveLeaf_Dry_Weight + Stem_Dry_Weight + Seed_Dry_Weight # Dry weight of live shoot (above-ground) organs
    
        TotalLiveOrgan_Dry_Weight = Shoot_Dry_Weight + LiveRoot_Dry_Weight # Dry weight of total live organs
    
        DeadLeaf_Dry_Weight = self.DeadLeaf_Carbon / self.CarbonFrac_Veg # Dry weight of dead leaves
    
        Shoot_Dry_Weight_ExcShedLeaves = Shoot_Dry_Weight + (DeadLeaf_Dry_Weight - self.ShedLeaf_Carbon / self.CarbonFrac_Veg) # Dry weight of shoot organs (excluding shedded leaves)
    
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
        Nitrogen_Shoot_ExcShedLeaves = Nitrogen_Shoot + (LiveLeaf_Dry_Weight - self.ShedLeaf_Carbon / self.CarbonFrac_Veg) * self.MinLeafN_Conc
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
     Carbon_Shoot = self.LiveLeaf_Carbon + self.Stem_Carbon + self.ReserveStem_Carbon + self.Seed_Carbon
     Carbon_Root = self.LiveRoot_Carbon + self.ReserveRoot_Carbon
     Carbon_Total = Carbon_Shoot + Carbon_Root
     self.Carbon_Total = Carbon_Total
     self.Carbon_Shoot = Carbon_Shoot
     self.Carbon_Root = Carbon_Root



    def Calculate_Respiration(self):
        Maintenance_Respiration_Extra = max(min(44./12.*0.218*(self.Nitrogen_Total-self.Shoot_Dry_Weight*self.MinLeafN_Conc-self.LiveRoot_Dry_Weight*self.MinRootN_Conc), self.Actual_Canopy_Photosynthesis-1.E-5-self.Respiration_Main), 0.)
        Maintenance_Respiration = max(0., min(self.Actual_Canopy_Photosynthesis-1.E-5, self.Respiration_Multiplier) + Maintenance_Respiration_Extra)
        self.Respiration_Main = Maintenance_Respiration

    def Calculate_Nitrogen_Fixation(self):
        Nitrogen_Fixed_Demand = max(0., self.Nitrogen_Demand - self.Nitrogen_Supply)
        Nitrogen_Fixed_Efficiency = max(0., self.Actual_Canopy_Photosynthesis - 1.E-5 - self.Respiration_Main) / self.CarbonCost_NFix * 12. / 44.
        Nitrogen_Fixed = self.Switch_Function(self.Is_Legume, 0., min(Nitrogen_Fixed_Efficiency, Nitrogen_Fixed_Demand))
        self.Nitrogen_Fixed = Nitrogen_Fixed

    def Calculate_Photo_Assimilates(self):
        Respired_Carbon = 44. / 12. * (self.CarbonCost_NFix * self.Nitrogen_Fixed)
        Photo_Assimilate = self.Actual_Canopy_Photosynthesis - self.Respiration_Main - Respired_Carbon
        self.Photo_Assimilate = Photo_Assimilate
        self.Respired_Carbon = Respired_Carbon

    


    def Calculate_Crop_Nitrogen_Demand(self, Specific_Leaf_Nitrogen):
        Shoot_Harvestable_Sugar_Allocation = 12.0 / 44.0 * self.YieldGrowth_Veg * (self.Actual_Canopy_Photosynthesis - self.Respiration_Main - self.Respired_Carbon) / self.Carbon_Shoot
        Nitrogen_Demand_Allocation = self.Carbon_Root * Shoot_Harvestable_Sugar_Allocation**2 
        High_Nitrogen_Conc_Crit = self.IniLeafN_Conc * np.exp(-0.4 * self.DevelopmentStage)
        Nitrogen_Demand_Dynamic = self.Switch_Function(self.DevelopmentStage - 1.0, self.Shoot_Dry_Weight * (High_Nitrogen_Conc_Crit - self.Shoot_Nitrogen_Conc) * (1.0 + self.Nitrogen_Root / self.Nitrogen_Shoot) / self.Model_TimeStep, 0.0)
        Nitrogen_Demand_Adjusted = self.Switch_Function(self.Leaf_Nitrogen_Conc - 1.5 * self.IniLeafN_Conc, max(Nitrogen_Demand_Allocation, Nitrogen_Demand_Dynamic), 0)
        Nitrogen_Demand = self.Switch_Function(self.MinSLN_Photosyn - Specific_Leaf_Nitrogen + 1.0e-5, min(self.MaxN_Uptake, Nitrogen_Demand_Adjusted), 0)  # Crop nitrogen demand
    
        print(self.Carbon_Root, Shoot_Harvestable_Sugar_Allocation)
        self.Nitrogen_Demand_Allocation = Nitrogen_Demand_Allocation
        self.Nitrogen_Demand = Nitrogen_Demand
    



    def Calculate_Nitrogen_Partitioning(self, Specific_Leaf_Nitrogen_Target):
        Nitrogen_Carbon_Ratio = self.Switch_Function(Specific_Leaf_Nitrogen_Target - self.MinSLN_Photosyn, 0, min(self.MaxN_Uptake, self.Nitrogen_Demand_Allocation)) / (self.YieldGrowth_Veg * (self.Actual_Canopy_Photosynthesis - self.Respiration_Main - self.Respired_Carbon) * 12 / 44)
        Shoot_Harvestable_Sugar_Allocation = 12.0 / 44.0 * self.YieldGrowth_Veg * (self.Actual_Canopy_Photosynthesis - self.Respiration_Main - self.Respired_Carbon) / self.Carbon_Shoot
        Fraction_Shoot_Nitrogen = 1 / (1 + Nitrogen_Carbon_Ratio * 1 / Shoot_Harvestable_Sugar_Allocation * self.Carbon_Shoot / self.Carbon_Root * self.Nitrogen_Root / self.Nitrogen_Shoot)
        
        print(Fraction_Shoot_Nitrogen)
        self.Shoot_Harvestable_Sugar_Allocation = Shoot_Harvestable_Sugar_Allocation
        self.Fraction_Shoot_Nitrogen = Fraction_Shoot_Nitrogen
        self.Nitrogen_Carbon_Ratio = Nitrogen_Carbon_Ratio



    def Calculate_Seed_Properties(self, Determinate, End_SeedFill_Indeterminate, Nitrogen_Removal_EarlyFlowering, Nitrogen_Removal_EndSeedFill, Fraction_SeedNitrogen_Remobilizable, Expected_SeedNitrogen_Concentration, Seed_Weight):
        End_SeedFill_Dynamic = self.Switch_Function(Determinate, End_SeedFill_Indeterminate, 1.)
        Nitrogen_Remobilization = Nitrogen_Removal_EarlyFlowering + (Nitrogen_Removal_EndSeedFill - Nitrogen_Removal_EarlyFlowering) * (End_SeedFill_Dynamic - 1.0) / self.Avoid_Zero_Division(min(self.DevelopmentStage, End_SeedFill_Dynamic) - 1)
        Total_Seed_Nitrogen = Nitrogen_Remobilization / Fraction_SeedNitrogen_Remobilizable / Expected_SeedNitrogen_Concentration / Seed_Weight
        Total_Seed_Weight = self.Seed_Carbon / self.CarbonFraction_Organs / self.Avoid_Zero_Division(Total_Seed_Nitrogen) * 1000
        
        self.End_SeedFill_Dynamic = End_SeedFill_Dynamic
        self.Total_Seed_Nitrogen = Total_Seed_Nitrogen
        self.Total_Seed_Weight = Total_Seed_Weight

    
    def Calculate_Senescence(self, LAI_Current, LAI_Next, Time_Step):
        LeafWeight_Loss_Max = (LAI_Current - min(LAI_Current, LAI_Next)) / self.SLA_Const / Time_Step
        LeafWeight_Loss = min(self.LiveLeaf_Dry_Weight - 1.e-5, LeafWeight_Loss_Max + self.FUN_REANOR(self.End_SeedFill_Dynamic - self.DevelopmentStage, LeafWeight_Loss_Max) * 0.03 * self.LiveLeaf_Dry_Weight)
        Nitrogen_Leaf_Loss = min(LeafWeight_Loss, LeafWeight_Loss_Max) * self.MinLeafN_Conc + (LeafWeight_Loss - min(LeafWeight_Loss, LeafWeight_Loss_Max)) * (self.Nitrogen_Leaf / self.LiveLeaf_Dry_Weight)
        Carbon_Leaf_Loss = LeafWeight_Loss_Max * self.CarbonFrac_Veg
        
        print(Nitrogen_Leaf_Loss)
        self.Nitrogen_Leaf_Loss = Nitrogen_Leaf_Loss
        self.Carbon_Leaf_Loss = Carbon_Leaf_Loss




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
        SeedFill_Day = self.Limit_Function(1, 2, self.DevelopmentStage) - 1 # Day of seed filling
    
        DailySeed_Fill_Rate = self.DevelopmentRate * (2 * SeedFill_Phase - SeedFill_Transition) * (SeedFill_Phase - SeedFill_Day) / SeedFill_Phase / (SeedFill_Phase - SeedFill_Transition) ** 2 * (SeedFill_Day / SeedFill_Phase) ** (SeedFill_Transition / (SeedFill_Phase - SeedFill_Transition))
    
        TotalCarbon_Seeds = self.Total_Seed_Nitrogen * self.Seed_Weight * self.CarbonFraction_Organs
        DailyCarbonDemand_SeedFill = self.Switch_Function(self.DevelopmentStage-SeedFill_Start, 0., TotalCarbon_Seeds/self.YieldGrowth_Output*DailySeed_Fill_Rate) 
        DailyCarbonSupply_Seed = DailyCarbonDemand_SeedFill + max(0, self.DailyCarbonDemand_SeedReserve) / self.Model_TimeStep
        FlowLimit_Carbon_Seed = min(DailyCarbonSupply_Seed, self.DailyCarbon_Supply_Shoot)
        
        self.DailyCarbonDemand_SeedFill = DailyCarbonDemand_SeedFill
        self.DailyCarbonSupply_Seed = DailyCarbonSupply_Seed
        self.FlowLimit_Carbon_Seed = FlowLimit_Carbon_Seed



    def Calculate_Carbon_Flow_Stem_Growth(self):
        DailyCarbon_Supply_Stem = self.DailyCarbon_Supply_Shoot - self.FlowLimit_Carbon_Seed  # Daily carbon supply from current photosynthesis for structural stem growth
    
        IntegralFactor_Stress_Height = self.Limit_Function(0, 1, DailyCarbon_Supply_Stem / self.Avoid_Zero_Division(self.DailyCarbonDemand_TotalPrev))  # Integral factor of stresses on plant height growth
        Transition_Height = self.MaxStemGrowth_DS * (1. + self.End_SeedFill_Dynamic) / 2. # Transition period for height growth
        Phase_Height = (1. + self.End_SeedFill_Dynamic) / 2.  # Phase length for height growth
        SeedFill_Day = min((1. + self.End_SeedFill_Dynamic) / 2., self.DevelopmentStage)  # Day of height growth phase
    
        DailyHeight_Growth_Rate = self.DevelopmentRate * (2 * Phase_Height - Transition_Height) * (Phase_Height - SeedFill_Day) / Phase_Height / (Phase_Height - Transition_Height) ** 2 * (SeedFill_Day / Phase_Height) ** (Transition_Height / (Phase_Height - Transition_Height))
    
        Start_Growth = 0 # Start of structural stem growth
        TotalCarbon_Stem = self.StemDW_Height * self.MaxPlant_Height * self.CarbonFrac_Veg
        DailyCarbonDemand_StemGrowth = self.Switch_Function(self.DevelopmentStage - Start_Growth, 0., TotalCarbon_Stem / self.YieldGrowth_Veg * DailyHeight_Growth_Rate * IntegralFactor_Stress_Height)
        DailyCarbon_Supply_Total = DailyCarbonDemand_StemGrowth + max(0, self.DailyCarbonDemand_StemReserve) / self.Model_TimeStep
        FlowLimit_Carbon_Stem = min(DailyCarbon_Supply_Total, DailyCarbon_Supply_Stem)
    
        # Daily carbon flow for structural stem growth
        Rate_Height = min(self.MaxPlant_Height - self.Plant_Height, DailyHeight_Growth_Rate * self.MaxPlant_Height * IntegralFactor_Stress_Height)  # Rate of plant height growth
        RateDaily_CarbonDemand_TotalPrev = (DailyCarbonDemand_StemGrowth - self.DailyCarbonDemand_TotalPrev) / self.Model_TimeStep  # Carbon demand for structural stem growth at the previous time step
        
        self.DailyCarbon_Supply_Stem = DailyCarbon_Supply_Stem
        self.FlowLimit_Carbon_Stem = FlowLimit_Carbon_Stem
        self.DailyCarbonDemand_StemGrowth = DailyCarbonDemand_StemGrowth
        self.DailyCarbon_Supply_Total = DailyCarbon_Supply_Total
        self.Rate_Height = Rate_Height
        self.RateDaily_CarbonDemand_TotalPrev = RateDaily_CarbonDemand_TotalPrev
    




    
    def Calculate_Carbon_Partitioning(self, LAI_Next, LAI_Current, CarbonSupply_RootTotalNew):
        Fraction_Carbon_Seed = self.FlowLimit_Carbon_Seed / self.DailyCarbon_Supply_Shoot
        Fraction_Carbon_Stem = self.Switch_Function(self.DevelopmentStage - (self.End_SeedFill_Dynamic + 0.2), self.FlowLimit_Carbon_Stem / self.DailyCarbon_Supply_Shoot, 0)
        Fraction_Carbon_Leaf = self.Indicator_Function(LAI_Next - LAI_Current, self.End_SeedFill_Dynamic - self.DevelopmentStage) * (1.0 - Fraction_Carbon_Seed - Fraction_Carbon_Stem)
        Fraction_Carbon_StemReserve = 1.0 - Fraction_Carbon_Leaf - Fraction_Carbon_Seed - Fraction_Carbon_Stem  # Fraction of new shoot carbon to stem reserves
        Fraction_Carbon_RootReserve = self.Switch_Function(CarbonSupply_RootTotalNew - self.LiveRoot_Carbon, 1.0, 0.0)  # Fraction of new root carbon to root reserves
    
        self.Fraction_Carbon_Seed = Fraction_Carbon_Seed
        self.Fraction_Carbon_Stem = Fraction_Carbon_Stem
        self.Fraction_Carbon_Leaf = Fraction_Carbon_Leaf
        self.Fraction_Carbon_StemReserve = Fraction_Carbon_StemReserve
        self.Fraction_Carbon_RootReserve = Fraction_Carbon_RootReserve

    def Update_Carbon_Reserve_Pools(self):
        Gap_Carbon = max(0, self.DailyCarbonDemand_Seed - self.DailyCarbon_Supply_Shoot)
        CarbonReserve_ShootIncrease = min(0.94 * self.CarbonReserve_Shoot, self.CarbonReserve_Shoot / self.Avoid_Zero_Division(self.CarbonReserve_Shoot + self.CarbonReserve_Root) * Gap_Carbon) / 0.94
        CarbonReserve_RootIncrease = min(0.94 * self.CarbonReserve_Root, self.CarbonReserve_Root / self.Avoid_Zero_Division(self.CarbonReserve_Shoot + self.CarbonReserve_Root) * Gap_Carbon) / 0.94
        CarbonReserve_Shoot = self.Switch_Function(self.DailyCarbonDemand_Seed - self.DailyCarbon_Supply_Shoot, 0, CarbonReserve_ShootIncrease)
        CarbonReserve_Root = self.Switch_Function(self.DailyCarbonDemand_Seed - self.DailyCarbon_Supply_Shoot, 0, CarbonReserve_RootIncrease)
        Reserve_Carbon_Shoot = self.Fraction_Carbon_StemReserve * self.DailyCarbon_Supply_Shoot - CarbonReserve_Shoot
        Reserve_Carbon_Root = self.Fraction_Carbon_RootReserve * self.DailyCarbon_Supply_Root - CarbonReserve_Root
    
        self.Reserve_Carbon_Shoot = Reserve_Carbon_Shoot
        self.Reserve_Carbon_Root = Reserve_Carbon_Root
        self.CarbonReserve_Shoot = CarbonReserve_Shoot
        self.CarbonReserve_Root = CarbonReserve_Root

        
    def Calculate_Carbon_Production_Rate(self, LastCarbonRoot):
        Carbon_Supply_RootTotal = 12 / 44 * self.Photo_Assimilate * (1 - self.Fraction_Shoot) * (1 - self.Fraction_Carbon_RootReserve) * self.YieldGrowth_Veg - LastCarbonRoot
        Carbon_Supply_StemTotal = 12.0 / 44.0 * self.Photo_Assimilate * self.Fraction_Shoot * self.Fraction_Carbon_Stem * self.YieldGrowth_Veg
        Carbon_Supply_Seed = 12.0 / 44.0 * self.Photo_Assimilate * self.Fraction_Shoot * self.Fraction_Carbon_Seed * self.YieldGrowth_Output + 0.94 * (self.CarbonReserve_Shoot + self.CarbonReserve_Root) * self.YieldGrowth_Output
        Carbon_Supply_Leaf = 12.0 / 44.0 * self.Photo_Assimilate * self.Fraction_Shoot * self.Fraction_Carbon_Leaf * self.YieldGrowth_Veg - self.Carbon_Leaf_Loss
    
        self.Carbon_Supply_RootTotal = Carbon_Supply_RootTotal
        self.Carbon_Supply_StemTotal = Carbon_Supply_StemTotal
        self.Carbon_Supply_Seed = Carbon_Supply_Seed
        self.Carbon_Supply_Leaf = Carbon_Supply_Leaf


    def Calculate_Biomass_Change_Rate(self):
        RootWeight_Rate = self.CarbonSupply_RootTotal / self.CarbonFrac_Veg + self.ReserveCarbon_Root / 0.444
        SeedWeight_Rate = self.CarbonSupply_Seed / self.CarbonFraction_Organs
        LeafWeight_Rate = self.CarbonSupply_Leaf / self.CarbonFrac_Veg
        StemWeight_Rate = self.CarbonSupply_StemTotal / self.CarbonFrac_Veg + self.ReserveCarbon_Shoot / 0.444
    
        self.RootWeight_Rate = RootWeight_Rate
        self.SeedWeight_Rate = SeedWeight_Rate
        self.LeafWeight_Rate = LeafWeight_Rate
        self.StemWeight_Rate = StemWeight_Rate

    
    def Calculate_Carbon_Rates(self):
        ReserveDemand_SeedRate = max(0.0, (self.DailyCarbonDemand_SeedFill - self.CarbonSupply_Seed / self.YieldGrowth_Output)) - (self.FlowLimit_Carbon_Seed - min(self.DailyCarbonDemand_SeedFill, self.DailyCarbon_Supply_Shoot))
        ReserveDemand_StemRate = max(0., (self.DailyCarbonDemand_StemGrowth - self.CarbonSupply_StemTotal / self.YieldGrowth_Veg)) - (self.FlowLimit_Carbon_Stem - min(self.DailyCarbonDemand_StemGrowth, self.DailyCarbon_Supply_Stem))
    
        self.ReserveDemand_SeedRate = ReserveDemand_SeedRate
        self.ReserveDemand_StemRate = ReserveDemand_StemRate
    


    def Calculate_Nitrogen_Dynamics(self, TotalCarbonPhotosynth, NitrogenUptake):
        NitrogenFixation_Rate = self.Nitrogen_Fixed - min(self.Nitrogen_Demand, self.NitrogenFixed_Rate / TotalCarbonPhotosynth)
        NitrogenDemand_Rate = (self.Nitrogen_Demand - self.NitrogenDemand_Previous) / self.Model_TimeStep
        NitrogenSupply_Rate = (NitrogenUptake - self.Nitrogen_Supply) / self.Model_TimeStep
    
        self.NitrogenFixation_Rate = NitrogenFixation_Rate
        self.NitrogenDemand_Rate = NitrogenDemand_Rate
        self.NitrogenSupply_Rate = NitrogenSupply_Rate



    def Calculate_Seed_Number_Rate(self, NitrogenUptake, LastCarbonRoot):
        NitrogenReserve_Rate = NitrogenUptake - (self.MinLeafN_Conc * (self.CarbonSupply_Leaf + self.Carbon_Leaf_Loss) + self.MinRootN_Conc * (self.CarbonSupply_RootTotal + LastCarbonRoot) + self.StemN_Conc * self.CarbonSupply_StemTotal) / self.CarbonFrac_Veg
        NitrogenRemoval_EndSeedFill_Rate = self.Switch_Function(self.DevelopmentStage - self.End_SeedFill_Dynamic, NitrogenReserve_Rate, 0.0)
        NitrogenRemoval_EarlyFlowering_Rate = self.Switch_Function(self.DevelopmentStage - 1.0, NitrogenReserve_Rate, 0.0)
    
        self.NitrogenRemoval_EndSeedFill_Rate = NitrogenRemoval_EndSeedFill_Rate
        self.NitrogenRemoval_EarlyFlowering_Rate = NitrogenRemoval_EarlyFlowering_Rate

    
    def Calculate_Nitrogen_Accumulation_Rate(self, NitrogenUptake, LeafNitrogenLoss_Rate):
        # Amount of N partitioned to shoot
        Nitrogen_Shoot_New = self.Fraction_NewNitrogen_Shoot * NitrogenUptake
    
        # Leaf N or root N available for remobilization
        Nitrogen_Leaf_Available = self.Switch_Function(self.MinLeafN_Conc - self.Leaf_Nitrogen_Conc, self.Nitrogen_Leaf - self.LiveLeaf_Dry_Weight * self.MinLeafN_Conc, 0.) / self.Model_TimeStep
        Nitrogen_Root_Available = self.Switch_Function(self.MinRootN_Conc - self.Root_Nitrogen_Conc, self.Nitrogen_Root - self.LiveRoot_Dry_Weight * self.MinRootN_Conc, 0.) / self.Model_TimeStep
        Nitrogen_Total_Available = Nitrogen_Leaf_Available + Nitrogen_Root_Available
    
        # Rate of N accumulation in stem
        Nitrogen_Stem_Rate = self.StemWeight_Rate * self.Switch_Function(-Nitrogen_Total_Available, self.StemN_Conc, 0.)
    
        # Expected N dynamics during seed filling
        Expected_SeedNitrogen_Conc =  self.Standard_SeedNitrogen_Conc
    
        # Rate of N accumulation in seed
        Nitrogen_Seed_Growth = Nitrogen_Shoot_New - Nitrogen_Stem_Rate - Expected_SeedNitrogen_Conc * self.SeedWeight_Rate
        Nitrogen_Organic_Compound_Conc = max(0., self.Switch_Function(Nitrogen_Total_Available + Nitrogen_Seed_Growth, (Nitrogen_Total_Available + Nitrogen_Shoot_New - Nitrogen_Stem_Rate) / self.Avoid_Zero_Division(self.SeedWeight_Rate), Expected_SeedNitrogen_Conc))
        Nitrogen_Seed_Rate = self.SeedWeight_Rate * Nitrogen_Organic_Compound_Conc
    
        # Rate of N accumulation in leaf
        Nitrogen_Leaf_Accumulation = self.Switch_Function(Nitrogen_Total_Available + Nitrogen_Seed_Growth, -Nitrogen_Leaf_Available - self.LeafNitrogenLoss, -Nitrogen_Leaf_Available / self.Avoid_Zero_Division(Nitrogen_Total_Available) * (-Nitrogen_Seed_Growth) - self.LeafNitrogenLoss)
        LeafNitrogen_Growth = self.Switch_Function(Nitrogen_Seed_Growth, Nitrogen_Leaf_Accumulation, Nitrogen_Shoot_New - Nitrogen_Stem_Rate - Nitrogen_Seed_Rate - self.LeafNitrogenLoss)
        Nitrogen_Leaf_Rate = max(-self.Nitrogen_Leaf + 1E-7, LeafNitrogen_Growth)
        Nitrogen_Leaf_Rate_Positive = max(0, Nitrogen_Leaf_Rate)
    
        # Rate of N accumulation in root
        Nitrogen_Root_Accumulation = self.Switch_Function(Nitrogen_Total_Available + Nitrogen_Seed_Growth, NitrogenUptake - Nitrogen_Shoot_New - Nitrogen_Root_Available - LeafNitrogenLoss_Rate, NitrogenUptake - Nitrogen_Shoot_New - Nitrogen_Root_Available / self.Avoid_Zero_Division(Nitrogen_Total_Available) * (-Nitrogen_Seed_Growth) - LeafNitrogenLoss_Rate)
        RootNitrogen_Growth = self.Switch_Function(Nitrogen_Seed_Growth, Nitrogen_Root_Accumulation, NitrogenUptake - Nitrogen_Shoot_New - LeafNitrogenLoss_Rate)
        Nitrogen_Root_Rate = max(-self.Nitrogen_Root + 5E-8, RootNitrogen_Growth)
    
        self.Nitrogen_Root_Rate = Nitrogen_Root_Rate
        self.Nitrogen_Stem_Rate = Nitrogen_Stem_Rate
        self.Nitrogen_Leaf_Rate = Nitrogen_Leaf_Rate
        self.Nitrogen_Leaf_Rate_Positive = Nitrogen_Leaf_Rate_Positive
        self.Nitrogen_Seed_Rate = Nitrogen_Seed_Rate

    def Calculate_Leaf_Area_ChangeRate(self, LAI_Current, SpecificLeafNitrogen_BottomCurrent, LeafNitrogen_ExtinctionCoefficient):
        Rate_SpecificLeafNitrogen_Bottom = (SpecificLeafNitrogen_BottomCurrent - self.SpecificLeafNitrogen_Bottom) / self.Model_TimeStep
    
        # Rate of LAI driven by carbon supply
        Rate_LeafAreaIndex = self.Switch_Function(self.LeafWeight_Rate, max(-LAI_Current + 1E-5, self.SpecificLeafArea_Const * self.LeafWeight_Rate), self.SpecificLeafArea_Const * self.LeafWeight_Rate)
    
        # Adjusting LAI rate based on nitrogen during the juvenile phase
        if LAI_Current < 1 and self.DevelopmentStage < 0.5:
            Rate_LeafAreaIndex = (self.SpecificLeafNitrogen_Bottom * self.Nitrogen_Leaf_Rate - self.Nitrogen_Leaf * Rate_SpecificLeafNitrogen_Bottom) / self.SpecificLeafNitrogen_Bottom / (self.SpecificLeafNitrogen_Bottom + LeafNitrogen_ExtinctionCoefficient * self.Nitrogen_Leaf)
    
        self.Rate_LeafAreaIndex = Rate_LeafAreaIndex
        self.Rate_SpecificLeafNitrogen_Bottom = Rate_SpecificLeafNitrogen_Bottom


    
    
    def Calculate_Respiration_Rates(self, NitrogenUptake_Nitrate, NitrogenUptake_Ammonium, LastCarbonRoot):
        MaintenanceLeaf_DarkRespiration = 0.06 * (1 - self.Fraction_Shoot) * self.Photo_Assimilate
        MaintenanceStem_UndergroundStorage = 0.06 * 0.05 / 0.454 * self.YieldGrowth_Veg * self.Photo_Assimilate
        MaintenanceNitrogen_UptakeNitrate = 44 / 12 * 2.05 * NitrogenUptake_Nitrate
        MaintenanceNitrogen_UptakeAmmonium = 44.0 / 12.0 * 0.17 * NitrogenUptake_Ammonium
        Rate_RespirationMultiplier = (MaintenanceNitrogen_UptakeNitrate + MaintenanceNitrogen_UptakeAmmonium + MaintenanceStem_UndergroundStorage + MaintenanceLeaf_DarkRespiration - self.Respiration_Main) / self.Model_TimeStep
        GrowthRespiration = 44.0 / 12.0 * (((1.0 - self.YieldGrowth_Veg) / self.YieldGrowth_Veg * (self.CarbonSupply_Leaf + self.CarbonSupply_StemTotal + self.CarbonSupply_RootTotal + self.Carbon_Leaf_Loss + LastCarbonRoot)) + ((1.0 - self.YieldGrowth_Output) / self.YieldGrowth_Output * self.CarbonSupply_Seed))
        TotalRespiration = self.Respiration_Main + self.Respired_Carbon + GrowthRespiration + 44.0 / 12.0 * 0.06 * (self.CarbonReserve_Shoot + self.CarbonReserve_Root)
    
        self.MaintenanceLeaf_DarkRespiration = MaintenanceLeaf_DarkRespiration
        self.MaintenanceStem_UndergroundStorage = MaintenanceStem_UndergroundStorage
        self.MaintenanceNitrogen_UptakeNitrate = MaintenanceNitrogen_UptakeNitrate
        self.MaintenanceNitrogen_UptakeAmmonium = MaintenanceNitrogen_UptakeAmmonium
        self.Rate_RespirationMultiplier = Rate_RespirationMultiplier
        self.GrowthRespiration = GrowthRespiration
        self.TotalRespiration = TotalRespiration

    def Calculate_Carbon_Nitrogen_Returns(self, LastCarbonRoot, LeafNitrogenLoss_Rate, AvgSoilTemp, PhotosynthNitrogenLimitation_Sensitivity):
        LeafDetritus_Carbon = (self.Carbon_DeadLeaf - self.Carbon_DeadLeafShed) / 10. * (AvgSoilTemp - self.BaseTemperature_Development) / (self.OptimumTemperature_Development - self.BaseTemperature_Development)
        LitterCarbon_Total = LastCarbonRoot + LeafDetritus_Carbon
        LitterNitrogen_Total = LeafNitrogenLoss_Rate + LeafDetritus_Carbon / self.CarbonFrac_Veg * self.MinLeafN_Conc * PhotosynthNitrogenLimitation_Sensitivity
        NitrogenReturns_Total = self.LitterNitrogen_TotalPrevious + self.Switch_Function(self.DevelopmentStage - 2, 0, self.Nitrogen_Leaf + self.Nitrogen_Stem + self.Nitrogen_Root + self.NitrogenFixed_Reserve + (self.Carbon_DeadLeaf - self.Carbon_DeadLeafShed) / self.CarbonFrac_Veg * self.MinLeafN_Conc * (1 + PhotosynthNitrogenLimitation_Sensitivity) / 2)
    
        self.LeafDetritus_Carbon = LeafDetritus_Carbon
        self.LitterCarbon_Total = LitterCarbon_Total
        self.LitterNitrogen_Total = LitterNitrogen_Total
        self.NitrogenReturns_Total = NitrogenReturns_Total


    def Check_Carbon_Nitrogen_Balance(self, TotalNitrogenUptake):
        CarbonCheck_Input = self.Carbon_Total + self.Carbon_DeadLeaf + self.Carbon_DeadRoot - self.Carbon_InitialLeaf - self.Carbon_InitialRoot
        CarbonCheck_Balance = (CarbonCheck_Input - (self.TotalPhotosynthesis_Canopy - self.TotalRespiration) * 12.0 / 44.0) / self.Avoid_Zero_Division(CarbonCheck_Input) * 100
        NitrogenCheck_Input = self.Nitrogen_Total + self.Nitrogen_DeadLeaf + self.Nitrogen_DeadRoot - self.Nitrogen_InitialLeaf - self.Nitrogen_InitialRoot
        NitrogenCheck_Balance = (NitrogenCheck_Input - TotalNitrogenUptake) / self.Avoid_Zero_Division(TotalNitrogenUptake) * 100
    
        self.CarbonCheck_Balance = CarbonCheck_Balance
        self.NitrogenCheck_Balance = NitrogenCheck_Balance




        
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
    
    def REANOR(self,x1, x2):
        if x1 <= 0 and x2 <= 0:
            out = 1.0
        else:
            out = 0.0
        return out



    def Update_State_Variables(self, LeafNitrogenLoss_Rate):
        self.DevelopmentStage += self.DevelopmentRate
    
        self.CarbonTotal_DailyUpdate += self.TotalDaily_ThermalUnit
        
        self.Carbon_Leaf += self.CarbonSupply_Leaf
        self.Carbon_DeadLeaf += self.Carbon_Leaf_Loss
        self.Carbon_Stem += self.CarbonSupply_StemTotal
        self.Carbon_Seed += self.CarbonSupply_Seed
        
        self.Carbon_Root += self.Carbon_InitialRoot
        
        self.Carbon_DeadLeafShed += self.LeafDetritus_Carbon
        
        self.Nitrogen_Root += self.Nitrogen_Root_Rate
        self.Nitrogen_Stem += self.Nitrogen_Stem_Rate
        self.Nitrogen_Leaf += self.Nitrogen_Leaf_Rate
        self.Nitrogen_Seed += self.Nitrogen_Seed_Rate
        self.TotalNitrogen_Leaf += self.Nitrogen_Leaf_Rate_Positive
        self.Nitrogen_DeadLeaf += self.Nitrogen_Leaf_Loss
    
        self.CarbonReserve_Shoot += self.ReserveCarbon_Shoot
        self.CarbonReserve_Root += self.ReserveCarbon_Root
        self.NitrogenRemoval_EndSeedFill += self.NitrogenRemoval_EndSeedFill_Rate
        self.NitrogenRemoval_EarlyFlowering += self.NitrogenRemoval_EarlyFlowering_Rate
        self.DailyCarbonDemand_SeedReserve += self.ReserveDemand_SeedRate
        self.DailyCarbonDemand_StemReserve += self.ReserveDemand_StemRate
        
        self.SpecificLeafNitrogen_Bottom += self.Rate_SpecificLeafNitrogen_Bottom
        self.LAI_Current += self.Rate_LeafAreaIndex
        
        self.RespirationMultiplier += self.Rate_RespirationMultiplier
        self.NitrogenDemand_Previous += self.NitrogenDemand_Rate
        self.Nitrogen_Supply += self.NitrogenSupply_Rate
        self.NitrogenFixed_Total += self.Nitrogen_Fixed
        self.NitrogenFixed_Rate += self.NitrogenFixation_Rate
        self.DailyCarbonDemand_TotalPrev += self.RateDaily_CarbonDemand_TotalPrev
        
        self.Plant_Height += self.Rate_Height
        
        self.TotalRespiration += self.TotalRespiration
        self.LitterNitrogen_TotalPrevious += self.LitterNitrogen_Total
        self.Total_CanopyPhotosyntheticallyActiveRadiation += self.Canopy_PhotosyntheticallyActiveRadiation
        
        self.Total_ActualCanopy_Photosynthesis += self.Actual_CanopyPhotosynthesis
        self.Total_ActualCanopy_Transpiration += self.Actual_CanopyTranspiration
        
        self.Nitrogen_DeadRoot += LeafNitrogenLoss_Rate

        
 
 
 













