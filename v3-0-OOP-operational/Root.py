import math
import numpy as np


class Root:
    def __init__(self, Soil_Layer_Property,Critical_root_weight_density, max_root_depth, Model_TimeStep):
        self.Critical_root_weight_density = Critical_root_weight_density
        self.max_root_depth = max_root_depth
        self.Model_TimeStep = Model_TimeStep
        self.nitrogen_determined_Root_Carbon = 0
        self.Root_carbon_loss_rate_senescence = 0
        self.weight_roots_living = 0
        self.Root_nitrogen_loss_rate_senescence = 0
        self.carbon_dead_roots = 0
        root_depth_ini = Soil_Layer_Property['rootdepth']*100 #cm
        self.Root_Depth = root_depth_ini

        
    def Switch_Function(self,x, y1, y2):

        if x < 0:
            out = y1
        else:
            out = y2

        return out


    def Calculate_Root_Senescence(self, Root_Carbon,CarbonFrac_Veg,MinRootN_Conc,ReserveRoot_Carbon,Nitrogen_Root):
        Extinction_coefficient_root_N  = -np.log(0.05) / 6.3424 / CarbonFrac_Veg / self.Critical_root_weight_density / self.max_root_depth
        nitrogen_determined_Root_Carbon = 1 / Extinction_coefficient_root_N * math.log(1.0 + Extinction_coefficient_root_N * max(0.0, (Nitrogen_Root * CarbonFrac_Veg - ReserveRoot_Carbon * MinRootN_Conc)) / MinRootN_Conc)
        Root_carbon_loss_rate_senescence = max(min(Root_Carbon - 1.0e-4, Root_Carbon - min(nitrogen_determined_Root_Carbon, Root_Carbon)), 0.0) / self.Model_TimeStep
        # print(Root_carbon_loss_rate_senescence)
        weight_roots_living = Root_carbon_loss_rate_senescence / CarbonFrac_Veg
        Root_nitrogen_loss_rate_senescence = weight_roots_living * MinRootN_Conc
        self.weight_roots_living = weight_roots_living
        self.Root_nitrogen_loss_rate_senescence = Root_nitrogen_loss_rate_senescence
        self.nitrogen_determined_Root_Carbon = nitrogen_determined_Root_Carbon
        self.Root_carbon_loss_rate_senescence = Root_carbon_loss_rate_senescence
        # print(nitrogen_determined_Root_Carbon , ReserveRoot_Carbon )


    def Calculate_Rooting_Depth(self, RootWeight_Rate, LiveRoot_Dry_Weight, DeadRoot_Dry_Weight):
        extinction_coefficient = -np.log(0.05) / self.max_root_depth
        root_depth_growth_rate = self.Switch_Function(self.Root_Depth - self.max_root_depth, min((self.max_root_depth - self.Root_Depth) / self.Model_TimeStep, (RootWeight_Rate + self.weight_roots_living) / (self.Critical_root_weight_density + extinction_coefficient * (LiveRoot_Dry_Weight + DeadRoot_Dry_Weight))), 0)
        self.root_depth_growth_rate = root_depth_growth_rate 
        # print(root_depth_growth_rate)
    def Update_State_Variables(self):
        self.carbon_dead_roots += self.Root_carbon_loss_rate_senescence 
        self.Root_Depth += self.root_depth_growth_rate 
