import math
import numpy as np


class Root:
    def __init__(self, carbon_fraction_veg, root_to_shoot_ratio, max_root_depth, root_nitrogen_content, carbon_reserves_roots, min_N_concentration_root, Model_TimeStep, Soil_Depth_1):
        self.carbon_fraction_veg = carbon_fraction_veg
        self.root_to_shoot_ratio = root_to_shoot_ratio
        self.max_root_depth = max_root_depth
        self.root_nitrogen_content = root_nitrogen_content
        self.carbon_reserves_roots = carbon_reserves_roots
        self.min_N_concentration_root = min_N_concentration_root
        self.Model_TimeStep = Model_TimeStep
        self.root_depth_initial = max(2.0, Soil_Depth_1)
        self.root_depth_current = self.root_depth_initial
        self.carbon_structural_roots = 0
        self.carbon_loss_rate_senescence = 0
        self.weight_roots_living = 0
        self.nitrogen_loss_rate_senescence = 0
        self.carbon_dead_roots = 0


        
    def Switch_Function(self,x, y1, y2):

        if x < 0:
            out = y1
        else:
            out = y2

        return out


    def calculate_root_senescence(self, carbon_structural_roots):
        kcrn = -np.log(0.05) / 6.3424 / self.carbon_fraction_veg / self.root_to_shoot_ratio / self.max_root_depth
        nitrogen_determined_csrt = 1 / kcrn * math.log(1.0 + kcrn * max(0.0, (self.root_nitrogen_content * self.carbon_fraction_veg - self.carbon_reserves_roots * self.min_N_concentration_root)) / self.min_N_concentration_root)
        carbon_loss_rate_senescence = max(min(carbon_structural_roots - 1.0e-4, carbon_structural_roots - min(nitrogen_determined_csrt, carbon_structural_roots)), 0.0) / self.delt
        weight_roots_living = carbon_loss_rate_senescence / self.carbon_fraction_veg
        nitrogen_loss_rate_senescence = weight_roots_living * self.min_N_concentration_root
        self.weight_roots_living = weight_roots_living
        self.nitrogen_loss_rate_senescence = nitrogen_loss_rate_senescence
        self.carbon_structural_roots = carbon_structural_roots
        self.carbon_loss_rate_senescence = carbon_loss_rate_senescence


    def calculate_rooting_depth(self, root_weight_total, weight_root_total, weight_root_top_down):
        extinction_coefficient = -np.log(0.05) / self.max_root_depth
        root_depth_growth_rate = self.switch_function(self.root_depth_current - self.max_root_depth, min((self.max_root_depth - self.root_depth_current) / self.delt, (root_weight_total + self.weight_roots_living) / (self.root_to_shoot_ratio + extinction_coefficient * (weight_root_total + weight_root_top_down))), 0)
        self.root_depth_growth_rate = root_depth_growth_rate

    def update_state_variables(self):
        self.carbon_dead_roots += self.carbon_loss_rate_senescence 
        self.root_depth_current += self.root_depth_growth_rate 
