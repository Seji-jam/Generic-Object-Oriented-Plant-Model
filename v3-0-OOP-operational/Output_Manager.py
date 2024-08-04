import os

class OutputManager:
    def __init__(self, filename='MODEL_OUTPUTS.csv'):
        self.current_dir = os.path.dirname(os.path.abspath('__file__'))
        self.output_file = os.path.join(self.current_dir, filename)
        self.header_written = False

    def write_header(self, header):
        if not self.header_written:
            with open(self.output_file, "w") as f:
                f.write(header + '\n') 
            self.header_written = True

    def append_data(self, data):
        with open(self.output_file, "a") as f:
            f.write(data)

    def format_header(self, Soil_object):
        base_header = "Year,doy,Development_Stage,LAI,Cumulative_Thermal_Unit,\
LiveRoot_Dry_Weight,Shoot_Dry_Weight,weight_roots_living,\
Fraction_Leaf_Carbon,Leaf_Carbon,\
Fraction_Seed_Carbon,Seed_Carbon,\
Fraction_Stem_Carbon,Stem_Carbon,\
Total_CarbonDemand_StemGrowth,DailyCarbon_Supply_Stem,Carbon_Flow_to_Stem,\
Root_Carbon,\
Nitrogen_Leaf,Nitrogen_Stem,Nitrogen_Seed,Nitrogen_Root,\
Nitrogen_uptake,\
Respiration_Uptakes,Maintenance_Respiration,Cumulative_Respiration,\
Potential_Canopy_Photosynthesis,Actual_Canopy_Photosynthesis,\
Potential_Canopy_Temp,Actual_Canopy_Temp,\
Actual_Daily_Evaporation,\
Potential_Canopy_Transpiration,Actual_Canopy_Transpiration,\
Plant_Height,\
Root_Depth,average_soil_temperature,\
soil_moisture_root_zone,soil_moisture_below_root_zone"
        
        # Add soil moisture columns based on the number of layers
        num_layers = Soil_object.number_of_layers
        for i in range(1, num_layers + 1):
            base_header += f",soil_moisture_{i}"
        
        return base_header

    def format_data(self, day_data, Canopy_object, Leaf_object, Root_object, Soil_object):
        base_data = f"{day_data['Year']},{day_data['Doy']},{float(Canopy_object.Development_Stage):.3f},{float(Leaf_object.Leaf_area_output['Leaf_Area_Index']):.4f},\
{float(Canopy_object.Cumulative_Thermal_Unit):.2f},\
{float(Canopy_object.LiveRoot_Dry_Weight):.3f},{float(Canopy_object.Shoot_Dry_Weight):.3f},{float(Root_object.weight_roots_living):.3f},\
{float(Canopy_object.Fraction_Leaf_Carbon):.3f},{float(Canopy_object.Leaf_Carbon):.3f},\
{float(Canopy_object.Fraction_Seed_Carbon):.3f},{float(Canopy_object.Seed_Carbon):.3f},\
{float(Canopy_object.Fraction_Stem_Carbon):.3f},{float(Canopy_object.Stem_Carbon):.3f},\
{float(Canopy_object.Total_CarbonDemand_StemGrowth):.3f},{float(Canopy_object.DailyCarbon_Supply_Stem):.3f},{float(Canopy_object.Carbon_Flow_to_Stem):.3f},\
{float(Canopy_object.Root_Carbon):.3f},\
{float(Canopy_object.Nitrogen_Leaf):.3f},{float(Canopy_object.Nitrogen_Stem):.3f},{float(Canopy_object.Nitrogen_Seed):.3f},{float(Canopy_object.Nitrogen_Root):.3f},\
{float(Soil_object.total_nitrogen_uptake):.3f},\
{float(Canopy_object.Respiration_Uptakes):.3f},{float(Canopy_object.Maintenance_Respiration):.3f},{float(Canopy_object.Cumulative_Respiration):.3f},\
{float(Canopy_object.Potential_Canopy_Photosynthesis):.5f},{float(Canopy_object.Actual_Canopy_Photosynthesis):.5f},\
{float(Canopy_object.Potential_Canopy_Temp):.5f},{float(Canopy_object.Actual_Canopy_Temp):.5f},\
{float(Soil_object.Actual_Daily_Evaporation):.5f},\
{float(Canopy_object.Potential_Canopy_Transpiration):.5f},{float(Canopy_object.Actual_Canopy_Transpiration):.5f},\
{float(Canopy_object.Plant_Height):.3f},\
{float(Root_object.Root_Depth):.3f},{float(Soil_object.average_soil_temperature):.3f},\
{float(Soil_object.average_root_zone_moisture):.4f},{float(Soil_object.average_below_root_zone_moisture):.4f}"
        
        # Add soil moisture values for each layer
        # layer = Soil_object.Soil_Layer_Property['layers']
        for layer in Soil_object.Soil_Layer_Property['layers']:
            value=layer['Soil_Moisture']
            base_data += f",{value:.4f}"
        
        return base_data + "\n"











