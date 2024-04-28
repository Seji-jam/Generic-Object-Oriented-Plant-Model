import os

class OutputManager:
    def __init__(self, filename='MODEL_OUTPUTS.csv'):
        self.current_dir = os.path.dirname(os.path.abspath('__file__'))
        self.output_file = os.path.join(self.current_dir, filename)
        self.header_written = False

    def write_header(self,header):
        if not self.header_written:
            with open(self.output_file, "w") as f:
                f.write(header + '\n') 
            self.header_written = True

    def append_data(self, data):
        with open(self.output_file, "a") as f:
            f.write(data )

    
    def format_header():
        header="Year,doy,Development_Stage,LAI,Cumulative_Thermal_Unit,\
LiveRoot_Dry_Weight,Shoot_Dry_Weight,weight_roots_living,\
Leaf_Carbon,Stem_Carbon,Seed_Carbon,Root_Carbon,\
Nitrogen_Leaf,Nitrogen_Stem,Nitrogen_Seed,Nitrogen_Root,\
Nitrogen_uptake,\
Respiration_Uptakes,Maintenance_Respiration,Cumulative_Respiration,\
Potential_Canopy_Photosynthesis,Actual_Canopy_Photosynthesis,\
Potential_Canopy_Temp,Actual_Canopy_Temp,\
potential_evap_daily,Actual_Daily_Evaporation,\
Potential_Canopy_Transpiration,Actual_Canopy_Transpiration,\
Plant_Height,\
available_soluble_N,Root_Depth,average_soil_temperature,\
Current_Soil_Moisture_Top_Layer"
        return header

    
    def format_data(day_data, Canopy_object,Leaf_object, Root_object, Soil_object):
        return f"{day_data['Year']},{day_data['Doy']},{float(Canopy_object.Development_Stage):.2f},{float(Leaf_object.Leaf_area_output['Leaf_Area_Index']):.3f},\
{float(Canopy_object.Cumulative_Thermal_Unit):.2f},\
{float(Canopy_object.LiveRoot_Dry_Weight):.3f},{float(Canopy_object.Shoot_Dry_Weight):.3f},{float(Root_object.weight_roots_living):.3f},\
{float(Canopy_object.Leaf_Carbon):.3f},{float(Canopy_object.Stem_Carbon):.3f},{float(Canopy_object.Seed_Carbon):.3f},{float(Canopy_object.Root_Carbon):.3f},\
{float(Canopy_object.Nitrogen_Leaf):.3f},{float(Canopy_object.Nitrogen_Stem):.3f},{float(Canopy_object.Nitrogen_Seed):.3f},{float(Canopy_object.Nitrogen_Root):.3f},\
{float(Soil_object.Nitrogen_uptake):.3f},\
{float(Canopy_object.Respiration_Uptakes):.3f},{float(Canopy_object.Maintenance_Respiration):.3f},{float(Canopy_object.Cumulative_Respiration):.3f},\
{float(Canopy_object.Potential_Canopy_Photosynthesis):.5f},{float(Canopy_object.Actual_Canopy_Photosynthesis):.5f},\
{float(Canopy_object.Potential_Canopy_Temp):.5f},{float(Canopy_object.Actual_Canopy_Temp):.5f},\
{float(Soil_object.potential_evap_daily):.5f},{float(Soil_object.Actual_Daily_Evaporation):.5f},\
{float(Canopy_object.Potential_Canopy_Transpiration):.5f},{float(Canopy_object.Actual_Canopy_Transpiration):.5f},\
{float(Canopy_object.Plant_Height):.3f},\
{float(Soil_object.available_soluble_N):.3f},{float(Root_object.root_depth_current):.3f},{float(Soil_object.average_soil_temperature):.3f},\
{float(Soil_object.Current_Soil_Moisture_Top_Layer):.4f}, \n"
    
    
    
    
    
    
    
    
    