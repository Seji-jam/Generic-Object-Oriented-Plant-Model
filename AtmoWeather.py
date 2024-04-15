import numpy as np
import pandas as pd



class WeatherMeteo:
    def __init__(self, filename, lat, Simulation_Start_DOY):
        self.filename = filename
        self.lat = lat
        self.Sun_Angle_Inclination = -2
        self.Start_Time = Simulation_Start_DOY
        self.rad = np.pi / 180

    def meteo_calculations(self, doy):
        if not -67 <= self.lat <= 67:
            raise ValueError('ERROR IN meteo_calculations: LAT out of bounds [-67, 67]')
        
        dec = np.arcsin(np.sin(23.45 * self.rad) * np.cos(2 * np.pi * (doy + 10) / 365))
        Sin_Solar_Declination = np.sin(self.rad * self.lat) * np.sin(dec)
        Cos_Solar_Declination = np.cos(self.rad * self.lat) * np.cos(dec)
        angle_factor = Sin_Solar_Declination / Cos_Solar_Declination

        Day_Length = 12.0 * (1 + 2 * np.arcsin(angle_factor) / np.pi)
        Photoperiod_Day_Length = 12.0 * (1 + 2 * np.arcsin((-np.sin(self.Sun_Angle_Inclination * self.rad) + Sin_Solar_Declination) / Cos_Solar_Declination) / np.pi)
        Daily_Sin_Beam_Exposure = 3600 * (Day_Length * (Sin_Solar_Declination + 0.4 * (Sin_Solar_Declination**2 + Cos_Solar_Declination**2 * 0.5)) +
                         12.0 * Cos_Solar_Declination * (2.0 + 3.0 * 0.4 * Sin_Solar_Declination) * np.sqrt(1.0 - angle_factor**2) / np.pi)
        Solar_Constant = 1367 * (1 + 0.033 * np.cos(2 * np.pi * (doy - 10) / 365))

        return Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Photoperiod_Day_Length, Daily_Sin_Beam_Exposure

    def process_data(self):
        nld = pd.read_excel(self.filename, sheet_name='Sheet1').values
        data_results = []

        for i, row in enumerate(nld):
            Station, Year, Doy, Solar_Radiation, Min_Temp, Max_Temp, Vapour_Pressure, Wind_Speed, Rain = row
            Solar_Radiation *= 1000  # Convert to the correct unit if necessary
            time = Doy

            dfs = time - self.Start_Time + 1

            Wind_Speed = max(0.1, Wind_Speed)

            Day_time_Temp = 0.29 * Min_Temp + 0.71 * Max_Temp
            Night_Time_Temp = 0.71 * Min_Temp + 0.29 * Max_Temp
            


            Solar_Constant, Sin_Solar_Declination, Cos_Solar_Declination, Day_Length, Photoperiod_Day_Length, Daily_Sin_Beam_Exposure=self.meteo_calculations(Doy)


            # Store the calculated values in a dictionary
            iteration_results = {
                'Station': Station, 'Year': Year, 'Doy': Doy,'dfs':dfs,
                'Max_Temp': Max_Temp, 'Min_Temp': Min_Temp, 'Solar_Radiation': Solar_Radiation,'Rain':Rain, 'Vapour_Pressure': Vapour_Pressure,
                'Wind_Speed': Wind_Speed, 'Day_time_Temp': Day_time_Temp, 'Night_Time_Temp': Night_Time_Temp,
                # Add meteorological calculations
                'Solar_Constant': Solar_Constant, 'Sin_Solar_Declination':Sin_Solar_Declination,'cosld':Cos_Solar_Declination,'Day_Length':Day_Length,
                'Photoperiod_Day_Length':Photoperiod_Day_Length,'Daily_Sin_Beam_Exposure':Daily_Sin_Beam_Exposure
            }

            data_results.append(iteration_results)

        return data_results



















