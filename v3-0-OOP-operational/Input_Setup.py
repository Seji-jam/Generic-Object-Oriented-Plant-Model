import pandas as pd
import os



class Import_Data:
    def __init__(self, file_path):
        self.file_path = file_path
        self.load_data()

    def load_data(self):
        """Loads data from an Excel file and sets attributes based on the file contents."""
        data = pd.read_excel(self.file_path)
        for index, row in data.iterrows():
            setattr(self, row['Parameter'], row['value'])

    def append_data(self, other):
        """Appends data from another CropData instance."""
        for key, value in other.__dict__.items():
            if key not in self.__dict__:  # This prevents overwriting existing attributes
                setattr(self, key, value)