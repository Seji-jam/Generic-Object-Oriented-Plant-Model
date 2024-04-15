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
            # Setting each parameter as an attribute of the class instance
            setattr(self, row['Parameter'], row['value'])