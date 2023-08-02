# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:03:50 2023

@author: jonas
"""

import pandas as pd
from pathlib import Path
from . import PINT_UNIT_REGISTRY

class FileFormatParser:
    
    def __init__(self):
        self._atomTableTemplate = pd.DataFrame(columns=["elementSymbol", "cartesian", "fixedInternalCoordinates"])
        self._atomTable = self._atomTableTemplate
    
    #
    #   DEFINE INTERACTION WITH DATA FRAME FOR STORING ATOM COORDINATES.
    #   ANY FILE FORMAT CONVERSION WILL BE CONVERTED INTO/FROM THIS DATA FRAME.
    #
    @property
    def atomTable(self):
        return self._atomTable
    @atomTable.setter
    def atomTable(self, value):
        # Check if input has valid format.
        assert isinstance(value, pd.DataFrame), "atomTable must be a pandas DataFrame!"
        assert "elementSymbol" in value, "atomTable must contain 'elementSymbol' column!"
        assert "cartesian" in value, "atomTable must contain 'cartesian' column!"

        # Empty existing atom table
        self._atomTable = self._atomTableTemplate
        # Add new input to now empty atom table. Keep all columns in atom table. Don't create new columns in atom table.
        commonColumns = set(self._atomTable.columns) & set(value.columns)        
        self._atomTable[commonColumns] = value[commonColumns]
    
    #
    #   DEFINE GETTER & SETTER FOR PASING AND CREATING XYZ-FILES
    #
    @property
    def xyz(self):
        # Convert table of atom coordinates into string in xyz file format.
        output = str(len(self.atomTable)) + "\n"
        for index, element in self.atomTable.iterrows():
            output += element.elementSymbol + "\t" + element.cartesian[0] + "\t" + element.cartesian[1] + "\t" + element.cartesian[2] + "\n"
        return output
    @xyz.setter
    def xyz(self, value):
        # If Path object was passed. Read the file and use the files content as input to the function.
        if isinstance(value, Path):
            assert value.exists(), f"xyz file {value.resolve()} does not exist!"
            assert value.is_file(), f"Can't read {value.resolve()}. It's not a file!"
            value = value.read_text()
        
        assert isinstance(value, str), "xyz file content must be a string!"
        
        # Split the file into lines
        lines = value.splitlines()
        # Remove first line and store it in seperate variable. First line contains number of atoms in the xyz file
        expectedAtomCount = int(lines.pop(0))
        # Check if xyz file is complete.
        assert len(lines) == expectedAtomCount, "Can't parse xyz file! Number of expected atoms does not match number of present atoms!"
        
        # Read element symbol and cartesian coordinates from xyz file. Store coordinates as angstrom.
        elementSymbol = [ line.split()[0] for line in lines ]
        x             = [ float(line.split()[1]) for line in lines ] * PINT_UNIT_REGISTRY.angstrom
        y             = [ float(line.split()[2]) for line in lines ] * PINT_UNIT_REGISTRY.angstrom
        z             = [ float(line.split()[3]) for line in lines ] * PINT_UNIT_REGISTRY.angstrom
                
        # Update table of atom coordinates.
        self.atomTable = pd.DataFrame({"elementSymbol":elementSymbol, "cartesian":zip(x,y,z)})