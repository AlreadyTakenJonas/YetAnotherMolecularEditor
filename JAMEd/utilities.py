# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 12:24:48 2023

@author: jonas
"""

import pandas as pd
import re


# Create periodic table of the elements (contains information about elements, used to create new atoms)
PERIODICTABLE = pd.read_csv("data/periodicTable.csv")

def parseChemicalSpeciesNotation(specie:str):
    """
    TODO ... Parse a string of a chemical specie to get elementSymbol, mass, charge and other information.

    Parameters
    ----------
    specie : str
        DESCRIPTION.

    Raises
    ------
    
        DESCRIPTION.

    Returns
    -------
    elementSymbol : TYPE
        DESCRIPTION.
    element : TYPE
        DESCRIPTION.
    mass : TYPE
        DESCRIPTION.
    charge : TYPE
        DESCRIPTION.

    """
    
    # Check input type
    assert isinstance(specie, str), "Specie must be a string!"
    
    # Extract element symbol from species (e.g. Fe or H)
    elementSymbol = re.findall("[A-Z][a-z]?", specie)[0]
    
    # Lookup the element symbol in the periodic table and get element properties
    element = PERIODICTABLE.loc[PERIODICTABLE.symbol.str.fullmatch(elementSymbol)]
    if len(element) == 0: raise  ValueError(f"I don't know the element {elementSymbol}!")
    
    # Get the mass from species string if given. Get default from periodic table
    mass = re.findall("^\d+", specie)
    if len(mass) == 0: mass = int(element.mass)
    else:              mass = int(mass)
    
    # Get charge from species string if given. Default is neutral charge (0).
    charge = re.findall("\d*[+,-]", specie)
    if len(charge) == 0: charge = 0
    else:
        charge = charge[1]
        number = charge[:-1]
        sign   = charge[-1]
        charge = int(f"{sign}{number}")
    
    # Return element properties in periodic table, mass and charge
    return element, mass, charge