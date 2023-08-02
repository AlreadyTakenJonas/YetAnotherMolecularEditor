# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 12:24:48 2023

@author: jonas
"""

import pandas as pd
import re
from . import PINT_UNIT_REGISTRY


# Create periodic table of the elements (contains information about elements, used to create new atoms)
PERIODICTABLE = pd.read_csv("data/periodicTable.csv")

def parseChemicalSpeciesNotation(specie:str):
    """
    TODO ... Parse a string of a chemical specie to get elementSymbol, mass, charge and other information.

    Details on parsing the atom-type/species:
        Write the mass (optional) in front of the element symbol and the charge (optional) behind the element symbol.
        The mass and charge need to be integer numbers. The charges +1 and -1 can be shortened to '+' and '-'.
        E.g. species="2H+" will add an deuterium ion with one positive charge to the molecule.
        E.g. species="56Fe3+" will add Fe with mass 56u and a charge of +3 to the molecule.
        E.g. species="C" will add one carbon atom without charge and a mass of 12u to the molecule.
        E.g. species="O2-" will add a double negatively charge oxygen atom of mass 16u to the atom.
        Be aware: This function does not check if the given mass or charge makes sense! It only checks the element symbol.

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
    
    # Check the specie-string for illegal characters: anything, but letters, digits, "+" and "+".
    assert re.search("[^( A-Z a-z \d + - )]", specie) == None, f"Illegal specie-string! Unallowed character in '{specie}'. Only digits, letters, '+' and '-' allowed!"
    assert re.search("\d$", specie) == None, "Charge ambigous! Add '+' or '-' at the end of the charge!"

    
    # Extract element symbol from species (e.g. Fe or H)
    elementSymbol = re.findall("[A-Z][a-z]?", specie)[0]
    
    # Lookup the element symbol in the periodic table and get element properties
    element = PERIODICTABLE.loc[PERIODICTABLE.symbol.str.fullmatch(elementSymbol)]
    if len(element) == 0: raise  ValueError(f"I don't know the element {elementSymbol}!")
    
    # Get the mass from species string if given. Get default from periodic table
    # Convert mass to atomic mass units with pint
    mass = re.findall("^\d+", specie)
    if len(mass) == 0: mass = int(element.mass) * PINT_UNIT_REGISTRY.amu
    else:              mass = int(mass) * PINT_UNIT_REGISTRY.amu
    
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


def isIterable(obj):
    try:
        _ = iter(obj)
    except TypeError:
        iterable = False
    else:
        iterable = True
    return iterable