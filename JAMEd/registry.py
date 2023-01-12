#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:29:26 2022

@author: jonas
"""

import pandas as pd
from pathlib import Path
import logging
import tempfile
import yaml
from zipfile import ZipFile
from typing import Union
from .version import __version__
from packaging import version
from pint import Quantity
from . import PINT_UNIT_REGISTRY
import re

# Create logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

#   PROPERTY LISTS
#
# Define properties saved in the tables of the registry
#
# What will be saved about every atom of a molecule?
ATOM_PROPERTY_LIST = [
    "moleculeID",   # To which molecule does the atom belong (unique identifier for molecule)?
    "atomID",       # Unique ID
    "label",        # Descriptive label
    "element",      # How many protons does the atom has?
    "mass",         # How much mass does the atom has (not the exact mass, rounded to integer)?
    "charge",       # How much charge does the atom has? 
    "cartesian",    # Where is the atom in 3D space (tupel of three coordinates)?
    "fixedInternalCoordinates", # Is the atom fixed in space (internal coordinates)?
    "color",        # What color does the atom have when shown on the screen?
    "radius"        # How big is the atom?
]
# What will be saved about every molecule?
MOLECULE_PROPERTY_LIST = [
    "moleculeID",       # Unique ID
    "name",             # Descriptive name
    "centerOfInertia",  # Coordinates of center of inertia in 3D space (tupel of three coordinates)?
    "symmetry"          # Molecules point group
]
# What will be saved about every bond in a molecule?
BOND_PROPERTY_LIST = [
    "bondID",       # Unique ID of bond
    "atomID1",      # Unique ID of first atom of bond
    "atomID2",      # Unique ID of second atom of bond
    "order"         # Order of the bond
]
# Which settings will be saved in the state registry?
DEFAULT_SETTINGS = {}
# Save meta information. Used to make sure loading files written by different version runs smoothly.
META_DATA = {
    "fileLocation": "",    # From which file was the StateRegistry read from?
    "version": __version__ # Current version of this package
}
SAVE_TO_FILENAME = {
    "atom"      : "atoms.csv",
    "bond"      : "bonds.csv",
    "meta"      : "meta.yaml", 
    "setting"   : "settings.yaml",
    "molecule"  : "molecules.csv"
}


class StateRegistry:
    
    
    
    def __init__(self):
        """
        Initialise empty State registry with default settings.

        Returns
        -------
        None.

        """
        # Create tables to store information on molecules and atoms
        self._atomTable     = pd.DataFrame(columns = ATOM_PROPERTY_LIST    )
        self._moleculeTable = pd.DataFrame(columns = MOLECULE_PROPERTY_LIST)
        self._bondTable     = pd.DataFrame(columns = BOND_PROPERTY_LIST    )
        # Create dictionary to save settings
        self._settings      = DEFAULT_SETTINGS
        # Create dictionary to save meta data
        self._metadata      = META_DATA
        # Create periodic table of the elements (contains information about elements, used to create new atoms)
        self._periodicTable = pd.read_csv("data/periodicTable.csv")
        
        
    def addMolecule(self, amount:int=1) -> list[int]:
        """
        Create new entries in the molecule registry.

        Parameters
        ----------
        amount : int, optional
            Number of atoms to create. The default is 1.

        Returns
        -------
        newMoleculeID : list of int
            The ids of all created molecules.

        """
        # Create a new ID for the molecule.
        # Set ID to 0 if there are no molecules in registry yet
        if len(self._moleculeTable) == 0: newMoleculeID = 0
        # Get the highest ID already taken and add 1 to get the new ID
        else:                             newMoleculeID = max(self._moleculeTable.moleculeID) + 1
        # Generate list of molecule IDs. Can generatet multiple IDs to create multiple molecules at once.
        newMoleculeID = [ newMoleculeID + i for i in range(amount) ]
        
        # Create a new row to insert into the self._moleculeTable
        newMolecule = pd.DataFrame({"moleculeID": newMoleculeID,
                                    "name": [ f"Molecule #{ID}" for ID in newMoleculeID ]})
        # Insert the new row into the molecule table
        self._moleculeTable = pd.concat([self._moleculeTable, newMolecule], ignore_index = True)
        
        # Return IDs of generated molecules.
        return newMoleculeID
    
    def destroyMolecule(self, moleculeID:list[int], recursive:bool=True):
        """
        Destroy molecule. TODO: DOCSTRING ...

        Parameters
        ----------
        moleculeID : list[int]
            DESCRIPTION.
        recursive : bool, optional
            DESCRIPTION. The default is True.

        Raises
        ------
        this
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # Check if the moleculeID is iterable.
        # Save the result in the variable iterable.        
        try:
            _ = iter(moleculeID)
        except TypeError:
            iterable = False
        else:
            iterable = True
            # moleculeID is iterable. Check if we want to continue. This option is needed to avoid recalling this function on nested iterables recursively.
            assert recursive, "Cannot destroy this iterable with moleculeIDs! Pass recursive=True to destroy every molecule in the iterable. Don't pass an iterable of iterables. It will always raise this exception."
        
        # Is atomID a iterable? If yes loop over it and call this function on each element.
        if iterable:
            for molecule in moleculeID:
                self.destroyMolecule(molecule, recursive=False)
        
        # If moleculeID is not iterable, delete the molecule
        else:
            # Is the moleculeID an integer?
            assert isinstance(moleculeID, int), "Destroying molecule failed! MoleculeID must be integer or itereable of integer!"
            
            # Get a list of all atomIDs belonging to the deleted molecule
            garbageAtoms = self._atomTable[ self._atomTable.moleculeID == moleculeID ].atomID
            # Destroy all atoms belonging to the molecule. After all atoms are destroyed, the molecule will be removed from the _moleculeTable.
            self.destroyAtom(garbageAtoms)
            
        
    def addAtom(self, moleculeID:int, specie:str, coordinates:Quantity):
        
        #
        # TODO: COMMENT, DOCSTRING
        #
        
        assert isinstance(moleculeID, int), f"Molecule IDs must be int not {type(moleculeID)}!"
        assert moleculeID in self._moleculeTable.moleculeID, f"Molecule with ID {moleculeID} not found!"
        assert isinstance(coordinates, Quantity), "Coordinates may not be unitless! Use pints unit registry defined in the packages __init__.py!"
        assert len(coordinates) == 3, f"You need to pass three coordinates! {len(coordinates)} was given!"
        
        elementSymbol = re.findall("[A-Z][a-z]?", specie)[0]
        
        element = self._periodicTable.loc[self._periodicTable.symbol.str.fullmatch(elementSymbol)]
        if len(element) == 0: raise  ValueError(f"I don't know the element {elementSymbol}!")
        
        mass = re.findall("^\d+", specie)
        if len(mass) == 0: mass = int(element.mass)
        else:              mass = int(mass)
        charge = re.findall("\d*[+,-]", specie)
        if len(charge) == 0: charge = 0
        else:
            charge = charge[1]
            number = charge[:-1]
            sign   = charge[-1]
            charge = int(f"{sign}{number}")
        
        # Convert coordinates to angstrom and then to float numbers
        coordinates = tuple([ float(coord.to("angstrom") / PINT_UNIT_REGISTRY.angstrom) for coord in coordinates ])
        
        if len(self._atomTable) == 0: newAtomID = 0
        else:                         newAtomID = max(self._atomTable.atomID) + 1
        newAtom = pd.DataFrame({"moleculeID": [moleculeID],
                                "atomID": [newAtomID],
                                "cartesian": [coordinates],
                                "label": [elementSymbol],
                                "element": element.element,
                                "charge": [charge],
                                "mass": [mass],
                                "color": element.color,
                                "radius": element.radius,
                                "fixedInternalCoordinates": [False]
                               })
        
        self._atomTable = pd.concat([self._atomTable, newAtom], ignore_index=True)
        
        return newAtomID
    
    def destroyAtom(self, atomID:list[int], recursive:bool = True):
        """
        Remove atoms from a molecule. TODO: DOCSTRING.

        Parameters
        ----------
        atomID : list[int]
            DESCRIPTION.
        recursive : bool, optional
            DESCRIPTION. The default is True.

        Raises
        ------
        this
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # Check if the atomID is iterable.
        # Save the result in the variable iterable.        
        try:
            _ = iter(atomID)
        except TypeError:
            iterable = False
        else:
            iterable = True
            # atomID is iterable. Check if we want to continue. This option is needed to avoid recalling this function on nested iterables recursively.
            assert recursive, "Cannot destroy this iterable with atomIDs! Pass recursive=True to destroy every atom in the iterable. Don't pass an iterable of iterables. It will always raise this exception."
        
        # Is atomID a iterable? If yes loop over it and call this function on each element.
        if iterable:
            for atom in atomID:
                self.destroyAtom(atom, recursive=False)
        
        # If atomID is not iterable, delete the atom from the molecule
        else:
            # Is the atomID an integer?
            assert isinstance(atomID, int), "Destroying atom failed! AtomID must be integer or itereable of integer!"
            
            # Get the moleculeID the deleted atom belongs to. Used to delete the molecule, if all its atoms are deleted.
            moleculeID = self._atomTable.loc[ self._atomTable["atomID"] == atomID ]["moleculeID"]
            # Remove the atom from the registry.
            self._atomTable = self._atomTable[self._atomTable.atomID != atomID]
            # Count the number of atoms left in the molecule.
            numberOfAtomsInMolecule = len(self._atomTable[self._atomTable.moleculeID == moleculeID].index)
            
            # Get IDs of bonds, that belong to destroyed atoms. Used to remove these bonds from registry.
            garbageBonds = self._bondTable[ self._bondTable.atomID1 == atomID or self._bondTable.atomID2 == atomID ]["bondID"]
            # Remove bonds to deleted atoms from registry.
            self.destroyBond(garbageBonds)
            
            # Delete the molecule, if it does not contain any atoms.
            if numberOfAtomsInMolecule == 0:
                self.destroyMolecule(moleculeID)
        
    def replaceAtom(self, atomID:int, specie:str):
        raise NotImplementedError("Replacing Atoms not implemented yet.")
        
    def addBond(self, atomID1:int, atomID2:int, order:int=1):
        pass
    
    def destroyBond(self, bondID:list[int], recursive:bool = True):
        """
        Remove a bond. TODO: DOCSTRING ...

        Parameters
        ----------
        bondID : list[int]
            DESCRIPTION.
        recursive : bool, optional
            DESCRIPTION. The default is True.

        Raises
        ------
        this
            DESCRIPTION.

        Returns
        -------
        None.

        """      
        
        # Check if the bondID is iterable.
        # Save the result in the variable iterable.        
        try:
            _ = iter(bondID)
        except TypeError:
            iterable = False
        else:
            iterable = True
            # bondID is iterable. Check if we want to continue. This option is needed to avoid recalling this function on nested iterables recursively.
            assert recursive, "Cannot destroy this iterable with bondIDs! Pass recursive=True to destroy every atom in the iterable. Don't pass an iterable of iterables. It will always raise this exception."
        
        # Is bondID a iterable? If yes loop over it and call this function on each element.
        if iterable:
            for bond in bondID:
                self.destroyBond(bond, recursive=False)
        # If bondID is not iterable, remove the entry with the given bondID
        else:
            # Check type of bondID
            assert isinstance(bondID, int), "Destroying bond failed! BondID must be integer or itereable of integer!"
            # Remove bond from registry.
            self._bondTable = self._bondTable[ self._bondTable.bondID != bondID ]
        
    def save(self, path:Union[str, Path], override:bool=False):
        """
        Save the instance of this class to a zip archive. Can be reinstantiated with load method.

        Parameters
        ----------
        path : Union[str, Path]
            Name of the file.
        override : bool, optional
            Override preexisting file with the same name? The default is False.

        Raises
        ------
        FileExistsError
            Can not save to file, because file with the same name already exists. Run with override=True to save anyways.
        TypeError
            Path was passed the wrong type.

        Returns
        -------
        None.

        """
        # Was the path passed as string or Path? Convert it to Path object if it is a string.
        if isinstance(path, str): path = Path(path)
        # Is the path now the correct type?
        if not isinstance(path, Path):
            log.error(f"Unexpected type of path parameter. Expected string or pathlib.Path not {type(path)}!")
            raise TypeError(f"Unexpected type of path parameter. Expected string or pathlib.Path not {type(path)}!")
        log.info(f"Saving molecule to {path.resolve()} ...")
        
        # Does the file already exist?
        if path.exists():
            log.info(f"{path.resolve()} already exists.")
            # May the existing file be overritten?
            if override == False:
                log.error(f"Can't save molecules to {path.resolve()}, because it already exists.")
                raise FileExistsError(f"Can't save molecules to {path.resolve()}, because it already exists.")
        
        # Is the suffix right? If not it is alright.
        if path.suffix != ".jam": log.warning(f"Unexpected file suffix '{path.suffix}'. Why not save it as .jam? Saving it anyways.")
        
        # Write the state of the StateRegistry to a temporary directory and move them into a zip archive.
        with tempfile.TemporaryDirectory() as tmp:
            # Write data on atoms, molecules and bonds to csv files
            self._atomTable.to_csv( tmp+"/"+SAVE_TO_FILENAME["atom"] )
            self._moleculeTable.to_csv( tmp+"/"+SAVE_TO_FILENAME["molecule"] )        
            self._bondTable.to_csv( tmp+"/"+SAVE_TO_FILENAME["bond"] )
            # Dump the settings dictionary to yaml file
            dumpedSettings = yaml.safe_dump(self._settings)
            Path(tmp+"/"+SAVE_TO_FILENAME["setting"]).write_text(dumpedSettings)
            # Dump the meta data dictionary to yaml file
            dumpedMetadata = yaml.sage_dump(self._metadata)
            Path(tmp+"/"+SAVE_TO_FILENAME["meta"]).write_text(dumpedMetadata)
            
            # Move all files in temporary director´y to zip archive
            with ZipFile(path.resolve(), "w") as archive:
                # Loop over all files in the temporary directory
                for file in Path(tmp).glob("*"):
                    # Write the temporary file to the zip archive
                    # Use just the name of the file in the archive, not the whole path (arcname=file.name)
                    archive.write(file.resolve(), arcname=file.name)   
    
    @classmethod
    def load(cls, path:Union[str, Path], forceLoad:bool=False) -> "StateRegistry":
        """
        Load StateRegistry saved to file previously by save method. Load position of atoms, bonds and metadata on molecules and settings.

        Parameters
        ----------
        path : Union[str, Path]
            File to load the StateRegistry from..
        forceLoad : bool, optional
            Load files that were written by a propably incompatible version of this package. The default is False.

        Raises
        ------
        TypeError
            Path was passed the wrong type.
        FileNotFoundError
            Given file does not exist.
        AssertionError
            You're trying to load a file, that was written by a propably incompatible version of this package. If you want to take your chances, use forceLoad=True.

        Returns
        -------
        register : StateRegistry
            Instance of this class.

        """
        
        # Was the path passed as string? If so convert it to a Path object.
        if isinstance(path, str): path = Path(path)
        # Is the path now the correct type?
        if not isinstance(path, Path):
            log.error(f"Unexpected type of path parameter. Expected string or pathlib.Path not {type(path)}!")
            raise TypeError(f"Unexpected type of path parameter. Expected string or pathlib.Path not {type(path)}!")
        log.info(f"Reading molecules from {path.resolve()} ...")
        
        # Does the file exist?
        if not path.exists():
            log.critical(f"{path.resolve()} doesn't exist.")
            raise FileNotFoundError(f"Can't read {path.resolve()}! File does not exist!")
        
        # Create an instance of StateRegistry.
        register = cls()
        
        # Read the meta data from the zip archive. To make sure the rest of the archive is compatible with the current version.
        log.info("Read meta data ...")
        with ZipFile(path.resolve(), "r") as archive:
            # Read the meta.yaml in the zip archive
            with archive.open(SAVE_TO_FILENAME["meta"]) as file:
                register._metadata = yaml.sage_load(file)
        
        # Check the meta data to make sure the rest of the archive is compatible with the current version.
        # This step can be skipped by using forceLoad=True
        if forceLoad == True:
            log.info("Ignoring meta info to load file. This may result in loading files that are incompatible with the current version of the package. Use forceLoad=False to safe load the file.")
        else:
            log.info("Checking versions of current package and version of package that saved the file.")
            # Get version object, that make the version easier to remember.
            currentVersion = version.parse(__version__)
            fileVersion = version.parse(register._metadata["version"])
            # Check the versions of current version and version that saved the file.
            if currentVersion.major != fileVersion:
                # Do the major version not match? -> Raise exception.
                log.critical(f"Major versions do not match! Refuse to load {path.resolve()}. Overrule this behaviour with forceLoad=True.")
                raise AssertionError(f"Major versions do not match! Refuse to load {path.resolve()}. Overrule this behaviour with forceLoad=True.")
            elif currentVersion < fileVersion:
                # Is the current version smaller than the file version? -> Raise exception.
                log.critical(f"Current version ({__version__}) is smaller than version which saved the file ({register._metadata['version']}). Refuse to load {path.resolve()}! Override this behaviour with forceLoad=True.")
                raise AssertionError(f"Major versions do not match! Refuse to load {path.resolve()}. Overrule this behaviour with forceLoad=True.")
        # Continue reading zip archive
        
        # Read the zip archive.
        with ZipFile(path.resolve(), "r") as archive:
            # Read the csv files and load them into the attributes of the StateRegistry.
            with archive.open(SAVE_TO_FILENAME["atom"]) as file:
                register._atomTable = pd.read_csv(file)
            with archive.open(SAVE_TO_FILENAME["molecule"]) as file:
                register._moleculeTable = pd.read_csv(file)
            with archive.open(SAVE_TO_FILENAME["bond"]) as file:
                register._bondTable = pd.read_csv(file)
            
            # Read the settings from yaml files and save them into StateRegistry.
            with archive.open(SAVE_TO_FILENAME["setting"]) as file:
                register._settings = yaml.safe_load(file)
        
        # Save the file this StateRegistry was read from in the settings dictionary.
        # Can be used to override the file with changes when user is done with the file.
        register._settings["fileLocation"] = path.resolve()
        
        # Return the instance of the StateRegistry.
        return register