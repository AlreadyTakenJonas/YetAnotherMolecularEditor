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
    "element",      # How many protons does the atom has?
    "mass",         # How much mass does the atom has?
    "charge",       # How much charge does the atom has? 
    "cartesian_x",  # Where is the atom in 3D space?
    "cartesian_y",  # Where is the atom in 3D space?
    "cartesian_z",  # Where is the atom in 3D space?
    "fixedInternalCoordinates", # Is the atom fixed in space (internal coordinates)?
    "color",        # What color does the atom have when shown on the screen?
    "radius"        # How big is the atom?
]
# What will be saved about every molecule?
MOLECULE_PROPERTY_LIST = [
    "moleculeID",   # Unique ID
    "name"          # Descriptive name
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
                log.error(f"Can't save molecules to {path.resolve}, because it already exists.")
                raise FileExistsError(f"Can't save molecules to {path.resolve()}, because it already exists.")
        
        # Is the suffix right? If not it is alright.
        if path.suffix != ".jam": log.warning(f"Unexpected file suffix '{path.suffix}'. Why not save it as .jam? Saving it anyways.")
        
        # Write the state of the StateRegistry to a temporary directory and move them into a zip archive.
        with tempfile.TemporaryDirectory() as tmp:
            # Write data on atoms, molecules and bonds to csv files
            self._atomTable.to_csv( tmp+"/atoms.csv" )
            self._moleculeTable.to_csv( tmp+"/molecules.csv" )        
            self._bondTable.to_csv( tmp+"/bonds.csv" )
            # Dump the settings dictionary to yaml file
            dumpedSettings = yaml.safe_dump(self._settings)
            Path(tmp+"/settings.yaml").write_text(dumpedSettings)
            # Dump the meta data dictionary to yaml file
            dumpedMetadata = yaml.sage_dump(self._metadata)
            Path(tmp+"/meta.yaml").write_text(dumpedMetadata)
            
            # Move all files in temporary directory to zip archive
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
            with archive.open("meta.yaml") as file:
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
            with archive.open("atoms.csv") as file:
                register._atomTable = pd.read_csv(file)
            with archive.open("molecules.csv") as file:
                register._moleculeTable = pd.read_csv(file)
            with archive.open("bonds.csv") as file:
                register._bondTable = pd.read_csv(file)
            
            # Read the settings from yaml files and save them into StateRegistry.
            with archive.open("settings.yaml") as file:
                register._settings = yaml.safe_load(file)
        
        # Save the file this StateRegistry was read from in the settings dictionary.
        # Can be used to override the file with changes when user is done with the file.
        register._settings["fileLocation"] = path.resolve()
        
        # Return the instance of the StateRegistry.
        return register