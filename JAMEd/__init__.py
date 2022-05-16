#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Export modules of this package
from .version import __version__
from . import *

# Setup a unit registry for the whole package
# Needed to work with units
from pint import UnitRegistry
PINT_UNIT_REGISTRY = UnitRegistry()
#QUANT = PINT_UNIT_REGISTRY.Quantity
