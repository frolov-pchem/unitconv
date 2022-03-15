#!/usr/bin/env python
# Andrey I. Frolov, Apr 2013, ISC RAS, Ivanovo, Russia
#
#   Copyright 2013 Andrey I. Frolov
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
import sys,math
from copy import copy
from numpy import isnan
from FactorToConvertUnits import FactorToConvertUnits
from FactorToConvertUnits import ShowHelp


# Info variables
init_date='13.04.2013'; init_author='Andrey Frolov'; init_place='Ivanovo, ISC RAS'
lastmod_date='25.07.2013'; version='1.5'; lastmod_author='Andrey Frolov'; lastmod_place='Massy, Paris'

### Version 1.5 release notes:
# 1) Changed special keys from upper case to lowercase: e.g. toGMX to togmx
# 2) Extend an error message to display the second possible reason of unit sequence missmatch error (e.g. mol/badunit^3 nm^-3):
# {
# unitconv.py mol/badunit^3 nm^-3
# !Error in FactorForUnits: the unit type sequences you specified are not equal (for shorter sequence: ['Number', 'none', 'Distance'] and for longer sequence: ['Number', 'none']). The missing unit is not of type number. Returning NaN.
# !Due to an error /home/frolov/Mysoft/unitconv/unitconv.py returned NaN.
# }
# Extended the error message: sys.stderr.write('!Error in FactorForUnits: the unit type sequences you specified WERE NOT equal AND still not equal after trying to deduce missing unit (we have: for shorter sequence: '+str(utls)+' and for longer sequence: '+str(utll)+'). This might be caused together or independently by the following reasons: 1) the missing unit is not of type "number", 2) some of the units are not known, but this could not be caugth automatically. Please, consult help for known units. Returning NaN.\n')


### Version 1.4 release notes:
# 1) Extended special keys TODO: write appropriate help for this
# 2) Added kBT energy unit
# 3) Bug fix: conversion to/from units "j" was wrong. Factor was 1e-3 now is 1000. 
# 4) "cm" and "dm" added


### Version 1.3 release notes:
# 1) The code splited into two parts: file FactorToConvertUnits.py - module which can be loaded from any python code, and the unitconv.py - program which does the conversion of units in command line using methods (fanction) build in FactorToConvertUnits.py. To use the module one should set the unitconv path to PYTHONPATH

### Version 1.1 and 1.2 release notes:
# 1) Bug fix: occasional errors unit1 is not known
# 2) Added special unitkeys feature. See help for more details.
# 3) Added special untkey: WHAM
# 4) Few other minor changes

#-------------------------------------------------------------------------------
# Program start

# Show help is required
try: lh=sys.argv[1]
except: lh=''

if lh=='-h' or lh=='--help':
    ShowHelp()
    sys.exit(0)

# Get the units for conversion
try:
    units1=sys.argv[1]
    units2=sys.argv[2]
except:
    sys.stderr.write("      EXAMPLE USAGE: "+sys.argv[0]+"  kJ/mol/rad^2  kcal/mol/rad^2 \n")
    sys.stderr.write("      EXAMPLE USAGE: "+sys.argv[0]+"  -h \n")
    sys.stderr.write("      EXAMPLE USAGE: "+sys.argv[0]+" kcal/mol/rad^2  toGMX \n")
    sys.stderr.write('      !Error in '+sys.argv[0]+': you must provide 2 arguments OR -h or --help. Returning NaN.\n')
    print 'NaN'
    sys.exit(1)


F = FactorToConvertUnits(units1,units2)
if isnan(F) != True:
    print F
else:
    sys.stderr.write('!Due to an error '+sys.argv[0]+' returned NaN.\n')
    print 'NaN'

