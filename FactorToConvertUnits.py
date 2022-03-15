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

# Info variables
init_date='13.04.2013'; init_author='Andrey Frolov'; init_place='Ivanovo, ISC RAS'
lastmod_date='25.07.2014'; version='1.5'; lastmod_author='Andrey Frolov'; lastmod_place='Massy, Paris'


# Universal constants
Na=6.02214129e23
pi=math.pi
kB=1.3806488e-23
T298=298.15

### Conversion factors
# Internal units:
#     Energy = 'kJ'
#     Number = 'mol'
#     Distance = 'nm'
#     Angle = 'degree'

# The list of conversion factors to internal units
# Extend this list if units/unit types you need are misiing
UnitFact = {'Energy'  : {'kj':1.0,'kcal':0.239005736, 'kbt':1.0/(kB*T298/1000.0),'j':1000, 'kb': 1000/kB}, 
            'Distance': {'nm':1.0,'angstr':10.0,'m':1.0e-9,'cm':1.0e-7,'dm':1e-8},
            'Angle'   : {'degree':1.0,'rad':1*pi/180},
            'Number'  : {'mol':1.0,'number':Na}
           }


kcal_mol_angstr_rad =  { 'Energy'  : 'kcal', 
                         'Distance': 'angstr',
                         'Angle'   : 'rad',
                         'Number'  : 'mol'
                       }

kj_mol_nm_rad = { 'Energy'  : 'kj', 
                  'Distance': 'nm',
                  'Angle'   : 'rad',
                  'Number'  : 'mol'
                }

GMXunits = copy(kj_mol_nm_rad) 
WHAMunits = copy(kcal_mol_angstr_rad)
AMBERunits = copy(kcal_mol_angstr_rad)
CHARMMunits = copy(kcal_mol_angstr_rad)
RISMMOLunits = copy(kcal_mol_angstr_rad)

# Mind lower case letters for the dictionary keys 
SpecialKeyUnits = { 'togmx':   GMXunits,
                    'towham':  WHAMunits,
                    'toamber': AMBERunits,
                    'tocharmm': CHARMMunits,
                    'torismmol': RISMMOLunits,
                    'kj_mol_nm_rad': kj_mol_nm_rad,
                    'kcal_mol_angstr_rad': kcal_mol_angstr_rad
                  }


#-------------------------------------------------------------------------------
# Internal functions:

# The function returns the factor to conver units
def FactorToConvertUnits(UnitIn,UnitOut):
    # Vars in:
    #       UnitIn - string with the units to convert from
    #       UnitOut - string with the units to convert to
    # Returning: the conversion factor, or nan if error occured

    # The function parses a string with units 
    def SplitUnitString(unitstring):
        # Vars in:
        #       unitstring - string with units   (e.g. "kJ*mol^-1")     
        # Returning: 
        #   ulv - list with the unit sequence in the given string (e.g. ['kJ','mol']) 
        #   ulp - list with the unit power sequence (e.g. [1,-1])

        ustr=unitstring.replace(' ','')
        ustr=ustr.replace('_','*')

        ul_init=ustr        
        
        ul=[] # list for units
        pf=[1.0] # list for unit power factors (e.g. if there is / sign, then the power should be inversed)
        tmp=ustr
        for l in ustr:
            if l == '*':
                val,tmp=tmp.split('*',1)
                ul.append(val)
                pf.append(1.0)
            elif l == '/':
                val,tmp=tmp.split('/',1)
                ul.append(val)
                pf.append(-1.0)
            else:
                continue
        ul.append(tmp)            

        ulv=[]
        ulp=[]
        for i,val in enumerate(ul):
            if '^' in val:
                tmp=val.split('^');
                if len(tmp) != 2: 
                    sys.stderr.write('!Error in SplitUnitString: more than 2 power signs [^] in a field. Returning NaN.\n')
                    return float('NaN')
                unit=tmp[0]
                try: power=float(tmp[1])*pf[i]
                except: 
                    sys.stderr.write('!Error in SplitUnitString: cant evaluate power for unit ['+unit+']. Either power is NaN ['+tmp[1]+'] or power factor is NaN ['+str(pf[i])+']. Returning NaN.\n')
                    return float('NaN')
            else:
                unit=val
                try: power=1.0*pf[i]
                except: 
                    sys.stderr.write('!Error in SplitUnitString: cant evaluate power for unit ['+unit+']. Either power is NaN ['+tmp[1]+'] or power factor is NaN ['+str(pf[i])+']. Returning NaN.\n')
                    return float('NaN')
            ulv.append(unit) 
            ulp.append(power) 
        

        ult=[]
        for v in ulv: uts = get_unit_type(v); ult.append(uts)
        return (ulv,ulp,ult)

    # Function calculates the conversion factor for an element in the list [ulv] having in mind he power [ulp]
    def FactorForUnit(unit1,p1, unit2,p2):
        # Vars in:
        #       unit1 - string unit to convert from (e.g. "kJ")     
        #       p1 - power for the unit (e.g. "1")     
        #       unit2 - string unit to convert to (e.g. "kcal")     
        #       p2 - power for the unit (e.g. "1")     
        # Returning: 
        #   factor to convert unit1^p1 to unit2^p1  (p1 must be equal to p2) 
        f1=float('NaN')
        f2=float('NaN')
        for UF in UnitFact.itervalues():
            unit1=unit1.lower()
            unit2=unit2.lower()

            #print unit1,UF.keys()
            if unit1 in UF.keys(): 
                f1=UF[unit1]
                if unit2 in UF.keys(): 
                    f2=UF[unit2]
                else:
                    sys.stderr.write('!Error in FactorForUnit: unit2 ['+unit2+'] is not known OR it is not in the same group as unit1 ['+unit1+']. Returning NaN.\n')
                    return float('NaN')
                #else:
                #break
                #continue
        if isnan(f1):
                sys.stderr.write('!Error in FactorForUnit: unit1 ['+unit1+'] is not known. Returning NaN.\n')
                return float('NaN')

        if p1 != p2:
                sys.stderr.write('!Error in FactorForUnit: power1 ['+str(p1)+'] is not equal to power2 ['+str(p2)+']. Returning NaN.\n')
                return float('NaN')
        else:
            try: p=float(p1)
            except: 
                sys.stderr.write('!Error in FactorForUnit: power1 ['+str(p1)+'] is not a number. Returning NaN.\n')
                return float('NaN')

        return (f2/f1)**p
   
    # The function determines the unit type for given unit string 
    def get_unit_type(unit):
        # Vars in: unit - string unit (e.g. "kJ")     
        # Returning: unit type (e.g. "Energy") 
        unit=unit.lower()
        for UF,UFname in zip(UnitFact.itervalues(), UnitFact):
            #pnrint UF,UFname
            if unit in UF.keys(): 
                return UFname
        return 'none'
  
  
    # Parse the given unit strings
    [ulv1,ulp1,ult1]=SplitUnitString(UnitIn)
    [ulv2,ulp2,ult2]=SplitUnitString(UnitOut)


    # Special cases
    if ulv2[0].lower() in SpecialKeyUnits.keys(): 
        ult2=ult1
        ulp2=ulp1
        key=ulv2[0].lower()
        ulv2=[]
        for v in ult1: 
            if v!='none': ulv2.append(SpecialKeyUnits[key][v])
            else: ulv2.append('none')            
    
    

    # Special cases 
    #if ulv2[0].lower() == 'togmx':
    #    ult2=ult1
    #    ulp2=ulp1
    #    ulv2=[]
    #    for v in ult1: 
    #        if v!='none': ulv2.append(GMXunits[v])
    #        else: ulv2.append('none')            

    #if ulv2[0].lower() == 'towham':
    #    ult2=ult1
    #    ulp2=ulp1
    #    ulv2=[]
    #    for v in ult1: 
    #        if v!='none': ulv2.append(WHAMunits[v])
    #        else: ulv2.append('none')            
    #print ulv1,ulp1,ult1
    #print ulv2,ulp2,ult2


    # Take care of the cases when the Number units are not specified explicitly (e.g. one converts "mol*m^-3" "nm^-3") 
    if len(ulv1)!=len(ulv2):
        if len(ulv1)>len(ulv2): longer_arr=ulv1; shorter_arr=ulv2 ; lap=ulp1; sap=ulp2; utll=ult1; utls=ult2
        else: longer_arr=ulv2; shorter_arr=ulv1; lap=ulp2; sap=ulp1 ; utll=ult2; utls=ult1 
        
        #utls=[]; utll=[]
        #for v in shorter_arr: uts = get_unit_type(v); utls.append(uts)
        #for v in longer_arr: utl = get_unit_type(v); utll.append(utl)

        for i,v in enumerate(longer_arr):
            ut1=get_unit_type(v)
            if i<=(len(shorter_arr)-1): ut2=get_unit_type(shorter_arr[i])
            else: ut2='none'
            
            if ut1!=ut2:
                shorter_arr.insert(i,"number") ; sap.insert(i,lap[i]); utls.insert(i,utll[i])
            
        if utls != utll:

            sys.stderr.write('!Error in FactorForUnits: the unit type sequences you specified WERE NOT equal AND still not equal after trying to deduce missing unit (we have: for shorter sequence: '+str(utls)+' and for longer sequence: '+str(utll)+'). This might be caused together or independently by the following reasons: 1) the missing unit is not of type "number", 2) some of the units are not known, but this could not be caugth automatically. Please, consult help for known units. Returning NaN.\n')
            return float('NaN')

        elif sap != lap:
            sys.stderr.write('!Error in FactorForUnits: the unit sequences equal but the sequence of unit powers mismatch (for shorter sequence: '+str(sap)+' and for longer sequence: '+str(lap)+'). Returning NaN.\n')
            return float('NaN')
        else:
            pass

    # Calculate the conversion factor
    factor=1.0
    #print ulv1, ulp1, ulv2, ulp2
    for u1,p1,u2,p2 in zip(ulv1,ulp1,ulv2,ulp2):
        f=FactorForUnit(u1,p1,u2,p2)
        factor=factor*f

    return factor



#---------------------------------------------------------------------
# Function to show the help
def ShowHelp():
    print '------------------------------------------------------------------------------------------'
    print '     The unit converter. Written on '+init_date+' by '+init_author+' in '+init_place
    print '     The current version '+version+' is released on '+lastmod_date+' by '+lastmod_author+' in '+lastmod_place
    print 
    print "     EXAMPLE USAGE: "+sys.argv[0]+" kJ/mol/rad^2  kcal/mol/rad^2 "
    print "     EXAMPLE USAGE: "+sys.argv[0]+" mol*m^-3  nm^-3 "
    print "     EXAMPLE USAGE: "+sys.argv[0]+" -h "
    print "     EXAMPLE USAGE: "+sys.argv[0]+" kcal/mol  toGMX "
    print 
    print "     The script takes exactly two command line arguments. It returns the factor to convert units in the first argument to units in the second argument. The special signs [*] and [_] mean multiplication, [/] means division, [^] means power. Case insensitive. If the first argument is -h or --help then the help message will be printed and nothing returned. If the second argument is a special unitkey (e.g. [toGMX]) then the units to convert to are generated from a special prestored list, based on the unit type sequence given in the first argument (e.g. [unitconv.py  kcal/mol toGMX] corresponds to [unitconv.py kcal/mol kJ/mol], because [kJ] and [mol] are intenal units of Gromacs code (GMX)). If the scripts fails the float(\'NaN\') is retuned."
    print  
    print '     The following units supported so far (easy to extend in the source code):'
    for UF,UFname in zip(UnitFact.itervalues(),UnitFact):
        print '        '+UFname.ljust(8)+':',UF.keys()

    
    unitkeys=['togmx','towham']
    unitkeysinfo=['convert to Gromacs (www.gromacs.org) internal units',
                  'convert to WHAM code by Grossfield lab (http://membrane.urmc.rochester.edu/content/wham/) units' ]
    unitkeydict=['GMXunits','WHAMunits']

    print  
    print '     The following special unitkeys supported so far (easy to extend in the source code):'
    for uk,uki,ukd in zip(unitkeys,unitkeysinfo,unitkeydict):
        sys.stdout.write('        '+uk.ljust(8)+'- '+uki+':  '); exec 'print '+ukd

    print
    
    return 

