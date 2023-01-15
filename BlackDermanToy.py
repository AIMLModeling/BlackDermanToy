from __future__ import division
import math
from scipy.optimize import fsolve

r = [0.1,0.11,0.12,0.125,0.13]
VolativityYield = [0.20,0.19,0.18,0.17,0.16]
precision = 0.001
InitValue=0.18
PV1 = 100/(1+r[0])**1  
PV2 = 100/(1+r[1])**2  
PV3 = 100/(1+r[2])**3
PV4 = 100/(1+r[3])**4
PV5= 100/(1+r[4])**5

def print_lattice(lattice, info = []):
    print ("")
    print ("")
    levels = len(lattice[-1])
    start_date = len(lattice[0]) - 1
    dates = levels - start_date 
    outlist = []
    col_widths = [0] * dates
    for j in range(levels):
        level = []
        for k in range(dates):
            try:
                point = "{:.2f}".format(lattice[k][levels - 1 - j]*100)
                esc_width = 0 
                if info != [] and info[k][levels - 1 - j] > 0:
                    point = (point, 'red')
                    esc_width += 9 
                level.append(point)
                col_widths[k] = max(col_widths[k], len(point) - esc_width)
            except IndexError:
                level.append('')
        outlist.append(level)
    separator = "|-".join(['--' * w for w in col_widths])
    formats = [ ]
    for k in range(dates):
        formats.append("     %%%ds" % col_widths[k])
    pattern = "  ".join(formats)
    print (pattern % tuple(str(start_date + time) for time in range(dates)))
    print (separator)
    for line in outlist:
        print (pattern % tuple(line))
    print ("")
    
def TreeOneStep(x):
    rd = x[0]
    ru = x[1]
    Nu = (100)/(1+ru)
    Nd = (100)/(1+rd)
    out = [(0.5*((Nu/(1+r[0])) + (Nd/(1+r[0])))-PV2)]
    out.append((0.5 * math.log(ru/rd) - VolativityYield[1] ))
    return out
    
resultTreeOneStep = fsolve(TreeOneStep,[0.1,0.1],xtol=precision)
ru = resultTreeOneStep[1]
rd = resultTreeOneStep[0]

def TreeTwoSteps(x):
    ruu = x[1]*x[1]/x[0]
    rud = x[1]
    rdd = x[0]
    N1 = (100)/(1+ruu)
    N2 = (100)/(1+rud)
    N3 = (100)/(1+rdd)
    ans1 = (0.5*N1 + 0.5*N2)/(1+ru)
    ans2 = (0.5*N2 + 0.5*N3)/(1+rd)
    yu = math.pow(100/ans1, 1/2) - 1
    yd = math.pow(100/ans2, 1/2) - 1
    out = [((0.5*ans1 + 0.5*ans2)/(1+r[0])-PV3)]
    out.append((0.5 * math.log(yu/yd) - VolativityYield[2]))
    return out

resultTreeTwoSteps = fsolve(TreeTwoSteps,[InitValue,InitValue],xtol=precision)
ruu = resultTreeTwoSteps[1]**2/resultTreeTwoSteps[0]
rud = resultTreeTwoSteps[1]
rdd = resultTreeTwoSteps[0]

def TreeThreeSteps(x):
    ruud = x[1]*x[1]/x[0]
    rdud = x[1]
    rddd = x[0]
    ruuu = ruud*ruud/rdud
  
    N1 = 100/(1+ruuu) 
    N2 = 100/(1+ruud)
    N3 = 100/(1+rdud)
    N4 = 100/(1+rddd)
    
    ans1 = ((0.5*N1 + 0.5*N2))/(1 + ruu)
    ans2 = ((0.5*N2 + 0.5*N3))/(1 + rud)
    ans3 = ((0.5*N3 + 0.5*N4))/(1 + rdd)
    
    Bu = ((0.5*ans1 + 0.5*ans2))/(1 + ru)
    Bd = ((0.5*ans2 + 0.5*ans3))/(1 + rd)
    
    yu = math.pow(100/Bu, 1/3)-1
    yd = math.pow(100/Bd, 1/3)-1
    out = [((0.5*Bu + 0.5*Bd)/(1+r[0]) - PV4)]
    out.append((0.5 * math.log((yu)/(yd)) - VolativityYield[3] ))
    return out

resultTreeThreeSteps = fsolve(TreeThreeSteps,[InitValue,InitValue],xtol=precision )

ruud = resultTreeThreeSteps[1] * resultTreeThreeSteps[1]/resultTreeThreeSteps[0]
rdud = resultTreeThreeSteps[1]
rddd = resultTreeThreeSteps[0]
ruuu = ruud*ruud/rdud

def TreeFourSteps(x):
    rudud = x[1]  * x[1]/x[0]
    ruddd = x[1]
    rdddd = x[0]  
    ruuud = rudud*rudud/ruddd
    ruuuu = ruuud*ruuud/rudud
  
    N1 = 100/(1+ruuuu) 
    N2 = 100/(1+ruuud)
    N3 = 100/(1+rudud)
    N4 = 100/(1+ruddd)
    N5 = 100/(1+rdddd)
    
    ans1 = ((0.5*N1 + 0.5*N2))/(1 + ruuu)
    ans2 = ((0.5*N2 + 0.5*N3))/(1 + ruud)
    ans3 = ((0.5*N3 + 0.5*N4))/(1 + rdud)
    ans4 = ((0.5*N4 + 0.5*N5))/(1 + rddd)
    
    fans1 = ((0.5*ans1 + 0.5*ans2))/(1 + ruu)
    fans2 = ((0.5*ans2 + 0.5*ans3))/(1 + rud)
    fans3 = ((0.5*ans3 + 0.5*ans4))/(1 + rdd)
    
    Bu = ((0.5*fans1 + 0.5*fans2))/(1 + ru)
    Bd = ((0.5*fans2 + 0.5*fans3))/(1 + rd)
    
    yu = math.pow(100/Bu, 1/4)-1
    yd = math.pow(100/Bd, 1/4)-1
    out = [((0.5*Bu + 0.5*Bd)/(1+r[0]) -PV5 )]
    out.append(((0.5 * math.log(yu/yd) - VolativityYield[4])))
    return out

resultTreeFourSteps = fsolve(TreeFourSteps,[InitValue,InitValue],xtol=precision)

rudud = resultTreeFourSteps[1] * resultTreeFourSteps[1]/resultTreeFourSteps[0]
ruddd = resultTreeFourSteps[1]
rdddd = resultTreeFourSteps[0] 
ruuud = rudud*rudud/ruddd
ruuuu = ruuud*ruuud/rudud

finalRate = [[0.1],[ru,rd],[ruu,rud,rdd],[ruuu,ruud,rdud,rddd],[ruuuu,ruuud,rudud,ruddd,rdddd]]
finalRate2 = [list(reversed(x)) for x in finalRate]
print_lattice(finalRate2, info = [])
