from scipy.io.idl import readsav
import os
from bokeh.plotting import figure
from bokeh.io import vform,output_file, show, vplot
import numpy as np
from elements import ELEMENTS as ele
try:
    from sympy.mpmath import *
except:
    from mpmath import *

'''
molDict = {
    "mols": {
        "H2O":1e-3, 
        "H2": 99.0, 
        "CH4": 1e-3, 
        "CO":1e-3, 
        "C2H2":1e-3, 
        "HCN": 1e-3}, 
    "T":1500.0, 
    "g":, 
    "Rp":, 
    "R*":, 
    "beta":, 
    "P": 1e-4, 
    "P0": 10.0, 
    "taueq": 0.56
}
'''
k = 1.380658e-16 #cgs 
rsun = 6.96e10 #cm 
barye = 1000000.0 #1 bar to 1e6 barye
amu =1.66e-24 #grams
color = ["red", "blue", "green", "purple", "black", "yellow", "orange", "pink","cyan","brown"]
__XSEC__ = os.path.join(os.path.dirname(__file__), "reference/CrossSections/")

def readXSecs(molDict):
    
    
    MOL = molDict["mols"].keys()
    

    files = os.listdir(__XSEC__)
    output_file("cross.html")
    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"

        
    xsec = {}
    waves = {}

    #compute concentration of background gas 
    molDict = addBackground(molDict)
    
    #create dictionary of wavelengths and cross sections
    for nmol in range(0,len(MOL)):
        if MOL[nmol] == "H2": 
            for i in files: 
                if i.find("H2H2") != -1:
                    readthis = i  
                    
        elif MOL[nmol] == "He": 
            for i in files: 
                if i.find("H2He") != -1:
                    readthis = i
                    
        else:
            for i in files: 
                if i.find(MOL[nmol]+'_') != -1:
                    readthis = i
                    
         
        cross = readsav(os.path.join(__XSEC__,readthis))
        Tind = find_nearest(cross.T, molDict["T"])
        Pind = find_nearest(cross.P, molDict["P"])

        actualT = cross.T[Tind]
        actualP = cross.P[Pind]
        
        xsec[MOL[nmol]] = cross.xsecarr[Pind, Tind]
        waves[MOL[nmol]] = cross.wno
        
        
    return {"waves":waves, "xsec":xsec}
    
def computeAlpha(molDict, plot = False):
    plot1 = False
    plot2 = False 
    
    MOL = molDict["mols"].keys()
    read = readXSecs(molDict)
    waves = read["waves"]
    xsec = read["xsec"]
    
    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"
    if plot:
        plot2 = figure(plot_width=800, plot_height=200,  tools=TOOLS, responsive =False, y_axis_type="log",
                                 x_axis_label='Wavelength [um]', x_axis_type="log",
                                 y_axis_label='Weighted Cross Section', y_range = [1e-29, 1e-17])

    #compute mean molecular weight and beta 
    mu = computeMu(molDict["mols"])*amu
    
    beta = molDict["P0"]*barye / molDict["taueq"] * np.sqrt(2.0*np.pi*molDict["Rp"]*rsun)

    H = k*molDict["T"]/mu/molDict["g"]
    
    sumprod = computeZSum(xsec, molDict["mols"], waves, plot2,plot)
    
    Z = H * np.log( 1.0 / np.sqrt(k*molDict["T"]*mu*molDict["g"]) * beta
                    * sumprod)
    
    C = 2.0*(molDict["Rp"]*rsun)/(molDict["R*"]*rsun)**2.0
    
    alpha = C*Z
    
    if plot:
        plot1 = figure(plot_width=800, plot_height=200,  tools=TOOLS, responsive =False, #y_axis_type="log",
                                 x_axis_label='Wavelength [um]', x_axis_type="log",
                                 y_axis_label='Alpha Lambda', y_range = [min(alpha[alpha> 0.0]), max(alpha[alpha> 0.0])]) 
        plot1.line(1e4/waves[MOL[0]][alpha> 0.0], alpha[alpha> 0.0], alpha = 0.5, line_width = 3)
    
   #"plot1":plot1,"plot2":plot2, 
    return {"w":1e4/waves[MOL[0]][alpha> 0.0],
             "Z":Z, "alpha":alpha[alpha> 0.0], "H":H, "mu":mu, "xsec":xsec,"mols":molDict["mols"]}
    
def find_nearest(array,value):
    #small program to find the a temp and pressure match in the cross section data 
    idx = (np.abs(array-value)).argmin()
    return idx

def computeZSum(xsec, squig, waves, plot2,plot):
    sumprod = 0.0
    num = -1

    for i in squig.keys():  
        if i=="H2":
            sumprod += squig[i]*xsec[i]*squig["H2"]
            num +=1
            if plot:
                plot2.line(1e4/waves[i],squig[i]*xsec[i]*squig["H2"], color = color[num], legend = i)

        elif i=="He":
            sumprod += squig[i]*xsec[i]*squig["H2"]
            num +=1
            if plot:
                plot2.line(1e4/waves[i],squig[i]*xsec[i]*squig["H2"], color = color[num], legend = i)

        else:
            sumprod += squig[i]*xsec[i]
            num +=1
            if plot:
                plot2.line(1e4/waves[i],squig[i]*xsec[i], color = color[num], legend = i)

            
    return sumprod

def computeMu(squig): 
    """
    Compute mean molecular weight of an atmosphere 

    Input: 
        Dict of molecules with mixing ratios 
        i.e. {"H2O": 0.98, "H2":0.2}
    Output: 
        mean molecular weight, float

    """
    
    mu = 0 
    for i in squig.keys():
        totmass = 0.0
        for j in range(0, len(i)):
            try: 
                check = float(i[j])        
            except: 
                if i[j].isupper(): elem = ele[i[j]]
                try: 
                    num = float(i[j+1])
                except: 
                    if i[j].islower():
                        elem = ele[i[j-1:j+1]]
                        num = 1.0
                    else: 
                        num = 1.0
                totmass += elem.mass * num
            
        mu += squig[i] * totmass

    return mu
    
            
def addBackground(molDict):
    #compute concentration of background gas 
    if (molDict["mols"]["H2"] == "bkg") & (molDict["mols"]["He"] == "bkg"): 
        molDict["mols"]["H2"] = 0.0
        molDict["mols"]["He"] = 0.0
        fH2He = 1.0 - sum(molDict["mols"].values())
        try: 
            check1 = np.sqrt(fH2He)
            check2 = 1.0/fH2He
        except: 
            raise Exception("Mixing Ratios Exceed total sum = 1, no background gas")
        frac = molDict["fracH2He"]
        molDict["mols"]["H2"] = fH2He/(1.+frac)
        molDict["mols"]["He"] = frac*molDict["mols"]["H2"]
    elif (molDict["mols"]["H2"] == "bkg") & (molDict["mols"]["He"] != "bkg"): 
        molDict["mols"]["H2"] = 0.0
        molDict["mols"]["H2"] = 1.0 - sum(molDict["mols"].values())
    elif (molDict["mols"]["H2"] != "bkg") & (molDict["mols"]["He"] == "bkg"): 
        molDict["mols"]["He"] = 0.0
        molDict["mols"]["He"] = 1.0 - sum(molDict["mols"].values())
    else: 
        for name in molDict["mols"].keys():
            if molDict["mols"][name] == "bkg": 
                molDict["mols"][name] = 0.0
                molDict["mols"][name] = 1.0 - sum(molDict["mols"].values())
        

    return molDict   


####### Should probably eventually separate into two classes starting here. #######

    
def computeJac(molDict): 
    """
    Function to construct Jacobian of equation alpha http://arxiv.org/pdf/1511.09443v2.pdf
    eqn. 4 
    
    Inputs: 
        T, mu, g, beta: float
        xsec: dict of cross sections
        squig: dict of mixing ratios
    """
    squig = molDict["mols"]
    T = molDict["T"]
    g = molDict["g"]
    
    read = readXSecs(molDict)
    xsec = read["xsec"]
    
    #xsec= {
    #    "H2O":[1.0,2.0], 
    #    "H2": [1.0,2.0], 
    #    "He": [1.0,2.0],
    #    "CH4": [1.0,2.0], 
    #    "CO":[1.0,2.0], 
    #    "C2H2":[1.0,2.0], 
    #    "HCN": [1.0,2.0]}
    
    #compute mean molecular weight and beta 
    mu = computeMu(squig)*amu
    #mu = 1.0
    beta = molDict["P0"]*barye / molDict["taueq"] * np.sqrt(2.0*np.pi*molDict["Rp"]*rsun)

    
    flag = [0]*int(len(squig.keys())+3.0)
    lenW = len(xsec["CO"])
    J = np.zeros((lenW,len(flag)))
    mp.dps = 15; mp.pretty = True
    for i in range(0,lenW):
        for j in range(0,len(flag)):
            flag[j] = 1
            aCO = xsec["CO"][i]
            aH2H2= xsec["H2"][i]
            aH2O= xsec["H2O"][i]
            aCH4= xsec["CH4"][i]
            aHCN= xsec["HCN"][i]
            aC2H2= xsec["C2H2"][i]
            aH2He= xsec["He"][i]
            
            val = diff(lambda Td,mud,gd, xCO, xH2,xH2O,xCH4,xHCN,xC2H2,xHe: k*Td/mud/gd * log(1.0 / sqrt(k*Td*mud*gd) * beta * 
                                                (xCO*aCO+xH2*xH2*aH2H2+xH2O*aH2O+xCH4*aCH4+
                                                  xHCN*aHCN+xC2H2*aC2H2+xH2*xHe*aH2He)), 
                 (T,mu,g,squig["CO"],squig["H2"],squig["H2O"],squig["CH4"],squig["HCN"],squig["C2H2"],squig["He"]),
                             (flag[0],flag[1],flag[2],flag[3],flag[4],flag[5],flag[6],flag[7],flag[8],flag[9]))
            
            #since this does not stable enough to handle mu derivative .. calculate analytically  instead
            if j == 1: 
                val = -k*T/mu**2.0/g*np.log(beta/np.sqrt(k*T*mu*g)*(squig["CO"]*aCO+squig["H2"]*squig["H2"]
                                                                    *aH2H2+squig["H2O"]*aH2O+squig["CH4"]*aCH4+
                                                  squig["HCN"]*aHCN+squig["C2H2"]*aC2H2+squig["H2"]*squig["He"]*aH2He)) - (
                    k*T/2.0/mu**2.0/g)
            
            flag[j] = 0
            
            #print(i, j, val, type(val))
            J[i,j] = float(val)

    return J
