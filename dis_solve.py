import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# ONLY BINARY MIXTURES

data = []
naturalAntione = False
plotGraph = True
finalPrint = False
print(f"Plotting graph is set to {str(plotGraph)}")

nTray = 0.8
Emv = 0.8
Ht = 0.1 #[m]

HvapSteam = 2100 #[kJ/kg]

# Operation settings
xF= 0.03
xD= 0.395
xB= 0.01
L= 99.045 #[mol/s] # feed stream liquid phase
F = 93 #[mol/s]
feedunit = "mol/s"
optReflux = 1.2


# Chemicals used
components = {
        1: # 1 [ethanol]
        {
        "name": "ethanol",
        "A": 5.24677,
        "B": 1598.673,
        "C": -46.424,
        "Vap_enthalpy": 40.2, #[kJ/mol]
        "Mol_weight": 46.07*10**-3, # [kg/mol]
        "Psat": 0},
        2: # 2 [cyclohexene]
        {
        "name": "cyclohexene",
        "A": 4.13983,
        "B": 1316.554,
        "C": -35.581,
        "Vap_enthalpy": 31.1, #[kJ/mol]
        "Mol_weight": 84.16*10**-3, # [kg/mol]
        "Psat": 0}
}

def SP_OpeTemp(start,end):
    # SP as function of operation temperature
    for Tope in range(start,end):
        # Tope given # T[K]

        if naturalAntione == True:
            Psat1=math.exp(components[1]["A"]-components[1]["B"]/(Tope+components[1]["C"]))
            Psat2=math.exp(components[2]["A"]-components[2]["B"]/(Tope+components[2]["C"]))
        else:
            Psat1=10**(components[1]["A"]-components[1]["B"]/(Tope+components[1]["C"]))
            Psat2=10**(components[2]["A"]-components[2]["B"]/(Tope+components[2]["C"]))
        pTOP = xD*Psat1+(1-xD)*Psat2
        pBOTTOM = xB*Psat1+(1-xB)*Psat2
        pAVG = (pTOP+pBOTTOM)/2
        # Seperation factor
            #binary
        SP=Psat1/Psat2

        # Parse values DF
        data.append(dict(zip(["Tope","SP","pAVG", "pTOP", "pBOTTOM", components[1]["name"], components[2]["name"]],[Tope,SP,pAVG, pTOP, pBOTTOM,Psat1,Psat2])))

    df = pd.DataFrame(data)

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(3, 5), constrained_layout = True)
    df.plot(x='Tope', y='SP', kind='line', ax=axes[0][0], xlabel="Operation Temperature [K]", ylabel="Seperation Factor", legend=None)
    df.plot(x='Tope', y='pAVG', kind='line', ax=axes[0][1], xlabel="Operation Temperature [K]", ylabel="Average Pressure [bar]", legend=None)
    df.plot(x='Tope', y='pTOP', kind='line', ax=axes[0][2], xlabel="Operation Temperature [K]", ylabel="Pressure at distillate [bar]", legend=None)
    df.plot(x='Tope', y='pBOTTOM', kind='line', ax=axes[1][0], xlabel="Operation Temperature [K]", ylabel="Pressure at bottom [bar]", legend=None)
    df.plot(x='Tope', y=components[1]["name"], kind='line', ax=axes[1][1], xlabel="Operation Temperature [K]")
    df.plot(x='Tope', y=components[2]["name"], kind='line', ax=axes[1][1], xlabel="Operation Temperature [K]", ylabel="Vapour pressure [bar]")
    fig.delaxes(axes[1][2])
    if plotGraph == True:
        plt.show(block=True)
    return

def vapLiqSP(constSP, importdata, xB, xD, xF, F, L, columnsFinalDF=False):

    D = ((xF-xB)/(xD-xB))*F
    B = ((xD-xF)/(xD-xB))*F

    

    # Minimum amount of trays:
    Nmin = (math.log((xD/(1-xD))/(xB/(1-xB))))/(math.log(constSP))

    # Minimum reflux ratio:
    Rmin = ((xD/xF)-constSP*((1-xD)/(1-xF)))/(constSP-1)
    
    # Optimum reflux ratio determined on costs:
    Ropt = optReflux*Rmin

    # Vapour fraction as a function of liquid fraction using seperation factor
    if importdata == False:
        intersectEqDiagonal = (constSP*nTray-1)/(constSP-1)
        if(xD>intersectEqDiagonal):
            print("xD is out of range, change xD")
            return
        liquidX = np.linspace(0,intersectEqDiagonal,100)
        vapourY = (constSP*liquidX)/(1+(constSP-1)*(liquidX))*nTray
        plt.plot(liquidX, vapourY, '--b')
        if nTray != 1:
            vapourY2 = (importdata[0] + importdata[1]*liquidX + importdata[2]*liquidX**2 + importdata[3]*liquidX**3)
            plt.plot(liquidX, vapourY2, '-b')
    else:
    # import data
        liquidX = np.linspace(0,1,100)
        vapourY = (importdata[0] + importdata[1]*liquidX + importdata[2]*liquidX**2 + importdata[3]*liquidX**3)*nTray
        plt.plot(liquidX, vapourY, '--b')
        if nTray != 1:
            vapourY2 = (importdata[0] + importdata[1]*liquidX + importdata[2]*liquidX**2 + importdata[3]*liquidX**3)
            plt.plot(liquidX, vapourY2, '-b')

    # Determination of q:
    q=L/F

    #Find cross point feed and rectifying
    if q==1: # not the clean way
        XintersectF_D = (xF/(1-0.999999999)-xD/(Ropt+1))/(0.999999999/(1-0.999999999)+Ropt/(Ropt+1))
    else:
        XintersectF_D = (xF/(1-q)-xD/(Ropt+1))/(q/(1-q)+Ropt/(Ropt+1))
    YitersectF_D = (1/(Ropt+1))*xD+(Ropt/(Ropt+1))*XintersectF_D

    # Rectifying section operation line:
    xRectfyingSpace = np.linspace(XintersectF_D,xD,100)
    yRect = (1/(Ropt+1))*xD+(Ropt/(Ropt+1))*xRectfyingSpace
    plt.plot(xRectfyingSpace, yRect, '-g')

    # Feed line:
    if q==1:
        plt.vlines(x = xF, ymin = xF, ymax = YitersectF_D, colors = 'purple')
    else:
        xFeedSpace = np.linspace(XintersectF_D,xF,100)
        yFeed = xF/(1-q)+(-q/(1-q))*xFeedSpace
        plt.plot(xFeedSpace, yFeed, '-m')
    # calculate Rmin from intersecting the feedline and the equillibrium curve
    if importdata != False:
        unfiltered = np.roots([importdata[3], importdata[2], importdata[1]+(q/(1-q)), importdata[0]-(xF/(1-q))])
        filtered = [i for i in unfiltered if i.imag == 0]
        filtered = [i for i in filtered if i >=0]
        filtered = [i for i in filtered if i <=1]
        intersectFeedEquillibriumX = filtered[0]
        intersectFeedEquillibriumY = xF/(1-q)+(-q/(1-q))*intersectFeedEquillibriumX
        RminNew = (intersectFeedEquillibriumY-xD)/(-1*(intersectFeedEquillibriumY+intersectFeedEquillibriumX))
        print(Rmin)
        print(RminNew)
    # Stripping section operation line:
        #Find slope and intersect
    slopeStrip = (YitersectF_D-xB)/(XintersectF_D-xB)
    intersectStrip = slopeStrip-1
        #Plotting
    xStrippingSpace = np.linspace(xB,XintersectF_D,100)
    yStrip = slopeStrip*xStrippingSpace-intersectStrip*xB
    plt.plot(xStrippingSpace, yStrip, '-g')

    # Straight line:
    x = np.linspace(0,1,100)
    plt.plot(x, x, '-r')
    plt.xlabel("Liquid fraction")
    plt.ylabel("Vapour fraction")
    # Number of stages lines:
    if importdata == False:
        start = xD/(constSP*nTray-xD*(constSP-1))
    else: 
        unfiltered = np.roots([importdata[3], importdata[2], importdata[1],importdata[0]-xD/nTray])
        filtered = [i for i in unfiltered if i.imag == 0]
        start = filtered[0]
    plt.hlines(y=xD, xmin=start, xmax=xD, color='k')
    YMAX = xD
    counter_var = 0
    for i in range(0,700): # number of plates in rectifying is end+1
        if start<XintersectF_D: # swap over to the stripping section
            if counter_var ==0:
                FeedTray = i+1 
                counter_var += 1
            YMIN = YMAX
            YMAX = slopeStrip*start-intersectStrip*xB
            plt.vlines(x=start,ymax=YMAX,ymin=YMIN, color='k')
            XMAX = start
            if importdata == False:
                start = YMAX/(constSP*nTray-YMAX*(constSP-1))
            else:
                unfiltered = np.roots([importdata[3], importdata[2], importdata[1],importdata[0]-YMAX/nTray])
                filtered = [i for i in unfiltered if i.imag == 0]
                start = filtered[0]
            if start<xB:
                break
            else:
                plt.hlines(y=YMAX, xmin=start, xmax=XMAX, color='k')
        else:
            YMIN = YMAX
            YMAX = (1/(Ropt+1))*xD+(Ropt/(Ropt+1))*start
            plt.vlines(x=start,ymax=YMAX,ymin=YMIN, color='k')
            XMAX = start
            if importdata == False:
                start = YMAX/(constSP*nTray-YMAX*(constSP-1))
            else:
                unfiltered = np.roots([importdata[3], importdata[2], importdata[1],importdata[0]-YMAX/nTray])
                filtered = [i for i in unfiltered if i.imag == 0]
                start = filtered[0]
            plt.hlines(y=YMAX, xmin=start, xmax=XMAX, color='k')
    if start>xB:
        YMIN = YMAX
        YMAX = (1/(Ropt+1))*xD+(Ropt/(Ropt+1))*start
        plt.vlines(x=start,ymax=YMAX,ymin=YMIN, color='k')
    if(i== (0 | 699)):
        print("Error with plotting graph, check specific options, ntray, xF, xB, xD and constSP")
        return
    else:
        Nt = (i+2)/Emv-1
        V_prim_prim = B/intersectStrip
        EnergyReboiler = xB*V_prim_prim*components[1]["Vap_enthalpy"]+(1-xB)*V_prim_prim*components[2]["Vap_enthalpy"]
        if finalPrint == True:
            print(f"Separation factor: {constSP}" + "\n" + f"Stage efficiency: {Emv}" + "\n" + f"q: {round(q,3)}" + "\n" + f"Stage height: {Ht} m" + "\n" + f"Minimum number of plates: {round(Nmin,2)}" + "\n" + f"Minimum reflux ratio: {round(Rmin,2)}" +"\n" + f"Optimal Reflux Ratio: {round(Ropt,2)}" +"\n" + f"Feed tray location: {FeedTray}" + "\n" + f"Number of theoretical plates: {i+2}" + "\n" + f"Actual number of stages: {round(Nt, 2)}" + "\n" + f"Column height: {round(Nt * Ht, 2)} m" + "\n" + f"Energy consumption reboiler: {round(EnergyReboiler, 2)} kJ/s" + "\n" + f"Steam consumption: {round(EnergyReboiler/HvapSteam, 2)} kg/s" "\n" + f"F: {round(F,2)} {feedunit}" + "\n" + f"B: {round(B,2)} {feedunit}" + "\n" + f"D: {round(D, 2)} {feedunit}" + "\n" + f"L': {round(optReflux*D,2)} {feedunit}" + "\n" + f"V': {round((optReflux+1)*D,2)} {feedunit}" + "\n" + f"V'': {round(V_prim_prim, 2)} {feedunit}" + "\n" + f"L'': {round(B/intersectStrip+B, 2)} {feedunit}")
        
        d = {'SeparationFactor': [constSP], 'xF': [xF], 'xB': [xB],'xD': [xD], 'q':[q], 'StageEfficiency': [Emv], 'StageHeight': [Ht], 'TrayEffieciency' :[nTray], "F": [F], 'Nmin' : [round(Nmin,2)], 'Rmin' : [round(Rmin,2)], 'Ropt':[round(Ropt,2)], 'NTheoretical': [i+2], 'NActual' : [round(Nt, 2)], "FeedTray": [FeedTray], "ColumnHeight": [round(Nt * Ht, 2)], "ReboilerEnergyConsumption" : [round(EnergyReboiler, 2)], "SteamConsumption" : [round(EnergyReboiler/HvapSteam, 2)], "B": [round(B,2)],"D": [round(D, 2)], "L'": [round(optReflux*D,2)], "V'": [round((optReflux+1)*D,2)], "V''": [round(V_prim_prim, 2)], "L''": [round(B/intersectStrip+B, 2)]}
        if columnsFinalDF:
            c = {}
            for i in columnsFinalDF:
                c[i] = [d[i][0]]
        try:
            df = pd.read_csv("raw_data.csv")
            if columnsFinalDF:
                df1= pd.DataFrame(c)
            else:
                df1= pd.DataFrame(d)
            df2 = pd.concat([df, df1], ignore_index=False)
            df2.to_csv('raw_data.csv', index=False)
        except FileNotFoundError:
            if columnsFinalDF:
                df = pd.DataFrame(data=c, index=[0])
            else:
                df = pd.DataFrame(data=d, index=[0])
            df.to_csv('raw_data.csv', index=False)
        if plotGraph == True:
            plt.show(block=True)
        else:
            plt.close('all') 
        return
def SP_X_Y():
    df = pd.read_csv("data_project\\output.csv")
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5, 5))
    df.plot.scatter(x='x1', y='y1', color="green", marker="v", ax=axes[0])
    df.plot.scatter(x='x1', y='y2', ax=axes[0], ylabel="y1")
    axes[0].legend(["y1", "y2"])
    df.plot.scatter(x='x1', y='a', ylabel='Seperation Factor', ax=axes[1])
    plt.show()
    #DATA range plotter
# for i in tqdm(np.arange(0.001,0.03,0.001)):
#     vapLiqSP(5.34, [0.07951,2.91250,-5.70358,3.49134], xB=i, xF=xF, xD=xD, F=F, L=L, columnsFinalDF=["xB", "D", "ReboilerEnergyConsumption"])
# fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(3, 5))
# df = pd.read_csv("raw_data.csv")
# df.plot(x='xB', y='diff', kind='scatter', legend=None)
# df.to_csv('raw_data.csv', index=False)
# plt.show()
# df.plot(x='xB', y='D', kind='scatter', ylabel="D (mol/s)", legend=None)
# plt.show()
# df = pd.read_csv("data_project\\raw_data_changeXb.csv")
# 
# plt.show()
# for i in tqdm(np.arange(0,93,0.5)):
#     vapLiqSP(1.09, False, xB=xB, xF=xF, xD=xD, F=F, L=i)
#SP_OpeTemp(323,367)

vapLiqSP(5.34, [0.07951,2.91250,-5.70358,3.49134], xB=xB, xF=xF, xD=xD, F=F, L=L)


