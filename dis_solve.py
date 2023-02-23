import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# ONLY LOOKING AT BINARY MIXTURES

data = []
nTray = 1

# Operation settings
xF= 0.5
xD= 0.9
xB= 0.1
L= 450 #[mol/s] # feed stream liquid phase
F = 500 #[mol/s]
optReflux = 1.2
D = ((xF-xB)/(xD-xB))*F
B = ((xD-xF)/(xD-xB))*F


# Chemicals used
components = {
        1: # 1 [Benzene]
        {"A": 9.2082,
        "B": 2755.64,
        "C": -54.00,
        "Psat": 0},
        2: # 2 [Toluene]
        {"A": 9.3716,
        "B": 3090.78,
        "C": -53.97,
        "Psat": 0}
}

def SP_OpeTemp(start,end):
    Ptot = 0 #[bar]
    # SP as function of operation temperature
    for Tope in range(start,end):
        # Tope given # T[K]

        Psat1=math.exp(components[1]["A"]-components[1]["B"]/(Tope+components[1]["C"]))
        Psat2=math.exp(components[2]["A"]-components[2]["B"]/(Tope+components[2]["C"]))
        pTOP = xD*Psat1+(1-xD)*Psat2
        pBOTTOM = xB*Psat1+(1-xB)*Psat2
        pAVG = (pTOP+pBOTTOM)/2
        # Seperation factor
            #binary
        SP=Psat1/Psat2

        # Parse values DF
        data.append(dict(zip(["Tope","SP","pAVG"],[Tope,SP,pAVG])))

    df = pd.DataFrame(data)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
    df.plot(x='Tope', y='SP', kind='line', ax=axes[0], xlabel="Operation Temperature [K]", ylabel="Seperation Factor", legend=None)
    df.plot(x='Tope', y='pAVG', kind='line', ax=axes[1], xlabel="Operation Temperature [K]", ylabel="Average Pressure [bar]", legend=None)
    fig.tight_layout()
    plt.show(block=True)
    return

def vapLiqSP(constSP, importdata):

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
        plt.plot(liquidX, vapourY, '-b')
    else:
    # import data
        dataEq = pd.read_csv("input.csv")
        df = pd.DataFrame(dataEq)
        df.plot(x='x', y='y', kind='line', legend=None)

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
    start = xD/(constSP*nTray-xD*(constSP-1))
    plt.hlines(y=xD, xmin=start, xmax=xD, color='k')
    YMAX = xD
    counter_var = 0
    for i in range(0,700): # number of plates in rectifying is end+1
        if start<XintersectF_D: # swap over to the stripping section
            if counter_var ==0:
                print(f"Feed tray location: {i+1}")
                counter_var += 1
            YMIN = YMAX
            YMAX = slopeStrip*start-intersectStrip*xB
            plt.vlines(x=start,ymax=YMAX,ymin=YMIN, color='k')
            XMAX = start
            start = YMAX/(constSP*nTray-YMAX*(constSP-1))
            if start<xB:
                break
            else:
                plt.hlines(y=YMAX, xmin=start, xmax=XMAX, color='k')
        else:
            YMIN = YMAX
            YMAX = (1/(Ropt+1))*xD+(Ropt/(Ropt+1))*start
            plt.vlines(x=start,ymax=YMAX,ymin=YMIN, color='k')
            XMAX = start
            start = YMAX/(constSP*nTray-YMAX*(constSP-1))
            plt.hlines(y=YMAX, xmin=start, xmax=XMAX, color='k')
    if start>xB:
        YMIN = YMAX
        YMAX = (1/(Ropt+1))*xD+(Ropt/(Ropt+1))*start
        plt.vlines(x=start,ymax=YMAX,ymin=YMIN, color='k')
    if(i== (0 | 699)):
        print("Error with plotting graph, check specific options, ntray, xF, xB, xD and constSP")
        return
    else:
        print(f"Minimum number of plates: {round(Nmin,2)}" +"\n" + f"Optimal Reflux Ratio: {round(Ropt,2)}" +"\n" + f"Number of theoretical plates: {i+1}")
        plt.show(block=True)
        return

#SP_OpeTemp(300,800)

vapLiqSP(2, False)


