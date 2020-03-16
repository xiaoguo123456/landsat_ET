from math import *
import numpy as np
import cupy as cp

def EA(NDVI,NDWI,Ta,RH,Rn,WS,DT):
    """
    :param LEcor:
    :param Ta:温度
    :param RH:相对湿度
    :param Rn: 净辐射
    :param WS:风速
    :param NDVI:Landsat-8 得到的NDVI
    :param DT:Tmax - Tmin 温差
    :return:
    """
    NDVI=cp.array(NDVI)
    NDWI=cp.array(NDWI)
    Ta=cp.array(Ta)
    RH=cp.array(RH)
    Rn=cp.array(Rn)
    WS=cp.array(WS)
    DT=cp.array(DT)
    sigma = 5.6697 * cp.power(10, -8)
    Cp = 1005.0
    gamma = 67.00
    Rd = 287.04
    kpar = 0.500
    z0 = 1.0000
    # RS-PM algorithm add()，subtract()，multiply() 和 divide()
    fc=cp.where(NDVI>=0.05,cp.divide(cp.subtract(NDVI,0.05),0.90),0)

    fc=cp.where(fc<=0,0,fc)
    fc=cp.where(fc>=1,1,fc)

    
    lai=cp.where((0<fc)&(fc < 1) ,cp.divide((-cp.log(1 - fc)),kpar),0)
    lai=cp.where(fc<=0,0.01,lai)
    lai=cp.where(fc>=1,6,lai)


    cl = 0.003  # unit: m s-1   from Yuan 2010
    Topt = 25.0
    VPDc = 2900
    VPDo = 650

    pa =cp.divide(101300.0,(cp.power(10.0,(cp.divide(z0,(cp.multiply(18400.0,(cp.add(1.0 , cp.divide(Ta,273.0))))))))))

    esat = cp.multiply(611.0 , cp.exp(cp.multiply(17.502 , cp.divide(Ta , (cp.add(Ta , 240.97))))))
    e0 = cp.add(esat , RH)
    mT = cp.exp(-cp.power((cp.divide(cp.subtract(Ta , Topt) , Topt)), 2))
    delta = 17.502 * 240.97 * esat / (Ta + 240.97) **2
    q_d = 0.622 * e0 / pa
    Tv = (Ta + 273.15) * (1 + 0.608 * q_d)
    rou = pa / (Rd * Tv)

    VPD = esat - e0
    mVPD=cp.divide(cp.subtract(VPDc,VPD),cp.subtract(VPDc,VPDo))

    gc = cl * mT * mVPD * lai#canopy conductance
    rc=cp.where(gc==0,0,1 / gc )   #canopy resistance

    rtot = 125.0  #from Yuan 2010
    rcorr=cp.where((pa==0)|(Ta == -273.15) ,0,1.0 / (((Ta + 273.15) / 293.15) ** 1.75 * (101300 / pa)))


    rtotc = rtot * rcorr
    rr=cp.where(Ta==-273.15,0,(rou * Cp) / (4.0 * sigma * (Ta + 273.15) **3.0))

    ra=cp.where(rtotc+rr==0,0,rtotc * rr / (rtotc + rr))

    Rnv = Rn * fc
    Rns = Rn * (1 - fc)

    LEcrop1=cp.where( ((delta+gamma*(1+rc/ra))==0.000)|(ra == 0.000),0,(delta * Rnv + rou * Cp * (esat - e0) / ra) / (delta + gamma * (1 + rc / ra)))

    evap_pot1=cp.where(((delta+gamma*rtotc/ra)==0.000)|(ra == 0.000),0,(delta * Rns + rou * Cp * (esat - e0) / ra) / (delta + gamma * rtotc / ra))

    soilc=cp.where(VPD<0.000,1,(RH) ** (VPD / 200.0))

    LEsoil1 = evap_pot1 * soilc
    LE_R = LEcrop1 + LEsoil1  #RS-PM algorithm
    LE_R=cp.where( (LE_R>=500.000)|(LE_R <= 0.000),0,LE_R)


    #SW algorithm
    ra1=cp.where(WS==0,0,(1 / (0.41 * 0.41 * WS)) * log((10.000 - 0.67 * 2.000) / (2.000 - 0.67 * 2.000)) * log((10.000 - 0.67 * 2.000) / (0.123 * 2.000)))

    dela1=cp.where(lai==-0.5,1,1 - (0.5 / (0.5 + lai)) * np.exp(-lai * lai / 8.000))

    uc = 0.83 * dela1 * WS + (1 - dela1) * WS
    raa=cp.where(WS==0,0,ra1 * (WS - uc) / WS)

    rac=cp.where(dela1*WS==0.000,0,uc * ra1 / (dela1 * WS))

    ras=cp.where((1-dela1)*WS==0,0,uc * ra1 / ((1 - dela1) * WS))


    G1 = 0.18 * (1 - fc) * Rn
    Rsc1 = rc
    Rss1 = cp.exp(8.206 - 4.255 * soilc)
    PMc=cp.where(delta*(raa+rac)+gamma*((raa+rac)+Rsc1)==0,0,(delta * (Rn - G1) * (raa + rac) + rou * Cp * (esat - e0) - delta * rac * (Rns - G1)) / (delta * (raa + rac) + gamma * ((raa + rac) + Rsc1)))
    PMs=cp.where(delta*(raa+ras)+gamma*((raa+ras)+Rss1)==0,0,(delta * (Rn - G1) * (raa + ras) + rou * Cp * (esat - e0) - delta * ras * Rnv) / (delta * (raa + ras) + gamma * ((raa + ras) + Rss1)))

    RRa1 = (delta + gamma) * raa
    RRs1 = (delta + gamma) * ras + gamma * Rss1
    RRc1 = (delta + gamma) * rac + gamma * Rsc1
    CC1=cp.where(RRs1*(RRc1+RRa1)+RRc1*RRa1==0,0,(RRs1 * (RRc1 + RRa1)) / (RRs1 * (RRc1 + RRa1) + RRc1 * RRa1))
    CS1=cp.where(RRc1*(RRs1+RRa1)+RRs1*RRa1==0,0,(RRc1 * (RRs1 + RRa1)) / (RRc1 * (RRs1 + RRa1) + RRs1 * RRa1))


    ET_SW = CC1 * PMc + CS1 * PMs
    ET_SW=cp.where((ET_SW<=0.000)|(ET_SW >= 500.000) ,0,ET_SW)

    #PT-JPL
    ag0 = 0.18
    lai = cp.where(lai == 0, 0.01, lai)
    fc4 = 1 - cp.exp(-0.5 * lai)
    G00 = Rn * (1 - fc4) * ag0
    Rns4 = Rn * cp.exp(0.9 * cp.log(1 - fc4))
    Rnv4 = Rn - Rns4

    m = 1.00
    b = -0.05
    beta = 1000
    Topt = 25
    PTc = 1.26
    fpar = cp.multiply(1.24 ,NDVI) - 0.168
    fparmax = 1.0000

    ft = cp.exp(-((Ta - Topt) / Topt) ** 2)
    fm =  fpar / fparmax
    fsm = RH ** (VPD / beta)
    fipar = m * NDVI + b
    fwet4 = RH ** 4

    fg = cp.where(fipar < 0, 0, fpar / fipar)
    fg = cp.where(fg > 1, 1, fg)
    LEsoil4 = cp.where(delta == -gamma, 0, (fwet4 + fsm * (1 - fwet4)) * PTc * delta * (Rns4 - G00) / (delta + gamma))
    LEcrop4 = cp.where(delta == -gamma, 0, (1 - fwet4) * fg * ft * fm * PTc * delta * Rnv4 / (delta + gamma))
    LEi4 = cp.where(delta == -gamma, 0, fwet4 * PTc * delta * Rnv4 / (delta + gamma))
    LE_JPL = LEsoil4 + LEcrop4 + LEi4
    LE_JPL = cp.where((LE_JPL <= 0.000) | (LE_JPL >= 500.000), 0, LE_JPL)


    # MS-PT
    DTmax = 40
    fsm = cp.where(DT<1.000,1,0)
    fsm = cp.where((DT>=1.000)&(DT/DTmax>=0.000), (1 / DT) ** (DT / DTmax), fsm)

    fsm = cp.where(DT/DTmax<0.000, 1, fsm)

    fwet = fsm**4
    ag = 0.18
    alpha = 1.26

    fcc = (NDVI - 0.05) / 0.900
    Rnvv = Rn * fcc
    Rnss = Rn * (1 - fcc)
    G0 = ag * (1 - fcc) * Rn
    LEs=cp.where((delta+gamma)==0.000,0,(1 - fwet) * alpha * fsm * (delta / (delta + gamma)) * (Rnss - G0))
    LEws=cp.where((delta+gamma)==0.000,0,fwet * alpha * (delta / (delta + gamma)) * (Rnss - G0))
    LEc=cp.where((delta+gamma)==0.000,0,(1 - fwet) * alpha * fcc * mT * (delta / (delta + gamma)) * Rnvv)
    LEic=cp.where((delta+gamma)==0.000,0,fwet * alpha * (delta / (delta + gamma)) * Rnvv)

    LE_MS = LEs + LEws + LEc + LEic
    LE_MS=cp.where((LE_MS<=0.000)|(LE_MS >= 500.000) ,0,LE_MS)



    # Wang2008
    ET_wang = Rn * (0.1440 + 0.6495 * NDVI + 0.009 * Ta - 0.0163 * DT)
    ET_wang=cp.where((ET_wang<=0.000)|(ET_wang >= 500.000) ,0,ET_wang)


    #water ET
    water_ET=1.26*delta/(delta+gamma)*0.74*Rn



    ET_products=0.201*LE_R+0.198*ET_SW+0.191*LE_JPL+0.205*LE_MS+0.205*ET_wang
    ET_products=cp.where(NDWI>0.05,water_ET,ET_products)
    ET_products = cp.where(ET_products<0,0, ET_products)


    return ET_products





