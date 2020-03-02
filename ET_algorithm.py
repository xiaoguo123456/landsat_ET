from math import *
import numpy as np


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
    sigma = 5.6697 * pow(10, -8)
    Cp = 1005.0
    gamma = 67.00
    Rd = 287.04
    kpar = 0.500
    z0 = 1.0000
    # RS-PM algorithm
    fc=np.where(NDVI>=0.05,(NDVI - 0.05) / 0.90,0)

    fc=np.where(fc<=0,0,fc)
    fc=np.where(fc>=1,1,fc)

    np.seterr(divide='ignore', invalid='ignore')
    lai=np.where((0<fc)&(fc < 1) ,(-np.log(1 - fc)) / kpar,0)
    lai=np.where(fc<=0,0.01,lai)
    lai=np.where(fc>=1,6,lai)


    cl = 0.003  # unit: m s-1   from Yuan 2010
    Topt = 25.0
    VPDc = 2900
    VPDo = 650

    pa = 101300.0 / (10.0 ** (z0 / (18400.0 * (1.0 + Ta / 273.0))))

    esat = 611.0 * np.exp(17.502 * Ta / (Ta + 240.97))
    e0 = esat * RH
    mT = np.exp(-((Ta - Topt) / Topt)** 2)
    delta = 17.502 * 240.97 * esat / (Ta + 240.97) **2
    q_d = 0.622 * e0 / pa
    Tv = (Ta + 273.15) * (1 + 0.608 * q_d)
    rou = pa / (Rd * Tv)

    VPD = esat - e0
    mVPD=np.where(VPDc-VPDo==0,0,(VPDc - VPD) / (VPDc - VPDo))

    gc = cl * mT * mVPD * lai#canopy conductance
    rc=np.where(gc==0,0,1 / gc )   #canopy resistance

    rtot = 125.0  #from Yuan 2010
    rcorr=np.where((pa==0)|(Ta == -273.15) ,0,1.0 / (((Ta + 273.15) / 293.15) ** 1.75 * (101300 / pa)))


    rtotc = rtot * rcorr
    rr=np.ma.where(Ta==-273.15,0,(rou * Cp) / (4.0 * sigma * (Ta + 273.15) **3.0))

    ra=np.ma.where(rtotc+rr==0,0,rtotc * rr / (rtotc + rr))

    Rnv = Rn * fc
    Rns = Rn * (1 - fc)

    LEcrop1=np.ma.where( ((delta+gamma*(1+rc/ra))==0.000)|(ra == 0.000),0,(delta * Rnv + rou * Cp * (esat - e0) / ra) / (delta + gamma * (1 + rc / ra)))

    evap_pot1=np.where(((delta+gamma*rtotc/ra)==0.000)|(ra == 0.000),0,(delta * Rns + rou * Cp * (esat - e0) / ra) / (delta + gamma * rtotc / ra))

    soilc=np.ma.where(VPD<0.000,1,(RH) ** (VPD / 200.0))

    LEsoil1 = evap_pot1 * soilc
    LE_R = LEcrop1 + LEsoil1  #RS-PM algorithm
    LE_R=np.ma.where( (LE_R>=500.000)|(LE_R <= 0.000),0,LE_R)


    #SW algorithm
    ra1=np.where(WS==0,0,(1 / (0.41 * 0.41 * WS)) * log((10.000 - 0.67 * 2.000) / (2.000 - 0.67 * 2.000)) * log((10.000 - 0.67 * 2.000) / (0.123 * 2.000)))

    dela1=np.where(lai==-0.5,1,1 - (0.5 / (0.5 + lai)) * np.exp(-lai * lai / 8.000))

    uc = 0.83 * dela1 * WS + (1 - dela1) * WS
    raa=np.where(WS==0,0,ra1 * (WS - uc) / WS)

    rac=np.where(dela1*WS==0.000,0,uc * ra1 / (dela1 * WS))

    ras=np.where((1-dela1)*WS==0,0,uc * ra1 / ((1 - dela1) * WS))


    G1 = 0.18 * (1 - fc) * Rn
    Rsc1 = rc
    Rss1 = np.exp(8.206 - 4.255 * soilc)
    PMc=np.where(delta*(raa+rac)+gamma*((raa+rac)+Rsc1)==0,0,(delta * (Rn - G1) * (raa + rac) + rou * Cp * (esat - e0) - delta * rac * (Rns - G1)) / (delta * (raa + rac) + gamma * ((raa + rac) + Rsc1)))
    PMs=np.where(delta*(raa+ras)+gamma*((raa+ras)+Rss1)==0,0,(delta * (Rn - G1) * (raa + ras) + rou * Cp * (esat - e0) - delta * ras * Rnv) / (delta * (raa + ras) + gamma * ((raa + ras) + Rss1)))

    RRa1 = (delta + gamma) * raa
    RRs1 = (delta + gamma) * ras + gamma * Rss1
    RRc1 = (delta + gamma) * rac + gamma * Rsc1
    CC1=np.where(RRs1*(RRc1+RRa1)+RRc1*RRa1==0,0,(RRs1 * (RRc1 + RRa1)) / (RRs1 * (RRc1 + RRa1) + RRc1 * RRa1))
    CS1=np.where(RRc1*(RRs1+RRa1)+RRs1*RRa1==0,0,(RRc1 * (RRs1 + RRa1)) / (RRc1 * (RRs1 + RRa1) + RRs1 * RRa1))


    ET_SW = CC1 * PMc + CS1 * PMs
    ET_SW=np.where((ET_SW<=0.000)|(ET_SW >= 500.000) ,0,ET_SW)

    #PT-JPL
    ag0 = 0.18
    lai = np.where(lai == 0, 0.01, lai)
    fc4 = 1 - np.exp(-0.5 * lai)
    G00 = Rn * (1 - fc4) * ag0
    Rns4 = Rn * np.exp(0.9 * np.log(1 - fc4))
    Rnv4 = Rn - Rns4

    m = 1.00
    b = -0.05
    beta = 1000
    Topt = 25
    PTc = 1.26
    fpar = 1.24 * NDVI - 0.168
    fparmax = 1.0000

    ft = np.exp(-((Ta - Topt) / Topt) ** 2)
    fm = np.where(fparmax == 0, 1, fpar / fparmax)
    fsm = RH ** (VPD / beta)
    fipar = m * NDVI + b
    fwet4 = RH ** 4

    fg = np.where(fipar < 0, 0, fpar / fipar)
    fg = np.where(fg > 1, 1, fg)
    LEsoil4 = np.where(delta == -gamma, 0, (fwet4 + fsm * (1 - fwet4)) * PTc * delta * (Rns4 - G00) / (delta + gamma))
    LEcrop4 = np.where(delta == -gamma, 0, (1 - fwet4) * fg * ft * fm * PTc * delta * Rnv4 / (delta + gamma))
    LEi4 = np.where(delta == -gamma, 0, fwet4 * PTc * delta * Rnv4 / (delta + gamma))
    LE_JPL = LEsoil4 + LEcrop4 + LEi4
    LE_JPL = np.where((LE_JPL <= 0.000) | (LE_JPL >= 500.000), 0, LE_JPL)


    # MS-PT
    DTmax = 40
    fsm = np.where(DT<1.000,1,0)
    fsm = np.where((DT>=1.000)&(DT/DTmax>=0.000), (1 / DT) ** (DT / DTmax), fsm)

    fsm = np.where(DT/DTmax<0.000, 1, fsm)

    fwet = fsm**4
    ag = 0.18
    alpha = 1.26

    fcc = (NDVI - 0.05) / 0.900
    Rnvv = Rn * fcc
    Rnss = Rn * (1 - fcc)
    G0 = ag * (1 - fcc) * Rn
    LEs=np.where((delta+gamma)==0.000,0,(1 - fwet) * alpha * fsm * (delta / (delta + gamma)) * (Rnss - G0))
    LEws=np.where((delta+gamma)==0.000,0,fwet * alpha * (delta / (delta + gamma)) * (Rnss - G0))
    LEc=np.where((delta+gamma)==0.000,0,(1 - fwet) * alpha * fcc * mT * (delta / (delta + gamma)) * Rnvv)
    LEic=np.where((delta+gamma)==0.000,0,fwet * alpha * (delta / (delta + gamma)) * Rnvv)

    LE_MS = LEs + LEws + LEc + LEic
    LE_MS=np.where((LE_MS<=0.000)|(LE_MS >= 500.000) ,0,LE_MS)



    # Wang2008
    ET_wang = Rn * (0.1440 + 0.6495 * NDVI + 0.009 * Ta - 0.0163 * DT)
    ET_wang=np.where((ET_wang<=0.000)|(ET_wang >= 500.000) ,0,ET_wang)


    #water ET
    water_ET=1.26*delta/(delta+gamma)*0.74*Rn



    ET_products=0.201*LE_R+0.198*ET_SW+0.191*LE_JPL+0.205*LE_MS+0.205*ET_wang
    ET_products=np.where(NDWI>0.05,water_ET,ET_products)
    ET_products = np.where(ET_products<0,0, ET_products)


    return ET_products





