import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


def KDR_COEFF(seas, lang):
   
    pi = np.pi
    se_angle = np.arcsin(seas)

    if seas >= np.sin(lang):
        par = seas * np.cos(lang)
    else:
        par = (2/pi) * (seas * np.cos(lang) * 
                        np.arcsin(np.tan(se_angle) / np.tan(lang)) +
                        ((np.sin(lang))**2 - seas**2)**0.5)

    lec = par / seas

    return lec

def INTERNAL_CO2(tlp, dpr, fvdp, ca2, cc4):
   
    svl = 0.611 * np.exp(17.4 * tlp / (tlp + 239.))
    vdl = max(0, svl - dpr)
    
    # Constants
    km2 = 650 
    ko2 = 450 
    o2c = 210
    evcm = 65330
    ekmc = 79430
    ekom = 36380
    erad = 46390
    rdv2 = 0.0089
    
    # Temperature corrections
    kmc2 = km2 * np.exp((1./298. - 1./(tlp + 273.)) * ekmc / 8.314)
    kom2 = ko2 * np.exp((1./298. - 1./(tlp + 273.)) * ekom / 8.314)
    
    # Calculating GAMMAX
    gmx2 = 0.5 * np.exp(-3.3801 + 5220./(tlp + 273.) / 8.314) * o2c * kmc2 / kom2
    
    # Corrected dark respiration rate
    rdvc2 = rdv2 * np.exp((1/298 - 1/(tlp + 273)) * (erad - evcm) / 8.314)
    
    # Calculating GAMMA0 and GAMMA
    g02 = (gmx2 + rdvc2 * kmc2 * (1 + o2c / kom2)) / (1 - rdvc2)
    gam2 = g02 / 10
    
    # Internal CO2 concentration adjustment
    rcic2 = 1 - (1 - gam2 / ca2) * (0.14 + fvdp * vdl)
    co2in = rcic2 * ca2
    
    return svl, co2in

def LEAF_CONDUCTANCE(plf, rdlf, tlf, ca, ci, rbw, rt):

    gc = (plf - rdlf) * (273. + tlf) / 0.53717 / (ca - ci)
    rs = max(1E-30, 1 / gc - rbw * 1.3 - rt) / 1.6
    return rs

def TEMP_DIFF(nrc, pt, rbh, rtt):

    lhv = 2.4E6  # Latent heat of vaporization
    vhca = 1200  # Volumetric heat capacity of air

    tdif = (nrc - lhv * pt) * (rbh + rtt) / vhca
    tdif = min(max(tdif, -25), 25)  # Ensuring tdif is within [-25, 25] range

    return tdif

def SINK_GROWTH(dvr, te, tx, ti):

    fd = dvr * (2 * te - tx) * (te - ti) / te / (te - tx) ** 2 * (ti / te) ** (tx / (te - tx))
    return fd

def WEATHER(doy, lat, insp):
    
    rad = np.pi / 180 

    dec = np.arcsin(np.sin(23.45 * rad) * np.cos(2 * np.pi * (doy + 10) / 365))

    sld = np.sin(rad * lat) * np.sin(dec)  # Sin of latitude * declination
    cld = np.cos(rad * lat) * np.cos(dec)  # Cos of latitude * declination
    aob = sld / cld

    dl = 12.0 * (1. + 2. * np.arcsin(aob) / np.pi)
    dlp = 12.0 * (1. + 2. * np.arcsin((-np.sin(insp * rad) + sld) / cld) / np.pi)
    dsb = 3600. * (dl * sld + 24. * cld * np.sqrt(1. - aob * aob) / np.pi)
    dsbe = 3600. * (dl * (sld + 0.4 * (sld * sld + cld * cld * 0.5)) +
                       12.0 * cld * (2.0 + 3.0 * 0.4 * sld) * np.sqrt(1.0 - aob * aob) / np.pi)

    sc = 1367. * (1. + 0.033 * np.cos(2. * np.pi * (doy - 10) / 365))

    return sc, sld, cld, dl, dlp, dsbe


def PHOTOSYNTHESIS(tl, ci, par, km25, ko25, eard, eakm,rdvx, eako, eavc, eaj, xv, xj, sjt, dej, o2c, hhc, phm, tht):

    up = 4.56 * par

    kmc = km25 * math.exp((1./298. - 1./(tl + 273.)) * eakm / 8.314)
    kom = ko25 * math.exp((1./298. - 1./(tl + 273.)) * eako / 8.314)

    gmx = 0.5 * math.exp(-3.3801 + 5220. / (tl + 273.) / 8.314) * o2c * kmc / kom

    vct = math.exp((1./298. - 1./(tl + 273.)) * eavc / 8.314)
    jt = math.exp((1./298. - 1./(tl + 273.)) * eaj / 8.314) * \
        (1. + math.exp(sjt / 8.314 - dej / 298. / 8.314)) / \
        (1. + math.exp(sjt / 8.314 - 1. / (tl + 273.) * dej / 8.314))

    vmx = xv * vct
    jmx = xj * jt

    a2 = (1 - (1 - ((0 * hhc * (0 + 3 * ci + 7 * gmx) / (4 * ci + 8 * gmx) + 1)) / \
                (hhc * (0 + 3 * ci + 7 * gmx) / (4 * ci + 8 * gmx) - 1))) / \
             (1 + (1 - (1 - ((0 * hhc * (0 + 3 * ci + 7 * gmx) / (4 * ci + 8 * gmx) + 1)) / \
                         (hhc * (0 + 3 * ci + 7 * gmx) / (4 * ci + 8 * gmx) - 1))) / phm)
    xp = a2 * up / max(1E-10, jmx)
    j2 = jmx * (1 + xp - ((1 + xp)**2 - 4 * xp * tht)**0.5) / 2 / tht

    vcr = vmx * ci / (ci + kmc * (o2c / kom + 1.))
    vjr = j2 * ci / (1 + 3 * ci + 7 * gmx)  

    alr = min(vcr, vjr)
    psyn = max(1E-10, (1E-6) * 44 * alr)

    rdr = math.exp((1/298 - 1/(tl + 273)) * eard / 8.314) 
    rd = (1E-6) * 44 * rdvx * rdr  

    return psyn, rd

def PHENOLOGY(dys, sldp, dlpr, spst, epst, psn, tdv, tdr, tday):
    
    if sldp < 0:
        mopt = 18 
        dlpt = min(mopt, dlpr)
    else:
        mopt = 11 

    if dys < spst or dys > epst:
        efpt = 1  
    else:
        efpt = max(0, 1 - psn * (dlpt - mopt)) 

    if 0 <= dys < 1.0:
        dvrt = 1 / tdv * tday * efpt  
    else:
        dvrt = 1 / tdr * tday  

    return dvrt
    
def SOIL_EVAP(day_temp, d_vp, rss, rts, rbws, rbhs, atrjs, atmtr, pt1, wsup1):
   
    svp = 0.611 * math.exp(17.4 * day_temp / (day_temp + 239.))
    vpd = max(0., svp - d_vp)
    slope = 4158.6 * svp / (day_temp + 239.) ** 2
    
    fpe, fnrads = PENMAN_MONTEITH(rss, rts, rbws, rbhs, atrjs, atmtr, 1, day_temp, d_vp, slope, vpd)
    fpesol = max(0, fpe)
    faesol = min(fpesol, fpesol / (pt1 + fpesol) * wsup1)
    fdifs = TEMP_DIFF(fnrads, faesol, rbhs, rts)
    tavs = day_temp + fdifs

    svps = 0.611 * math.exp(17.4 * tavs / (tavs + 239.))
    slopes = (svps - svp) / AVOID_ZERO_DIVISION(fdifs)  
    pe, radi = PENMAN_MONTEITH(rss, rts, rbws, rbhs, atrjs, atmtr, 1, tavs, d_vp, slopes, vpd)
    soil_e = max(0, pe)

    return soil_e, radi

    
def ACTIVE_N(sn_tot, sn_min, leaf_area_index, k_n, k_beam):

    n_upper_canopy = sn_tot * (1. - np.exp(-k_n * leaf_area_index)) / k_n - sn_min * leaf_area_index

    n_sunlit = sn_tot * (1. -  np.exp(-(k_n + k_beam) * leaf_area_index)) / (k_n + k_beam) \
               - sn_min * (1. -  np.exp(-k_beam * leaf_area_index)) / k_beam
    n_shaded = n_upper_canopy - n_sunlit

    return n_sunlit, n_shaded

def LIGHT_ABSORB(soil_cover, k_b, k_b_prime, k_diffuse_prime, p_cb, p_cd, i_beam_0, i_diffuse_0, leaf_area_index):
    
    i_canopy = (1. - p_cb) * i_beam_0 * (1. - np.exp(-k_b_prime * leaf_area_index)) \
               + (1. - p_cd) * i_diffuse_0 * (1. - np.exp(-k_diffuse_prime * leaf_area_index))

    i_sunlit = (1 - soil_cover) * i_beam_0 * (1 - np.exp(-k_b * leaf_area_index)) \
               + (1 - p_cd) * i_diffuse_0 / (k_diffuse_prime + k_b) * k_diffuse_prime \
               * (1 - np.exp(-(k_diffuse_prime + k_b) * leaf_area_index)) \
               + i_beam_0 * ((1 - p_cb) / (k_b_prime + k_b) * k_b_prime \
                             * (1 - np.exp(-(k_b_prime + k_b) * leaf_area_index)) \
                             - (1 - soil_cover) * (1 - np.exp(-2 * k_b * leaf_area_index)) / 2)

    i_shaded = i_canopy - i_sunlit

    return i_sunlit, i_shaded

def KDF_COEFF(lai, leaf_angle, scp):

    k_b_15 = KDR_COEFF(np.sin(15. * pi / 180.), leaf_angle)
    k_b_45 = KDR_COEFF(np.sin(45. * pi / 180.), leaf_angle)
    k_b_75 = KDR_COEFF(np.sin(75. * pi / 180.), leaf_angle)

    k_diffuse_prime = -1 / lai * np.log(0.178 * np.exp(-k_b_15 * (1.0 - scp)**0.5 * lai) +
                                     0.514 * np.exp(-k_b_45 * (1.0 - scp)**0.5 * lai) +
                                     0.308 * np.exp(-k_b_75 * (1.0 - scp)**0.5 * lai))

    return k_diffuse_prime



def CANOPY_PHOTO_TRANSPIRATION(SF,sC, sldL, csdL, dLg, dSBe, dDTR, tmx, tmn, vDP, wnM, cTyp, lAI, tLAI, hT, lWd, rDp, sD1, rSS, bLA, kN, kW, sLN, sLT, sLNN, sLMN, wSUP, c2A, lStr, eJM, xVN, xJN, thT, wCL, fvPD):
    ig = 5
    xG = np.array([0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899])
    wG = np.array([0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635])
    ppC, apCS, apCN, apC, ptC, atC, peS, aeS, dfS, dfSU, dfSH, dpR = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    for i in range(ig):
        srT = 12 - 0.5 * dLg
        hr = srT + dLg * xG[i]
        sinB = max(0., sldL + csdL * np.cos(2. * np.pi * (hr - 12.) / 24.))
        dTR = dDTR * (sinB * sC / 1367) / dSBe
        dTp = tmn + (tmx - tmn) * np.sin(np.pi * (hr + dLg / 2 - 12) / (dLg + 3))
        wSP = wSUP * (sinB * sC / 1367) / dSBe
        wSP1 = wSP * sD1 / rDp
        wND = wnM  
        pAR = 0.5 * dTR
        nIR = 0.5 * dTR
        aTR = pAR / (0.5 * sC * sinB)
        frDF = 1 if aTR < 0.22 else 1 - 6.4 * (aTR - 0.22) ** 2 if 0.22 < aTR <= 0.35 else 1.47 - 1.66 * aTR
        frDF = max(frDF, 0.15 + 0.85 * (1 - np.exp(-0.1 / sinB))) if aTR < 0.22 else max(1 - 6.4 * (aTR - 0.22) ** 2, 0.15 + 0.85 * (1 - np.exp(-0.1 / sinB))) if 0.22 < aTR <= 0.35 else max(1.47 - 1.66 * aTR, 0.15 + 0.85 * (1 - np.exp(-0.1 / sinB)))
        pDF = pAR * frDF
        pDR = pAR - pDF
        nDF = nIR * frDF
        nDR = nIR - nDF
        bL = bLA * np.pi / 180  
        kB = KDR_COEFF(sinB, bL)
        scPR = 0.2 
        scNR = 0.8 
        kDPR = KDF_COEFF(tLAI, bL, scPR)
        kDNR = KDF_COEFF(tLAI, bL, scNR)
        kbPR, pcbPR = REFLECTION_COEFF(scPR, kB)
        kbNR, pcbNR = REFLECTION_COEFF(scNR, kB)
        pcdPR = 0.057 
        pcdNR = 0.389 
        rT = 0.74 * (np.log((2. - 0.7 * hT) / (0.1 * hT))) ** 2 / (0.4 ** 2 * wND)
        rTS = 0.74 * (np.log(56.)) ** 2 / (0.4 ** 2 * wND)
        frSU = 1. / kB / lAI * (1. - np.exp(-kB * lAI))
        frSH = 1. - frSU
        gbHLF = 0.01 * np.sqrt(wND / lWd)
        gbHC = (1 - np.exp(-0.5 * kW * lAI)) / (0.5 * kW) * gbHLF
        gbHSU = (1 - np.exp(-(0.5 * kW + kB) * lAI)) / (0.5 * kW + kB) * gbHLF
        gbHSH = gbHC - gbHSU
        rbHSU = 1. / gbHSU 
        rbWSU = 0.93 * rbHSU
        rbHSH = 1. / gbHSH  
        rbWSH = 0.93 * rbHSH 
        rbHS = 172. * np.sqrt(0.05 / max(0.1, wND * np.exp(-kW * tLAI)))
        rbWS = 0.93 * rbHS
        npSU, npSH = ACTIVE_N(sLT, sLMN, lAI, kN, kB)
        npSUN, npSNH = ACTIVE_N(sLNN, sLMN, lAI, kN, kB)
        apRSU, apRSH = LIGHT_ABSORB(scPR, kB, kbPR, kDPR, pcbPR, pcdPR, pDR, pDF, lAI)
        anRSU, anRSH = LIGHT_ABSORB(scNR, kB, kbNR, kDNR, pcbNR, pcdNR, nDR, nDF, lAI)
        apR = apRSU + apRSH
        atJSU = apRSU + anRSU
        atJSH = apRSH + anRSH
        psPAR = 0.1 
        psNIR = INSWTICH(wCL - 0.5, 0.52 - 0.68 * wCL, 0.18)
        atJS = (1 - psPAR) * (pDR * np.exp(-kbPR * tLAI) + pDF * np.exp(-kDPR * tLAI)) + \
               (1 - psNIR) * (nDR * np.exp(-kbNR * tLAI) + nDF * np.exp(-kDNR * tLAI))
        plFSU, ptSU, rWSU, nrSU, slSU = POTENTIAL_PHOTO_TRANSPIRATION(frSU, dTp, vDP, c2A, cTyp, fvPD, apRSU, npSU, rbWSU, rbHSU, rT * frSU, atJSU, aTR, eJM, xVN, xJN, thT)
        plFSH, ptSH, rWSH, nrSH, slSH = POTENTIAL_PHOTO_TRANSPIRATION(frSH, dTp, vDP, c2A, cTyp, fvPD, apRSH, npSH, rbWSH, rbHSH, rT * frSH, atJSH, aTR, eJM, xVN, xJN, thT)
        iPp = plFSU + plFSH
        iPpt = ptSU + ptSH
        pT1 = iPpt * sD1 / rDp
        iPe, nrDS = SOIL_EVAP(dTp, vDP, rSS, rTS, rbWS, rbHS, atJS, aTR, pT1, wSP1)
        iAe = min(iPe, iPe / (pT1 + iPe) * wSP1)
        iAt = min(iPpt, pT1 / (pT1 + iPe) * wSP1 + wSP - wSP1)
        aTSU = ptSU / iPpt * iAt
        aTSH = ptSH / iPpt * iAt
        aDfS = TEMP_DIFF(nrDS, iAe, rbHS, rTS)
        pASU, pANSU, aDfSU = SF*POTENTIAL_PHOTO_TRANSPIRATION(dTp, apRSU, vDP, c2A, cTyp, fvPD, nrSU, aTSU, ptSU, rT * frSU, rbHSU, rbWSU, rWSU, slSU, npSU, npSUN, eJM, xVN, xJN, thT)
        pASH, pANSH, aDfSH = SF*POTENTIAL_PHOTO_TRANSPIRATION(dTp, apRSH, vDP, c2A, cTyp, fvPD, nrSH, aTSH, ptSH, rT * frSH, rbHSH, rbWSH, rWSH, slSH, npSH, npSNH, eJM, xVN, xJN, thT)
        iAPS = pASU + pASH
        iAPN = pANSU + pANSH
        aSVP, aCO2i = INTERNAL_CO2(dTp + aDfSU, vDP, fvPD, c2A, cTyp)
        iPPl, iRDl = PHOTOSYNTHESIS(cTyp, (1. - scPR) * pAR, dTp + aDfSU, aCO2i, sLN - sLMN, eJM, xVN, xJN, thT)
        aRSWU = (ptSU - aTSU) * (slSU * (rbHSU + rT * frSU) + 0.067 * (rbWSU + rT * frSU) / aTSU / 0.067 + ptSU / aTSU * rWSU)
        iAPL = (1.6 * rWSU + 1.3 * rbWSU + rT * frSU) / (1.6 * aRSWU + 1.3 * rbWSU + rT * frSU) * (iPPl - iRDl) + iRDl * (1. - np.exp(-lAI))
        iAP = min(iAPS, (1 - lStr) * iAPS + lStr * iAPL)
        iAPNN = min(iAPN, (1 - lStr) * iAPN + lStr * iAPL * iAPN / iAPS)
        ppC += iPp * wG[i]
        apCS += iAPS * wG[i]
        apCN += iAPNN * wG[i]
        apC += iAP * wG[i]
        ptC += iPpt * wG[i]
        atC += iAt * wG[i]
        peS += iPe * wG[i]
        aeS += iAe  * wG[i]
        dfS += aDfS * wG[i]
        dfSU += aDfSU * wG[i]
        dfSH += aDfSH * wG[i]
        dpR += apR * wG[i]
    ppC *= dLg * 3600
    apCS *= dLg * 3600
    apCN *= dLg * 3600
    apC *= dLg * 3600
    ptC *= dLg * 3600
    atC *= dLg * 3600
    peS *= dLg * 3600
    aeS *= dLg * 3600
    dfS *= dLg * 3600
    dfSU *= dLg * 3600
    dfSH *= dLg * 3600
    dpR *= dLg * 3600
    return ppC, apCS, apCN, apC, ptC, atC, peS, aeS, dfS, dfSU, dfSH, dpR

def SINK_SOURCE(stg, sgf, c_tot, yld_g, f_d, c_dmd_r, c_sup, delta_t):
    c_chg = INSWTICH(stg - sgf, 0., c_tot / yld_g * f_d)
    c_dem = c_chg + max(0, c_dmd_r) / delta_t
    c_flow = min(c_dem, c_sup)
    return c_chg, c_dem, c_flow

def NITROGEN_DYNAMIC(pls, ups, grwt, stnc, minlc, minrc, leafnc, rootnc, leafv, rootv, wleaf, wroot, delta, basec, maxc, tempm, stage, sncd, wstor, lnleaf, lnroot):
    shootn = pls * ups
    leafr = INSWTICH(minlc - leafnc, leafv - wleaf * minlc, 0.) / delta
    rootr = INSWTICH(minrc - rootnc, rootv - wroot * minrc, 0.) / delta
    totalr = leafr + rootr
    stemr = grwt * INSWTICH(-totalr, stnc, 0.)
    cdif = basec + (maxc - basec) * (4 - tempm - stage) / (2 - tempm) * (stage - 1) ** (1 / (2 - tempm))
    esnc = LIMIT_FUNCTION(basec, maxc, cdif) * sncd
    seedsn = shootn - stemr - esnc * wstor
    nextn = max(0., INSWTICH(totalr+seedsn, (totalr+shootn-stemr) / AVOID_ZERO_DIVISION(wstor), esnc))
    snout = wstor * nextn
    leafn = INSWTICH(totalr+seedsn, -leafr-lnleaf, -leafr/AVOID_ZERO_DIVISION(totalr) * (-seedsn) -lnleaf)
    gleaf = INSWTICH(seedsn, leafn, shootn-stemr-snout-lnleaf)
    rleaf = max(-leafv+1E-7, gleaf)
    rtleaf = max(0, rleaf)
    rootn = INSWTICH(totalr+seedsn, ups-shootn-rootr-lnroot, ups-shootn-rootr/AVOID_ZERO_DIVISION(totalr) * (-seedsn) -lnroot)
    groot = INSWTICH(seedsn, rootn, ups-shootn-lnroot)
    rroot = max(-rootv + 5E-8, groot)
    return rroot, stemr, rleaf, rtleaf, snout


def LAI_RATE(dstage, init_sla, rwgt_lv, cur_lai, k_fact, n_lv, n_chg_lv, base_nlv, rate_nlv):
    lai_chg = INSWTICH(rwgt_lv, max(-cur_lai + 1E-5, init_sla * rwgt_lv), init_sla * rwgt_lv)
    if cur_lai < 1 and dstage < 0.5:
        lai_chg = (base_nlv * n_chg_lv - n_lv * rate_nlv) / base_nlv / (base_nlv + k_fact * n_lv)
    return lai_chg

def REFLECTION_COEFF(sc, kbeam):
    kmod = kbeam * (1 - sc)**0.5
    refl = (1 - (1 - sc)**0.5) / (1 + (1 - sc)**0.5)
    pc_b = 1 - np.exp(-2 * refl * kbeam / (1 + kbeam))
    return kmod, pc_b


def PENMAN_MONTEITH(rsw, r_t, r_bw, r_bh, atr_j, atm_tr, fr, t_lf, d_vp, slp, v_pd):
    bzm = 5.668E-8
    l_vap = 2.4E6
    v_ca = 1200
    psyc = 0.067

    clr = max(0, min(1, (atm_tr - 0.25) / 0.45))
    bb_r = bzm * (t_lf + 273)**4
    rlw_n = bb_r * (0.56 - 0.079 * np.sqrt(d_vp * 10)) * (0.1 + 0.9 * clr) * fr

    nr_dc = atr_j - rlw_n

    p_sr = psyc * (r_bw + r_t + rsw) / (r_bh + r_t)

    p_tr = nr_dc * slp / (slp + p_sr) / l_vap

    p_td = (v_ca * v_pd / (r_bh + r_t)) / (slp + p_sr) / l_vap

    p_t = max(1E-10, p_tr + p_td)

    return p_t, nr_dc



def POTENTIAL_PHOTO_TRANSPIRATION(f, dtmp, vp, ca, plant_type, fv, lgt, npk, r_bw, r_bh, r_t, aj, atr, eaj, xv, xj, th):
    sv, fc_i = INTERNAL_CO2(dtmp, vp, fv, ca, plant_type)
    p_lf, l_rd = PHOTOSYNTHESIS(plant_type, lgt, dtmp, fc_i, npk, eaj, xv, xj, th)

    vpd = max(0, sv - vp)
    slp = 4158.6 * sv / (dtmp + 239.)**2
    
    lsw = LEAF_CONDUCTANCE(p_lf, l_rd, dtmp, ca, fc_i, r_bw, r_t)
    pt, nadc = PENMAN_MONTEITH(lsw, r_t, r_bw, r_bh, aj, atr, f, dtmp, vp, slp, vpd)
    tdif = TEMP_DIFF(nadc, pt, r_bh, r_t)
    
    tlf = dtmp + tdif
    
    svl, ci = INTERNAL_CO2(tlf, vp, fv, ca, plant_type)
    lf, rd = PHOTOSYNTHESIS(plant_type, lgt, tlf, ci, npk, eaj, xv, xj, th)
    
    slpl = (svl - sv) / AVOID_ZERO_DIVISION(tlf - dtmp)
    
    lrsw = LEAF_CONDUCTANCE(lf, rd, tlf, ca, ci, r_bw, r_t)
    lpt, lnadc = PENMAN_MONTEITH(lrsw, r_t, r_bw, r_bh, aj, atr, f, tlf, vp, slpl, vpd)
    
    return lf, lpt, lrsw, lnadc, slpl


def THERMAL_TIME(dstage, tmx, tmn, dlt, dl, base, opt, ceil, curv):
    rise = 12. - 0.5 * dl
    set = 12. + 0.5 * dl
    avg_temp = (tmx + tmn) / 2
    tt_acc = 0
    for hr in range(1, 25):
        if rise <= hr <= set:
            temp_hr = avg_temp + dlt + 0.5 * abs(tmx - tmn) * math.cos(0.2618 * (hr - 14))
        else:
            temp_hr = avg_temp + 0.5 * abs(tmx - tmn) * math.cos(0.2618 * (hr - 14))
        temp_hr = min(temp_hr, opt) if dstage > 1 else temp_hr
        if temp_hr < base or temp_hr > ceil:
            temp_unit = 0
        else:
            temp_unit = (((ceil - temp_hr) / (ceil - opt)) * ((temp_hr - base) / (opt - base)) ** ((opt - base) / (ceil - opt))) ** curv
        tt_acc += temp_unit / 24
    daily_temp_util = tt_acc
    return daily_temp_util





def SWTICH3(x, y1, y2, y3):
    return y1 if x < 0 else y2 if x == 0 else y3

def INSWTICH(x, y1, y2):
    return y1 if x < 0 else y2

def LIMIT_FUNCTION(xl, xh, x):
    return max(xl, min(x, xh))

def AVOID_ZERO_DIVISION(x):
    return x if x != 0 else 1.0


def SWITCH(x1, x2):
    return 1.0 if x1 <= 0 and x2 <= 0 else 0.0

