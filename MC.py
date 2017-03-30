import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy.stats import poisson
import matplotlib.pyplot as plt


def mT(A):
    return (A * m_p)


def reduced(m1, m2):
    return (m1 * m2 / (m1 + m2))


def r(mH, A):
    return (4 * A * m_p * mH / (A * m_p + mH) ** 2)


def F2_2(eR, A):
    qF = np.sqrt(2 * mT(A) * eR)
    cF = 1.23 * A ** (1 / 3) - 0.6
    rF = np.sqrt(cF ** 2 + 7 * ((np.pi * aF) ** 2) / 3 - 5 * sF ** 2)
    qrF = qF * rF / 0.197
    if qrF <= 2:
        return (np.exp(-(0.2 + (0.9 / rF) ** 2) * qrF ** 2))
    else:
        return (9 * 0.81 * np.exp(-(0.9 * qrF / rF) ** 2) / qrF ** 4)


def F(eR, A):
    if eR == 0:
        return (1)
    else:
        qF = np.sqrt(2 * mT(A) * eR)
        cF = 1.23 * A ** (1 / 3) - 0.6
        rF = np.sqrt(cF ** 2 + 7 * ((np.pi * aF) ** 2) / 3 - 5 * sF ** 2)
        qrF = qF * rF / 0.197
        return (3 * np.exp((qF * sF / 0.197) ** 2 / 2) * (np.sin(qrF) - np.cos(qrF) * qrF) / qrF ** 3)


def si(mH, si0, A, eR):
    return (si0 * (A * reduced(A * m_p, mH) * F(eR, A) / reduced(m_p, mH)) ** 2)


def f(A):
    if A == Si:
        return (0.83145)
    elif A == Hg:
        return (0.10308)
    else:
        return (0.06547)


def lambdainv(mH, si0, A, eR):
    return (N0 * rhoTarget * si(mH, si0, A, eR) * f(A) / 1000 / A)


def lambdaeff(lambdainvSi, lambdainvHg, lambdainvTe):
    return ((lambdainvSi + lambdainvHg + lambdainvTe) ** -1)


def pA(lambdainv, leff):
    return (lambdainv * leff)


def rE(mH, A, randomCos, v):
    return (mH * r(mH, A) * (1 - randomCos) * v ** 2 / 4)


def rCos():
    return (2 * np.random.random_sample() - 1)


def phi():
    return (2 * np.pi * np.random.random_sample())


def r01():
    return (np.random.random_sample())


def a_sel(ra, pASi, pATe):
    if ra >= 0 and ra < pASi:
        return (Si)
    elif ra >= pASi and ra < 1 - pATe:
        return (Hg)
    else:
        return (Te)


def O(i):
    return (bin[i][2])


def E(i):
    nSi = nHg = nTe = 0
    if i == 0:
        fe = 0.3815
    elif i == 1:
        fe = 0.5083
    else:
        fe = 1
    sSi = []
    sHg = []
    sTe = []
    for j in range(length):
        if save[j][4] >= bin[i][0] and save[j][4] < bin[i][1]:
            if save[j][5] == Si:
                nSi += 1
                sSi.append(save[j])
            elif save[j][5] == Hg:
                nHg += 1
                sHg.append(save[j])
            else:
                nTe += 1
                sTe.append(save[j])
    sigmavSi = sum(sSi[i][0] * sSi[i][3] for i in range(nSi))
    sigmavHg = sum(sHg[i][0] * sHg[i][3] for i in range(nHg))
    sigmavTe = sum(sTe[i][0] * sTe[i][3] for i in range(nTe))
    return (3e+7 * eTime * N0 * fe * rhoDM / nj / mH * (
    sigmavSi * mSi / Si / pASi + sigmavHg * mHg / Hg / pAHg + sigmavTe * mTe / Te / pATe))


def scattering(v_ini):
    if v_ini[zenith_column_index] < zenith_angle_limit:
        return (0)
    else:
        ra = r01()
        A = a_sel(ra, pASi, pATe)
        CosXiCM = rCos()
        rec = rE(mH, A, CosXiCM, v_ini[0])
        if rec <= 4e-6 and rec >= 2.505e-6:
            return (0)
        elif rec <= 29e-9:
            return (0)
        else:
            return ([v_ini[0], mH, si0, si(mH, si0, A, rec), rec, A])


def getKey(custom):
    return (custom[0])


if __name__ == '__main__':
    Si = 28
    Hg = 200
    Te = 127
    rhoSi = 2.329
    rhoHgTe = 8.1
    volSi = 34 * (14 * 0.5 * 2 + 12 * 0.245 ** 2 + 7 * 0.25) * 1e-6
    volHgTe = 34 * (0.96 * 0.5 * 2) * 1e-6
    mSi = rhoSi * volSi
    mHg = rhoHgTe * volHgTe * Hg / (Te + Hg)
    mTe = rhoHgTe * volHgTe * Te / (Te + Hg)
    rhoTarget = (rhoSi * volSi + rhoHgTe * volHgTe) / (volHgTe + volSi)
    eTime = 100.7
    m_p = 0.938
    nj = int(1e+6)
    rhoDM = 0.3
    aF = 0.52
    sF = 0.9
    N0 = 6.02 * 1e+26
    savet = []
    bin = [[29e-9, 36e-9, 0], [36e-9, 128e-9, 11], [128e-9, 300e-9, 129], [300e-9, 540e-9, 80], [540e-9, 700e-9, 90],
           [700e-9, 800e-9, 32], [800e-9, 945e-9, 48], [945e-9, 1100e-9, 31], [1100e-9, 1310e-9, 30],
           [1310e-9, 1500e-9, 29], [1500e-9, 1810e-9, 32], [1810e-9, 2505e-9, 15], [4000e-9, 1, 60]]

    shielding_scenario = str(input("<body> or <earth> shielding ="))
    if shielding_scenario == 'earth':
        v_fn = 'viE_sortedLaunch.csv'
        zenith_column_index = 3
        zenith_angle_limit = .245
        output_fn = 'Ch2_earth.dat'
    elif shielding_scenario == 'body':
        v_fn = 'viE_sortedDet.csv'
        zenith_column_index = 4
        zenith_angle_limit = -.848
        output_fn = 'Chi2_body.dat'

    v_df = pd.read_csv(v_fn)
    v = [v_row[0] for v_row in v_df.values]

    try:
        mHSigma = np.loadtxt(output_fn).tolist()
        for i in range(len(mHSigma)):
            eRec = 0
            mH = mHSigma[i][0]
            si0 = mHSigma[i][2]
            print("Expected spectrum calculation for " + format(mH, '.1f') + " GeV is initiated!")
            l_inv_Si = lambdainv(mH, si0, Si, eRec)
            l_inv_Hg = lambdainv(mH, si0, Hg, eRec)
            l_inv_Te = lambdainv(mH, si0, Te, eRec)
            leff = lambdaeff(l_inv_Si, l_inv_Hg, l_inv_Te)
            pASi = pA(l_inv_Si, leff)
            pATe = pA(l_inv_Te, leff)
            pAHg = 1 - pASi - pATe
            pool1 = mp.Pool(4)
            save = []
            s = pool1.map(scattering, v)
            for j in range(nj):
                if s[j] != 0:
                    save.append(s[j])
            length = len(save)
            pool2 = mp.Pool(4)
            Ei = pool2.map(E, range(13))
            savet.append([mH, si0] + Ei)
    except FileNotFoundError:
        mH_list = np.logspace(np.log10(0.3), 2, 27)
        for i in range(27):
            mH = mH_list[i]
            si0 = 1e-29
            eRec = 0
            print("Expected spectrum calculation for " + format(mH, '.1f') + " GeV is initiated!")
            l_inv_Si = lambdainv(mH, si0, Si, eRec)
            l_inv_Hg = lambdainv(mH, si0, Hg, eRec)
            l_inv_Te = lambdainv(mH, si0, Te, eRec)
            leff = lambdaeff(l_inv_Si, l_inv_Hg, l_inv_Te)
            pASi = pA(l_inv_Si, leff)
            pATe = pA(l_inv_Te, leff)
            pAHg = 1 - pASi - pATe
            pool1 = mp.Pool(4)
            save = []
            s = pool1.map(scattering, v)
            for j in range(nj):
                if s[j] != 0:
                    save.append(s[j])
            length = len(save)
            pool2 = mp.Pool(4)
            Ei = pool2.map(E, range(13))
            scale = max(Ei)
            save.append([mH, 5 * si0 / scale] + [5 * Ei[i] / scale for i in range(13)])

    np.savetxt(filename, savet)
