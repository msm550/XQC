import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy.stats import poisson
import matplotlib.pyplot as plt


def mT(A):
    return (A * m_p)


def reduced(m1, m2):
    return (m1 * m2 / (m1 + m2))


def vmin(eR, A):
    return (np.sqrt(mT(A) * eR / 2) / reduced(A * m_p, mH))


def F(eR, A):
    if eR == 0:
        return (1)
    else:
        qF = np.sqrt(2 * mT(A) * eR)
        cF = 1.23 * A ** (1 / 3) - 0.6
        rF = np.sqrt(cF ** 2 + 7 * ((np.pi * aF) ** 2) / 3 - 5 * sF ** 2)
        qrF = qF * rF / 0.197
        return (3 * np.exp((qF * sF / 0.197) ** 2 / 2) * (np.sin(qrF) - np.cos(qrF) * qrF) / qrF ** 3)


def O(i):
    return (bin[i][2])


def E(i):
    if i == 0:
        fe = 0.3815
    elif i == 1:
        fe = 0.5083
    else:
        fe = 1
    if i == 12:
        start = int(bin[i][0] * 1e+9) - 1523
        end = len(s) - 1
    else:
        start = int(bin[i][0] * 1e+9) - 29
        end = int(bin[i][1] * 1e+9) - 29
    return (1.52e+15 * eTime * fe * (sum(s[j][1] for j in range(start, end + 1)) - (s[start][1] + s[end][1]) / 2))


def integ(rec):
    vmSi = vmin(rec, Si)
    vmTe = vmin(rec, Te)
    vmHg = vmin(rec, Hg)
    iSi = iTe = iHg = 0
    i = 0
    while v[i][zenith_column_index] < zenith_angle_limit:
        vnorm = v[i][0]
        if vnorm > vmSi:
            iSi += (1 / vnorm)
        if vnorm > vmTe:
            iTe += (1 / vnorm)
        if vnorm > vmHg:
            iHg += (1 / vnorm)
        i += 1
    return ([rec, 1.11276e+4 * (si0 * rhoDM / 2 / mH / reduced(mH, m_p) ** 2) * (
    mSi * iSi * (Si * F(rec, Si)) ** 2 + mTe * iTe * (Te * F(rec, Te)) ** 2 + mHg * iHg * (Hg * F(rec, Hg)) ** 2)])


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
    rhoDM = 0.3
    aF = 0.52
    sF = 0.9
    N0 = 6.02 * 1e+26
    bin = [[29e-9, 36e-9, 0], [36e-9, 128e-9, 11], [128e-9, 300e-9, 129], [300e-9, 540e-9, 80], [540e-9, 700e-9, 90],
           [700e-9, 800e-9, 32], [800e-9, 945e-9, 48], [945e-9, 1100e-9, 31], [1100e-9, 1310e-9, 30],
           [1310e-9, 1500e-9, 29], [1500e-9, 1810e-9, 32], [1810e-9, 2505e-9, 15], [4000e-9, 1, 60]]
    eRint = [i * 1e-9 for i in range(29, 2506)] + [i * 1e-9 for i in range(4000, 10001)]
    eRExtention = [i * 1e-9 for i in range(10001, 30001)]

    shielding_scenario = str(input("<body> or <earth> shielding ="))
    if shielding_scenario == 'earth':
        v_fn = 'viE_sortedLaunch.csv'
        zenith_column_index = 3
        zenith_angle_limit = .245
        output_fn = 'earth_integralMethod.dat'
    elif shielding_scenario == 'body':
        v_fn = 'viE_sortedDet.csv'
        zenith_column_index = 4
        zenith_angle_limit = -.848
        output_fn = 'body_integralMethod.dat'

    v_df = pd.read_csv(v_fn)
    v = [v_row[0] for v_row in v_df.values]
    save = []
    try:
        mHSigma = np.loadtxt(output_fn).tolist()
        for i in range(len(mHSigma)):
            mH = mHSigma[i][0]
            si0 = mHSigma[i][2]
            print("Expected spectrum calculation for " + format(mH, '.1f') + " GeV is initiated!")
            pool1 = mp.Pool(4)
            if mH > 10:
                s = pool1.map(integ, eRint + eRExtention)
            else:
                s = pool1.map(integ, eRint)
            s = sorted(s, key=getKey)
            pool2 = mp.Pool(4)
            Ei = pool2.map(E, range(13))
            save.append([mH, si0] + Ei)

    except FileNotFoundError:
        mH_list = np.logspace(np.log10(0.3), 2, 27)
        for i in range(27):
            mH = mH_list[i]
            si0 = 1e-29
            print("Expected spectrum calculation for " + format(mH, '.1f') + " GeV is initiated!")
            pool1 = mp.Pool(4)
            if mH > 10:
                s = pool1.map(integ, eRint + eRExtention)
            else:
                s = pool1.map(integ, eRint)
            s = sorted(s, key=getKey)
            pool2 = mp.Pool(4)
            Ei = pool2.map(E, range(13))
            scale = max(Ei)
            save.append([mH, 5 * si0 / scale] + [5 * Ei[i] / scale for i in range(13)])
    np.savetxt(output_fn, save)