import numpy as np
import multiprocessing as mp
from scipy.stats import poisson
import matplotlib.pyplot as plt
from scipy.stats import maxwell


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


def bounds(v0):
    print("The spectrum calculation for v0=" + str(int(v0 * 300000)) + " km/s is initiated!")
    s = []
    for rec in eRint:
        vmSi = vmin(rec, Si)
        vmTe = vmin(rec, Te)
        vmHg = vmin(rec, Hg)
        Norm = v0 ** 3 * np.sqrt(np.pi) / 4 * maxwell.cdf(vesc * np.sqrt(2) / v0, loc=0, scale=1)
        if vmSi < vesc:
            iSi = fr * v0 ** 2 / 2 * (np.exp(-(vmSi / v0) ** 2) - np.exp(-(vesc / v0) ** 2)) / Norm
        else:
            iSi = 0
        if vmTe < vesc:
            iTe = fr * v0 ** 2 / 2 * (np.exp(-(vmTe / v0) ** 2) - np.exp(-(vesc / v0) ** 2)) / Norm
        else:
            iTe = 0
        if vmHg < vesc:
            iHg = fr * v0 ** 2 / 2 * (np.exp(-(vmHg / v0) ** 2) - np.exp(-(vesc / v0) ** 2)) / Norm
        else:
            iHg = 0
        s.append([rec, 1.11276e+10 * (si0 * rhoDM / 2 / mH / reduced(mH, m_p) ** 2) * (
        mSi * iSi * (Si * F(rec, Si)) ** 2 + mTe * iTe * (Te * F(rec, Te)) ** 2 + mHg * iHg * (Hg * F(rec, Hg)) ** 2)])

    Ei = []
    for i in range(13):
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
        Ei.append(1.52e+15 * eTime * fe * (sum(s[j][1] for j in range(start, end + 1)) - (s[start][1] + s[end][1]) / 2))
    scale = max(Ei)
    return([300000 * v0, 5 * si0 / scale] + [5 * Ei[i] / scale for i in range(13)])


if __name__ == '__main__':
    n_cores = int(input("# of cores = "))
    shielding_scenario = str(input('body or earth shielding? '))
    if shielding_scenario == 'body':
        fr = (1 - 0.848) / 2
    elif shielding_scenario == 'earth':
        fr = (1 + 0.245) / 2
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
    eRint = [i * 1e-9 for i in range(29, 2506)] + [i * 1e-9 for i in range(4000, 10001)] + [i*1e-9 for i in range(10001,30001)]
    mH = float(input("DM mass in GeV="))
    si0 = 5.5e-31
    vesc = 584 / 300000
    vpeak = [(5 + 5 * i) / 300000 for i in range(20)]+[10 * (i+11) / 300000 for i in range(20)]
    pool1 = mp.Pool(n_cores)
    savet = pool1.map(bounds,vpeak)
    np.savetxt('sigma-v0_' + shielding_scenario +'(mH=' + format(mH, '.1f') + ').dat', savet)
