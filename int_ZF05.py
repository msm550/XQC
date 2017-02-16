import numpy as np
import multiprocessing as mp
from scipy.stats import poisson
import matplotlib.pyplot as plt


def mT(A):
    return (A * m_p)


def reduced(m1, m2):
    return (m1 * m2 / (m1 + m2))


def vmin(eR, A):
    return (np.sqrt(mT(A) * eR / 2) / reduced(A * m_p, mH))


def rCos():
    return(2*np.random.random_sample()-1)


def F(eR, A):
    if eR == 0:
        return (1)
    else:
        qF = np.sqrt(2 * mT(A) * eR)
        cF = 1.23 * A ** (1 / 3) - 0.6
        rF = np.sqrt(cF ** 2 + 7 * ((np.pi * aF) ** 2) / 3 - 5 * sF ** 2)
        qrF = qF * rF / 0.197
        return (3 * np.exp((qF * sF / 0.197) ** 2 / 2) * (np.sin(qrF) - np.cos(qrF) * qrF) / qrF ** 3)


def E():
    start = 0
    end = len(eRint)-1
    return (1.52e+15 * eTime * fe * (sum(s[j][1] for j in range(start, end + 1)) - (s[start][1] + s[end][1]) / 2))


def integ(rec):
    vmSi = vmin(rec, Si)
    vmTe = vmin(rec, Te)
    vmHg = vmin(rec, Hg)
    iSi = iTe = iHg = 0
    i = 0
    while v[i][3] < .25:
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
    return(custom[0])


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
    fe=0.5
    N0 = 6.02 * 1e+26
    eRint = [i * 1e-9 for i in range(25, 61)]
    v = np.loadtxt('viE650_sortedLaunch.dat').tolist()
    save = []
    mHSigma = np.loadtxt('Earth_ZF05.dat').tolist()

    for i in range(len(mHSigma)):
        mH=mHSigma[i][0]
        si0=mHSigma[i][1]
        print("Expected spectrum calculation for "+format(mH,'.1f')+" GeV is initiated!")
        pool= mp.Pool(4)
        s=pool.map(integ, eRint)
        s=sorted(s,key=getKey)
        Ei=E()
        save.append(mHSigma[i]+[4.61*si0/Ei])
    np.savetxt("Earth_ZF05.dat",save)
'''
    X2=0
    for i in range(13):
        if O(i)<Ei[i]:
            X2+=(Ei[i]-O(i))**2/Ei[i]
    print("chi_square=", X2)
    po=np.random.poisson(Ei,[100000,13])
    X2list=[]
    X2_counts=0
    for j in range(100000):
        Xx=0
        for i in range(13):
            if po[j,i]<Ei[i]:
                Xx+=(Ei[i]-po[j,i])**2/Ei[i]
        if Xx>=X2:
            X2_counts+=1
    CL=100*(1-X2_counts/100000)
    print("CL=",CL)
'''