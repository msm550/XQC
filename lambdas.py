import numpy as np

def mT(A):
    return(A*m_p)

def reduced(m1,m2):
    return(m1*m2/(m1+m2))

def r(mH,A):
    return(4*A*m_p*mH/(A*m_p+mH)**2)

def F(eR,A):
    if eR==0:
        return(1)
    else:
        qF=np.sqrt(2*mT(A)*eR)
        cF=1.23*A**(1/3)-0.6
        rF=np.sqrt(cF**2+7*((np.pi*aF)**2)/3-5*sF**2)
        qrF=qF*rF/0.197
        return(3*np.exp((qF*sF/0.197)**2/2)*(np.sin(qrF)-np.cos(qrF)*qrF)/qrF**3)

def si(mH,si0,A,eR):
    return(si0*(A*reduced(A*m_p,mH)*F(eR,A)/reduced(m_p,mH))**2)

def f(A):
    if A==Si:
        return(0.83145)
    elif A==Hg:
        return(0.10308)
    else:
        return(0.06547)

def lambdainv(mH,si0,A,eR):
    return(N0*rhoTarget*si(mH,si0,A,eR)*f(A)/1000/A)

def lambdaeff(lambdainvSi,lambdainvHg,lambdainvTe):
    return((lambdainvSi+lambdainvHg+lambdainvTe)**-1)
    
if __name__ == '__main__':
    Si=28
    Hg=200
    Te=127
    rhoSi=2.329
    rhoHgTe=8.1
    volSi=34*(14*0.5*2+12*0.245**2+7*0.25)*1e-6
    volHgTe=34*(0.96*0.5*2)*1e-6
    mSi=rhoSi*volSi
    mHg=rhoHgTe*volHgTe*Hg/(Te+Hg)
    mTe=rhoHgTe*volHgTe*Te/(Te+Hg)
    rhoTarget=(rhoSi*volSi+rhoHgTe*volHgTe)/(volHgTe+volSi)
    eTime=100.7
    m_p=0.938
    rhoDM=0.3
    aF=0.52
    sF=0.9
    N0=6.02*1e+26
    file_name="body_integralMethod.dat"
    mHSigma=np.loadtxt(file_name).tolist()
    s=[]
    
    for i in range(len(mHSigma)):
        eRec=0        
        mH=mHSigma[i][0]
        si0=mHSigma[i][1]
        l_inv_Si=lambdainv(mH,si0,Si,eRec)
        l_inv_Hg=lambdainv(mH,si0,Hg,eRec)
        l_inv_Te=lambdainv(mH,si0,Te,eRec)
        leff=lambdaeff(l_inv_Si,l_inv_Hg,l_inv_Te)
        s.append(mHSigma[i]+[leff])
    np.savetxt(file_name,s)

    