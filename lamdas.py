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
    if A==Ox:
        return(0.465)
    elif A==Si:
        return(0.289)
    elif A==Al:
        return(0.089)
    else:
        return(0.048)
    
if __name__ == '__main__':
    Si=28
    Ox=16
    Al=27
    Fe=56
    Pb=207
    Cu=63
    C=12
    H=1
    rhoAl=2.7
    rhoFe=7.87
    aF=0.52
    sF=0.9
    m_p=0.938
    rhoE=2.7
    rhoDM=0.3
    N0=6.02*1e+26
    s= np.loadtxt('XQCE_13.dat').tolist()
    sE=np.loadtxt('XQCE_read.dat').tolist()
    l=[]    
    for i in range(len(s)):
        mH=s[i][0]
        sigmap=s[i][1]
        l_Al=(N0*rhoAl*si(mH,sigmap,Al,0)/1000/Al)**-1
        l.append([mH,sigmap,l_Al,1-np.exp(-1/l_Al),r(mH,Al)/2])
    lE=[]
    for i in range(len(sE)):
        mH=sE[i][0]
        sigmap=sE[i][1]
        l_Al=(N0*rhoAl*si(mH,sigmap,Al,0)/1000/Al)**-1
        lE.append([mH,sigmap,l_Al,1-np.exp(-1/l_Al),r(mH,Al)/2])

        