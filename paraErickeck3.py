import numpy as np
import multiprocessing as mp
from scipy.stats import poisson
import matplotlib.pyplot as plt

def mT(A):
    return(A*m_p)

def reduced(m1,m2):
    return(m1*m2/(m1+m2))

def r(mH,A):
    return(4*A*m_p*mH/(A*m_p+mH)**2)

def F2_2(eR,A):
    qF=np.sqrt(2*mT(A)*eR)
    cF=1.23*A**(1/3)-0.6
    rF=np.sqrt(cF**2+7*((np.pi*aF)**2)/3-5*sF**2)
    qrF=qF*rF/0.197
    if qrF<=2:
        return(np.exp(-(0.2+(0.9/rF)**2)*qrF**2))
    else:
        return(9*0.81*np.exp(-(0.9*qrF/rF)**2)/qrF**4)
        
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

def pA(lambdainv,leff):
    return(lambdainv*leff)

def rE(mH,A,randomCos,v):
    return(mH*r(mH,A)*(1-randomCos)*v**2/4)

def rCos():
    return(2*np.random.random_sample()-1)

def phi():
    return(2*np.pi*np.random.random_sample())

def r01():
    return(np.random.random_sample())

def a_sel(ra,pASi,pATe):
    if ra>=0 and ra<pASi:
        return(Si) 
    elif ra>=pASi and ra<1-pATe:
        return(Hg)
    else:
        return(Te)

def O(i):
    return(bin[i][2])

def E(i):
    nSi=nHg=nTe=0
    if i==0:
        fe=0.3815
    elif i==1:
        fe=0.5083
    else:
        fe=1
    sSi=[]
    sHg=[]
    sTe=[]
    for j in range(length):
        if save[j][4]>=bin[i][0] and save[j][4]<bin[i][1]:
            if save[j][5]==Si:
                nSi+=1
                sSi.append(save[j])
            elif save[j][5]==Hg:
                nHg+=1
                sHg.append(save[j])
            else:
                nTe+=1
                sTe.append(save[j])
    sigmavSi=sum(sSi[i][0]*sSi[i][3] for i in range(nSi))
    sigmavHg=sum(sHg[i][0]*sHg[i][3] for i in range(nHg))
    sigmavTe=sum(sTe[i][0]*sTe[i][3] for i in range(nTe))
    return(3e+7*eTime*N0*fe*rhoDM/nj/mH*(sigmavSi*mSi/Si+sigmavHg*mHg/Hg+sigmavTe*mTe/Te))
       
def scattering(v_ini):
    if v_ini[3]<-.25:
        return(0)
    else:    
        ra=r01()
        A=a_sel(ra,pASi,pATe)
        CosXiCM=rCos()
        rec=rE(mH,A,CosXiCM,v_ini[0])
        if rec<=4e-6 and rec>=2.505e-6:
            return(0)
        elif rec<=29e-9:
            return(0)
        else:
            return([v_ini[0],mH,si0,si(mH,si0,A,rec),rec,A,v_ini[1]])
    
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
    mH=float(input("DM's mass (GeV)="))
    si0=float(input("cross-section (ubarn)="))*1e-30
    rhoDM=0.3
    aF=0.52
    sF=0.9
    N0=6.02*1e+26
    eRec=0
    bin=[[29e-9,36e-9,0],[36e-9,128e-9,11],[128e-9,300e-9,129],[300e-9,540e-9,80],[540e-9,700e-9,90],[700e-9,800e-9,32],[800e-9,945e-9,48],[945e-9,1100e-9,31],[1100e-9,1310e-9,30],[1310e-9,1500e-9,29],[1500e-9,1810e-9,32],[1810e-9,2505e-9,15],[4000e-9,1,60]]
    l_inv_Si=lambdainv(mH,si0,Si,eRec)
    l_inv_Hg=lambdainv(mH,si0,Hg,eRec)
    l_inv_Te=lambdainv(mH,si0,Te,eRec)
    leff=lambdaeff(l_inv_Si,l_inv_Hg,l_inv_Te)
    pASi=pA(l_inv_Si,leff)
    pATe=pA(l_inv_Te,leff)
    nj=int(1e+6)
    v = np.loadtxt('vi_Erickeck_3.dat').tolist()
    pool1= mp.Pool(4)    
    save=[]
    s=pool1.map(scattering, v)
    for j in range(nj):
        if s[j]!=0:
            save.append(s[j])
    length=len(save)
    pool2 = mp.Pool(4)
    Ei=pool2.map(E,range(13))
    X2=0
    c_X2=0
    Likelihood=0
    for i in range(13):
        if O(i)<Ei[i]:
            X2+=(Ei[i]-O(i))**2/Ei[i]
            c_X2+=1
            Likelihood+=-poisson.logpmf(O(i), Ei[i], loc=0)
    print("chi_square=", X2)
    print("chi_square_weighted=", X2/(c_X2+1))
    print("Likelihood=",Likelihood)
    po=np.random.poisson(Ei,[100000,13])
    X2list=[]
    X2list_weighted=[]
    X2_counts=0
    X2_counts_weighted=0
    Li_counts=0
    for j in range(100000):
        Xx=0
        c_Xx=0
        Likelihoodx=0
        for i in range(13):
            if po[j,i]<Ei[i]:
                Xx+=(Ei[i]-po[j,i])**2/Ei[i]
                c_Xx+=1
                Likelihoodx+=-poisson.logpmf(po[j,i], Ei[i], loc=0)
        if Xx>=X2:
            X2_counts+=1
        if Xx/(c_Xx+1)>=X2/(c_X2+1):
            X2_counts_weighted+=1
        if Likelihoodx>=Likelihood:
            Li_counts+=1
    CL=100*(1-X2_counts/100000)
    CL_weighted=100*(1-X2_counts_weighted/100000)
    CL_Li=100*(1-Li_counts/100000)
    print("CL=",CL)
    print("CL_weighted=",CL_weighted)
    print("CL_L=",CL_Li)
            
         
