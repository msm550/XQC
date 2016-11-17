import numpy as np

def reduced(m1,m2):
    return(m1*m2/(m1+m2))

def r(mH,A):
    return(4*A*m_p*mH/(A*m_p+mH)**2)

def rE(mH,A,randomCos,v):
    return(mH*r(mH,A)*(1-randomCos)*v**2/4)

def rCos():
    return(2*np.random.random_sample()-1)
    
if __name__ == '__main__':  
    Si=28
    Hg=200
    Te=127
    m_p=0.938
    nj=int(1e+6)
    v = np.loadtxt('vi_Erickeck.dat').tolist()
    Count_Si=0
    Count_Te=0
    Count_Hg=0
    mH=float(input("DM's mass (GeV)="))            

    for i in range(nj):
        rcos=rCos()
        E_Si=rE(mH,Si,rcos,v[i])
        E_Te=rE(mH,Te,rcos,v[i])
        E_Hg=rE(mH,Hg,rcos,v[i])
        if E_Si>=2.9e-8 and E_Si<=1.28e-7:
            Count_Si+=1
        if E_Te>=2.9e-8 and E_Te<=1.28e-7:
            Count_Te+=1
        if E_Hg>=2.9e-8 and E_Hg<=1.28e-7:
            Count_Hg+=1
    print("fSi*=",Count_Si/nj)
    print("fTe*=",Count_Te/nj)
    print("fHg*=",Count_Hg/nj)
    
           
            