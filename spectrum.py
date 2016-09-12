import numpy as np
import multiprocessing as mp

def reduced(m1,m2):
    return(m1*m2/(m1+m2))

def r(mH,A):
    return(4*A*m_p*mH/(A*m_p+mH)**2)

def rE(mH,A,randomCos,v):
    return(mH*r(mH,A)*(1-randomCos)*v**2/4)

def rCos():
    return(2*np.random.random_sample()-1)

def recoil(v_i):
    rcos=rCos()
    E_Si=1e+9*rE(mH,Si,rcos,v_i)
    E_Te=1e+9*rE(mH,Te,rcos,v_i)
    E_Hg=1e+9*rE(mH,Hg,rcos,v_i)
    return([E_Si,E_Te,E_Hg])

def which_bin(i):
    nSi=nHg=nTe=0
    for j in 100000:
        if save[j][4]>=bin[i][0] and save[j][4]<bin[i][1]:
            if save[j][5]==Si:
                nSi+=1
            elif save[j][5]==Hg:
                nHg+=1
            else:
                nTe+=1  
    return([nSi/100000,nTe/100000,nHg/100000])
if __name__ == '__main__':  
    Si=28
    Hg=200
    Te=127
    m_p=0.938
    nj=int(1e+6)
    v = np.loadtxt('vi_Erickeck.dat').tolist()
    bin=[[29e-9,36e-9,0],[36e-9,128e-9,11],[128e-9,300e-9,129],[300e-9,540e-9,80],[540e-9,700e-9,90],[700e-9,800e-9,32],[800e-9,945e-9,48],[945e-9,1100e-9,31],[1100e-9,1310e-9,30],[1310e-9,1500e-9,29],[1500e-9,1810e-9,32],[1810e-9,2505e-9,15],[4000e-9,1,60]]
    mH=float(input("DM's mass (GeV)="))            
    pool= mp.Pool(4)    
    save=pool.map(recoil, v[:100000])
    np.savetxt("10_spectrum.dat",np.array(save))   
    pool2=mp.Pool(4)
    Ei=pool2.map(which_bin,range(13))    
    
           
            