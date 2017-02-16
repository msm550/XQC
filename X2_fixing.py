import numpy as np
import multiprocessing as mp

def O(i):
    return(bin[i][2])
    
def CL_calc(Ei):
    X2=0
    for i in range(13):
        if O(i)<Ei[i]:
            X2+=(Ei[i]-O(i))**2/Ei[i]
    po=np.random.poisson(Ei,[100000,13])
    X2_counts=0
    for j in range(100000):
        Xx=0
        for i in range(13):
            if po[j,i]<Ei[i]:
                Xx+=(Ei[i]-po[j,i])**2/Ei[i]
        if Xx>=X2:
            X2_counts+=1
    CL=100*(1-X2_counts/100000)
    return(CL)

def search(Ei):
    CL=CL_calc(Ei[2:])
    print("The search to fix the CL for mass "+format(Ei[0],'.1f')+" GeV is initiated!")    
    while CL>91 or CL<89:
        Ei=[Ei[0]]+[(1+(90-CL)*(180-CL)/100/180)*Ei[i] for i in range(1,15)]
        CL=CL_calc(Ei[2:])
    print("mass "+format(Ei[0],'.1f')+" GeV's cross-section is "+format(Ei[1],'.2g')+" cm^2")
    return(Ei)
    
if __name__ == '__main__':
    bin=[[29e-9,36e-9,0],[36e-9,128e-9,11],[128e-9,300e-9,129],[300e-9,540e-9,80],[540e-9,700e-9,90],[700e-9,800e-9,32],[800e-9,945e-9,48],[945e-9,1100e-9,31],[1100e-9,1310e-9,30],[1310e-9,1500e-9,29],[1500e-9,1810e-9,32],[1810e-9,2505e-9,15],[4000e-9,1,60]]
    filename = 'body_integralMethod.dat'
    E= np.loadtxt(filename).tolist()
    pool1=mp.Pool(4)
    save=pool1.map(search,E)
    np.savetxt(filename,save)
