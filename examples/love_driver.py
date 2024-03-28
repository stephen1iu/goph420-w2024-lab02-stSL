from goph420_lab02.root import root_newton_raphson
import numpy as np
import matplotlib.pyplot as plt

def main():
    p1=1800
    p2=2500
    B1=1900
    B2=3200
    H=4000

    cl=[]
    wavelength=[]
    frequencies=[0.01, 0.1, 0.5, 1, 5, 10]
    zeta_arr=[]

    zeta_max=np.sqrt(H**2 * (B1**-2 - B2**-2)) #1.6939991038925966

    for f in frequencies:
        k=0
        x0=(2*k+1) / (4*f)
        zeta_a=[]

        def g(x):
            return (p2/p1)*(np.sqrt(H**2*(zeta_max**2-x**2))/x)-np.tan(2*np.pi*f*x)
        def dgdx(x):
            return -((p2*(zeta_max**2))/(p1*(x**2)*np.sqrt((zeta_max**2-x**2)))) - (2*np.pi*f*((1/np.cos(2*np.pi*f*x))**2))
        print (f"freq {f}")
        
        while x0 < zeta_max:
            zeta_a.append(x0)
            k+=1
            x0=(2*k+1) / (4*f) - 1e-5

        print (zeta_a)

        for zeta in zeta_a:
            zeta_0 = zeta
            x, _, eps_a = root_newton_raphson(zeta_0, g, dgdx)
            zeta_arr.append(x)
        # cl.append(1/np.sqrt(B1**-2-(x/H)**2))
        # wavelength.append(cl/f)

if __name__=="__main__":
    main()
