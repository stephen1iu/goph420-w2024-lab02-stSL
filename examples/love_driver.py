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
    frequencies=[0.01, 0.1, 0.5, 1, 2, 5]

    zeta_max=np.sqrt(H**2*(B1**-2-B2**-2))

    #find zeta max
    #based on asymptotes
    #zeta**2 cant be > H^2*(B1^-2-B2^-2) 
    #tan(2 pi f)zeta has asymptotes
    #2k+1 for asymptotes    2k+1(pi/2)=tan(2 pi f)zeta
    #find asymptotes between 0 and zeta max
    #while zeta_a < zeta_max
    #always use primary mode, mode closest to 0
    #use small offset of asymptote for initial guess (guess=asymptote-1e-5)
    k=0
    while zeta_a < zeta_max:
        zeta_a=((2*k+1)*(np.pi/2))/np.tan(2*np.pi*f)
        k+=1

    for f in frequencies:
        def g(x,f):
            return (p1/p2)*(np.sqrt(H**2*(B1**-2-B2**-2)-x**2)/x)-np.tan(2*np.pi*f*x)
        def dgdx(x,f):
            return -(np.sqrt((B1**-2-B2**-2)*H**2-x**2)/x**2)-(1/np.sqrt((B1**-2-B2**-2)*H**2-x**2))-(2*np.pi*f*np.arccos(2*np.pi*f*x)**2)
        
        x, _, eps_a = root_newton_raphson(x, g, dgdx)
        cl.append(1/np.sqrt(B1**-2-(x/H)**2))
        wavelength.append(cl/f)

if __name__=="__main__":
    main()
