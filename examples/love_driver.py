from src.root import root_newton_raphson
import numpy as np
import matplotlib.pyplot as plt

def main():
    p1=1800
    p2=2500
    B1=1900
    B2=3200
    H=4000

    """def f(x,f):
        return (np.tan(2*np.pi*f*x)*x*(p1/p2))**2+x**2-H**2*(B1**-2-B2**-2)
    def f2(x):
        return (np.tan(2*np.pi*2*x)*x*(p1/p2))**2+x**2-H**2*(B1**-2-B2**-2)
    def feq2(x):
        return H*np.sqrt(B1**-2-x**-2)
    def dfdx1(x):
        return 2*x*(np.tan(2*np.pi*x)**2+2*np.pi*x*np.arcsin(2*np.pi*x)**2*np.tan(2*np.pi*x))+1
    def dfdx2(x):
        return"""
    
    def g(x,f):
        return (p1/p2)*(np.sqrt(H**2*(B1**-2-B2**-2)-x**2)/x)-np.tan(2*np.pi*f*x)
    
    def dgdx(x,f):
        return -(np.sqrt((B1**-2-B2**-2)*H**2-x**2)/x**2)-(1/np.sqrt((B1**-2-B2**-2)*H**2-x**2))-(2*np.pi*f*np.arccos(2*np.pi*f*x)**2)
    
    for f in range (0,6):
        g(x,f)
        dgdx(x,f)
    


if __name__=="__main__":
    main()