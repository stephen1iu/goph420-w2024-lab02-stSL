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
    cl_arr=[]
    wavelength_arr=[]

    zeta_max=np.sqrt(H**2 * (B1**-2 - B2**-2)) #1.6939991038925966

    for f in frequencies:
        k=0
        zeta_a=[]
        zeta_arr_new = []
        cl=[]
        wavelength=[]
        zeta_a_new = (2*k+1) / (4*f)

        plt.figure()

        def g(x):
            return (p2/p1)*np.sqrt(zeta_max**2-x**2)/x - np.tan(2*np.pi*f*x)
        def dgdx(x):
            return -(p2/p1)/np.sqrt(zeta_max**2-x**2) - (p2/p1)*np.sqrt(zeta_max**2-x**2)/x**2 - 2*np.pi*f*((1/np.cos(2*np.pi*f*x))**2)
        
        while zeta_a_new < zeta_max:
            zeta_a.append(zeta_a_new)
            k+=1
            zeta_a_new = (2*k+1) / (4*f)
        if g(zeta_max) < 0:
            zeta_a.append(zeta_max)

        i=0
        for zeta in zeta_a:
            zeta_0 = zeta - 1e-5
            x, _, eps_a = root_newton_raphson(zeta_0, g, dgdx)
            zeta_arr_new.append(x)
            cl.append(1/np.sqrt(B1**-2-(x/H)**2))
            wavelength.append(cl[i]/f)
            i+=1
        zeta_arr.append(zeta_arr_new)
        cl_arr.append(cl)
        wavelength_arr.append(wavelength)

        modes_wave=[[], [], []]
        for wavelength in wavelength_arr:
            length=len(wavelength)
            if length > 3:
                length = 3
            for i in range(length):
                modes_wave[i].append(wavelength[i])

        modes_cl=[[], [], []]
        for cl in cl_arr:
            length=len(cl)
            if length > 3:
                length = 3
            for i in range(length):
                modes_cl[i].append(cl[i])

        def tan(zeta, f):
            return np.tan(2*np.pi*f*zeta)

        for i, zeta in enumerate (zeta_arr):
            x = np.linspace(1e-4, zeta_max, 10)
            y = tan(x, f)
            #print (y)
            for z in zeta:
                plt.plot(z, tan(z, f), "or") #plots point of root
                plt.plot(x, y, "-b") #plots graph of tan function
                plt.ylim((-10, 10))
                plt.savefig(f"../figures/modes for frequency {f}.png")

if __name__=="__main__":
    main()