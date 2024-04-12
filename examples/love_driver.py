from goph420_lab02.root import root_newton_raphson
import numpy as np
import matplotlib.pyplot as plt

def main():
    p1=1800
    p2=2500
    B1=1900
    B2=3200
    H=4000
    frequencies=[0.01, 0.1, 0.5, 1, 5, 10]

    zeta_arr=[]
    cl_arr=[]
    wavelength_arr=[]

    zeta_max=np.sqrt(H**2 * (B1**-2 - B2**-2)) #1.6939991038925966

    for _, f in enumerate(frequencies):
        k=0
        zeta_a=[]
        zeta_arr_new = []
        cl=[]
        wavelength=[]
        freq=[]
        zeta_a_new = (2*k+1) / (4*f)

        def g(x):
            return (p2/p1)*np.sqrt(zeta_max**2-x**2)/x - np.tan(2*np.pi*f*x)
        def dgdx(x):
            return -(p2/p1)/np.sqrt(zeta_max**2-x**2) - (p2/p1)*np.sqrt(zeta_max**2-x**2)/x**2 - 2*np.pi*f*((1/np.cos(2*np.pi*f*x))**2)
        
        while zeta_a_new < zeta_max:
            zeta_a.append(zeta_a_new)
            k+=1
            zeta_a_new = (2*k+1) / (4*f) #finding asymptotes
        if g(zeta_max) < 0:
            zeta_a.append(zeta_max)

        i=0
        for zeta in zeta_a:
            zeta_0 = zeta - 1e-5
            x, _, eps_a = root_newton_raphson(zeta_0, g, dgdx)
            zeta_arr_new.append(x)
            cl.append(1/np.sqrt(B1**-2-(x/H)**2))
            wavelength.append(cl[i]/f)
            freq.append(f)
            i+=1
        zeta_arr.append(zeta_arr_new)
        cl_arr.append(cl)
        wavelength_arr.append(wavelength)

        modes_wave=[[], [], []]
        modes_cl=[[], [], []]
        modes_f=[[],[],[]]

        for wavelength, cl in zip(wavelength_arr, cl_arr):
            length = len(wavelength)
            if length > 3:
                length = 3
            for i in range (length):
                modes_cl[i].append(cl[i])
                modes_wave[i].append(wavelength[i])

        for i in range (len(modes_f)):
            mode_length = -len(modes_cl[i])
            for freq in frequencies[mode_length:]:
                modes_f[i].append(freq)
    
    plt.figure(figsize=(6, 8))
    i=0
    for f, cl in zip(modes_f, modes_cl):
        plt.plot(f, cl, label=f"mode{i}")
        i+=1
    plt.legend()
    plt.savefig("../figures/cl_modes.png")
    plt.close("all")

    i=0
    for f, wave in zip(modes_f, modes_wave):
        plt.plot(f, wave, label=f"mode{i}")
        i+=1
    plt.legend()
    plt.savefig("../figures/wave_modes.png")

if __name__=="__main__":
    main()