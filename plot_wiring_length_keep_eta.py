from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
from PIL import Image
from io import BytesIO

if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    cmap = mpl.colormaps['Blues']
    color1 = cmap(0.3)
    color2 = cmap(0.6)
    color3 = cmap(0.9)

    growth = np.arange(1,1.23125,0.00125)

    with open('Job-K1-100-K2-100-total.npy', 'rb') as f:
        total1 = np.load(f)
    with open('Job-K1-200-K2-200-total.npy', 'rb') as f:
        total2 = np.load(f)
    with open('Job-K1-300-K2-300-total.npy', 'rb') as f:
        total3 = np.load(f)

    plt.figure()
    plt.plot(growth,total1/total1[0],color=color1,linestyle='solid')
    plt.plot(growth,total2/total2[0],color=color2,linestyle='solid')
    plt.plot(growth,total3/total3[0],color=color3,linestyle='solid')

    plt.legend([ r'$K_{\rm 1} = K_{\rm 2} = 100$  \rm  N/m', r'$K_{\rm 1} = K_{\rm 2} = 200$  \rm  N/m', r'$K_{\rm 1} = K_{\rm 2} = 300$  \rm  N/m'])

    plt.xlabel(r"$\vartheta_{g}$ [-]")
    plt.ylabel(r"normalized length $\overline{L}(t)$ [-]")

    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500

    # save figure to TIFF if needed
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_total_wiring_length_keep_eta.tiff')
    png1.close()

    plt.savefig("figure_total_wiring_length_keep_eta.png")
    plt.show()