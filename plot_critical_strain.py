import pandas as pd
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib as mpl
from PIL import Image
from io import BytesIO
mpl.rcParams.update(mpl.rcParamsDefault)
plt.rcParams["font.family"] = "Times New Roman"

if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    plt.figure()
    plt.axis([0, 30, 0, 0.5])
    plt.xlabel(r"stiffness ratio $\beta$ [-]")
    plt.ylabel(r"critical strain $\varepsilon_{\rm crit}$ [-]")

    # name colors
    cmap = matplotlib.colormaps['Blues']
    color_0 = cmap(0.2)
    color_10 = cmap(0.4)
    color_150 = cmap(0.6)
    color_300 = cmap(0.8)

    # read in and plot data from Holland 2018
    data = pd.read_csv("data_critical-strain.csv")
    betas_theory = data['beta'].tolist()
    eps_theory = data['strain'].tolist()

    plt.plot(betas_theory, eps_theory, color=color_0, linestyle="-")

    betas = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
    eps_10 = [0.0833,0.0556,0.0455,0.0400,0.0361,0.0338,0.0316,0.0304,0.0293,0.0287]
    eps_150 = [0.0755,0.0498,0.0383,0.0316,0.0270,0.0235,0.0218,0.0200,0.0182,0.0165]
    eps_300 = [0.0695,0.0433,0.0333,0.0253,0.0218,0.0176,0.0165,0.0135,0.0129,0.0123]

    plt.plot(betas, eps_10, color=color_10, marker='o', linestyle='')
    plt.plot(betas, eps_150, color=color_150, marker='o', linestyle='')
    plt.plot(betas, eps_300, color=color_300, marker='o', linestyle='')

    plt.legend(['theoretical', r'$K = 10$ \rm  N/m', r'$K = 150$ \rm  N/m', r'$K = 300$ \rm  N/m'])
    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500


    # save figure to TIFF if needed
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_critical-strain.tiff')
    png1.close()

    plt.savefig("figure_critical-strain.png")
    plt.show()