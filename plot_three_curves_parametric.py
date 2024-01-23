import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
from io import BytesIO
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
plt.rcParams["font.family"] = "Times New Roman"

def heat_map(csv_file):
 
    # read in data and plot
    data = pd.read_csv(csv_file,header=None)

    # specify primary and secondary axon tracts stiffness
    data.columns = [10,120,230,340,450,560,670,780,890,1000]
    data.index = [1000,890,780,670,560,450,340,230,120,10]
    g = sns.heatmap(data, cmap='Blues', square=True, linewidths=1, linecolor='white', 
                    cbar_kws ={'label': 'mean squared displacement $\psi$ [mm]','location': 'bottom','shrink':0.51,'ticks': [0.0,0.4,0.8,1.2,1.6]} )
    
    # axes
    g.tick_params(left=False, bottom=False)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_xticklabels(g.get_xticklabels(), rotation=0)
    g.set_xlabel("secondary axon tract stiffness $K_{\r 2}$ [N/m]")
    g.set_ylabel(r"primary axon tract stiffness $K_{\rm 1}$ [N/m]")
    plt.tight_layout()
    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500

    # save figure to TIFF if needed
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_heatmap_three_curves.tiff')
    png1.close()
    
    plt.savefig("figure_heatmap_three_curves.png")
    plt.show()


if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    heat_map("psi_array_three_curves.csv")