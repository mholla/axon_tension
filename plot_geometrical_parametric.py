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
    
    # specify span and tangent
    data.index = [25,24,23,22,21,20,19,19,17,16]
    data.columns = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]

    # plot heatmap
    g = sns.heatmap(data, cmap='Blues', square=True, linewidths=1, linecolor='white', 
                    cbar_kws ={'label': 'mean squared displacement $\psi$ [mm]','location': 'bottom','shrink':0.5,'ticks': [0,0.5,1.0,1.5,2,2.5]})

    # axes
    g.tick_params(left=False, bottom=False)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_xticklabels(g.get_xticklabels(), rotation=0)
    g.set_xlabel(r"slope of the parabola root "'tan$( \\theta)$' r" [-]")
    g.set_ylabel(r"span of the parabola s [mm]")
    plt.tight_layout()
    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500

    # save figure to TIFF if needed
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_heatmap_geometrical.tiff')
    png1.close()
    
    plt.savefig("figure_heatmap_geometrical.png")
    plt.show()


if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    heat_map("psi_array_geometry.csv")