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
    
    # specify thickness perturbation and axon tract stiffness
    data.index = [5,4.75,4.5,4.25,4,3.75,3.5,3.25,3,2.75,2.5]
    data.columns = [10,20,30,40,50,60,70,80,90,100]
    
    # plot heatmap
    g = sns.heatmap(data, cmap='Blues', square=True, linewidths=1, linecolor='white', 
                    cbar_kws ={'label': 'mean squared displacement $\psi$ [mm]','location': 'bottom','shrink':0.51,'ticks': [0,0.5,1.0,1.5,2,2.5]} )
    
    # axes
    g.tick_params(left=False, bottom=False)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_xticklabels(g.get_xticklabels(), rotation=0)
    g.set_xlabel("axon tract stiffness K [N/m]")
    g.set_ylabel(r'hickness perturbation  $\xi/H_{\rm c}$ [\%]')
    plt.tight_layout()
    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500

    # save figure to TIFF if needed
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_psi_heatmap_thickness_perturbation.tiff')
    png1.close()
    
    plt.savefig("figure_psi_heatmap_thickness_perturbation.png")
    plt.show()

if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    heat_map("psi_array_perturbation.csv")