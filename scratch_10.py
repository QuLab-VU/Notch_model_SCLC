from pysb import *
from pysb.simulator.scipyode import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator
import numpy as np
import matplotlib.pyplot as plt
import copy
import seaborn as sns
import matplotlib.colors
from mpl_toolkits.mplot3d import Axes3D


kahp_multiplier = np.arange(0.5, 1.0, 0.5)

for m, mult in enumerate(kahp_multiplier):
    print('This is m: {0:.2f} This is mult: {1:.2f}'.format(m, mult))

    plt.subplot(1, len(kahp_multiplier), m + 1, title='KaHp = {0:.2f}* 1 * 0.02535'.format(mult))

    fig, ax = plt.subplots()
    print('In loop, after: fig type: {0}'.format(type(fig)))
    #example_plot(ax, fontsize=24)
    #ax.text(0.65 ,0.35, 'example text', xycoords='axes fraction')
    plt.tight_layout()

    #plt.text(0.5 ,0.5,'Example text.')
    #matplotlib.text('XYZ', (1,1))
    #ax.text('example', (0.7,0.3), xycoords='axes fraction')
    #text('KaHp = {0:.2f} * 0.02535'.format(mult), (1,1), xycoords = 'axes fraction')
    #text('XYZ', (1,1))

print(len(kahp_multiplier))












