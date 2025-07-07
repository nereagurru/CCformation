# Standarized format for plots

import matplotlib.pyplot as plt


def init_plot(ax, title='', xlabel='', ylabel='', font = 40,
              length1=10, length2=5, tickwidth=1.75):
    """
    Function to create a plot.

    Parameters
    ----------
    ax : Axes
        Axes object to plot in.
    title : str
        Title of the plot
    xlabel : str
        Label of the x axis
    ylabel : str
        Label of the y axis
    font : float 
        Font-size. 40 by default
    Returns
    -------
    None.

    """

    plt.rcParams['axes.linewidth'] = 2.25
    plt.rcParams.update({'font.size': font})
    
    ax.set_title(title, pad=20)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] 
                 + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(font)

    ax.minorticks_on()
    ax.tick_params(direction='in', length=length1, width=tickwidth,
                which = 'major', top=True, right=True, pad=7)
    ax.tick_params(direction='in', length=length2, width=tickwidth,
                which = 'minor', top=True, right=True)
