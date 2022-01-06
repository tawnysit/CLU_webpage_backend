"""
Author: Andy Tzanidakis (atzanida@caltech.edu//atzanida@uw.edu)
Edits: Tawny Sit (tsit@caltech.edu)

Refined methods & for Generating Figures for CLU.
Diagnostic Figures & Survey Performance
"""

from astropy.io import ascii
from astropy.time import Time
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from datetime import date, timedelta
from tqdm import tqdm
from calendar import monthrange
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import warnings
warnings.filterwarnings('ignore')

global path, today

# Date today; If you want your website to display older day data, adjust timedelta (by days)
today = date.today() - timedelta(days=0)

# Set save directory as the global path
#path = "/Users/andytzanidakis/Desktop/astro_research/fritz_tools/scripts/CLU_filters/CLU_scanning/clu_scanning_page_backend/query_data/"
path = "C:/Users/tawny/OneDrive/Documents/ZTF/clu_scanning_page_backend/clu_scanning_page_backend_beta/"

def figure_1(data, years=['2021', '2022']):
    """ Will generate & store figure for:

    Number of saved sources to CLU per month:
    1) All saved sources (for all sources in CLU volume (i.e z<0.05))
    2) All souces that are: M>=-17 Abs. mag
    3) All sources that are: M<-17 Abs. mag
    4) & 5) All sources that meet 2&3 and are d_CLU<=150 Mpc & d_CLU<30kpc

    Input:
    data (astropy.Table): Astropy.Table that contains all important CLU information
    year (str): Year you want to conduct your query
    """
    # Select all transients that have CLU redshift <0.05
    data = data[data['z']<0.05]
    
    fig, ax = plt.subplots(len(years), sharex=True, figsize=(15,0.5+(5*len(years))))
    for ind, year in enumerate(years):
        for ij in range(1, 13):
            max_days = monthrange(int(year), ij)[1] # maximum days of each month
            T1, T2 = Time(f"{year}-{ij}-01"), Time(f"{year}-{ij}-{max_days}") # start and end of Time query

            # Query Data:

            w = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2)) # all saved TRANSIENTS (within z<0.05)
            N_total = len(w[0]) # total number of saved souces

            # Luminosity cut
            w_special_sublum = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2) & (data['luminous_event']=="False")) # all transients Mabs>=-17
            w_special_lum = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2) & (data['luminous_event']=="True")) # all saved transients Mabs<-17
            N_total_lum = len(w_special_lum[0]) # total number of luminous events
            N_total_sublum = len(w_special_sublum[0]) # total number of subluminous events

            # Luminosity & distance cut:
            w_special_sublum_dcut = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2) & (data['luminous_event']=="False") & (data['z']<=0.0331) & (data['clu_d_kpc']<=30) & (data['clu_d_host']<100)) # Mabs>=-17 & distance<=150Mpc & <=30kpc & within 100"
            w_special_lum_dcut = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2) & (data['luminous_event']=="True") & (data['z']<=0.0331) & (data['clu_d_kpc']<=30) & (data['clu_d_host']<100)) # Mabs<-17 & distance<=150Mpc & <=30kpc & within 100"
            N_total_lum_dcut = len(w_special_lum_dcut[0]) # total number of luminous & d<=150 Mpc
            N_total_sublum_dcut = len(w_special_sublum_dcut[0]) # total number of subluminous & d<=150 Mpc

            if N_total>0:
                ax[ind].bar(ij, N_total, width=0.2, color='#166666') # total
                ax[ind].text(ij, N_total+3, f"{N_total}", color='#166666') # total text

                ax[ind].bar(ij+0.2, N_total_lum, alpha=1, width=0.2, color='limegreen') # lumimnous
                ax[ind].text(ij+0.2, N_total_lum+5,f"{N_total_lum}", color='limegreen') # luminous text
                ax[ind].bar(ij+0.2, N_total_lum_dcut, alpha=1, width=0.2, color='forestgreen') # luminous and <=150Mpc
                ax[ind].text(ij+0.15, N_total_lum_dcut+3,f"({N_total_lum_dcut})", color='forestgreen') # luminous and <=150Mpc text

                ax[ind].bar(ij+0.4, N_total_sublum, width=0.2, color='salmon') # subluminous
                ax[ind].text(ij+0.4, N_total_sublum+5, f"{N_total_sublum}", color='salmon') # subluminous text
                ax[ind].bar(ij+0.4, N_total_sublum_dcut, width=0.2, color='Crimson') # subluminous & <=150Mpc
                ax[ind].text(ij+0.45, N_total_sublum_dcut+3, f"({N_total_sublum_dcut})", color='Crimson') # subluminous & <=150Mpc text

        ax[ind].set_xlim(0.5, 13)
        ax[ind].set_ylim(-.01, 160)
        if ind==0:
            ax[ind].bar(-1, 0, color='#166666', label='All sources')
            ax[ind].bar(-1, 0, color='salmon', label='M$\geq$-17')
            ax[ind].bar(-1, 0, color='limegreen', label='M$<$-17')
            ax[ind].bar(-1, 0, color='Crimson', label='M$\geq$-17, z$\leq$0.0331, $\delta_{CLU}$ $\leq$30kpc')
            ax[ind].bar(-1, 0, color='forestgreen', label='M$<$-17, z$\leq$0.0331, $\delta_{CLU}$ $\leq$30kpc')
            ax[ind].legend()
        if ind==len(years)-1:
            ax[ind].set_xlabel(r"Month", fontsize=16)
        ax[ind].set_ylabel(r"Number of Sources Saved", fontsize=16)
        ax[ind].set_title(year, fontsize=18)
        
    fig.suptitle(f"ZTF CLU Summary as of {today.year}-{today.month}-{today.day}", fontsize=22, y=0.95)
    plt.savefig(f"static/figstables/Ns_{today.year}-{today.month}-{today.day}.jpeg", format='jpeg', dpi=100, bbox_inches='tight')
    #plt.close()

def clu_cc(d0, cond, mode='app', bw=0.5):
    """Takes the data and a condition (numpy.where condition) and returns the cumulative completeness (CC).

    Input:
    -----
    d0 (astropy.table) original data
    cond (np.array): condition you want to query

    Supporting modes: app (apparent mag) OR abs (absolute magnitude)
    """

    if mode=='app':
        mag_bin = np.arange(14, 22, step=bw)
    elif mode=='abs':
        mag_bin = np.arange(-25, -8, step=bw)
    else:
        logging.warning('Please add your mode type (i.e app or abs)')

    # New data structure
    D = d0[cond]

    X, Y = [], []
    for i in range(len(mag_bin)):
        if mode=='app':
            w = np.where((D['peak_app_mag']<=mag_bin[i]))
            DD = D[w]
            if len(DD)>0:
                N_all = len(DD)
                N_clf = len(np.where(DD['classification'].data.data!='0')[0])
                f0 = N_clf/N_all
                X.append(mag_bin[i])
                Y.append(f0)
            else:
                continue

        elif mode=='abs':
            w = np.where((D['peak_abs_mag']<=mag_bin[i]))
            DD = D[w]
            if len(DD)>0:
                N_all = len(DD)
                N_clf = len(np.where(DD['classification'].data.data!='0')[0])
                f0 = N_clf/N_all

                X.append(mag_bin[i])
                Y.append(f0)
            else:
                continue

    return np.array(X), np.array(Y)

def figure_2(data):
    """ Will generate & store figure for:

    Cumulative completeness for bins in apparent magnitude. Default is 0.5 mag bins in apparent magnitude.

    Input:
    data (astropy.Table): Astropy.Table that contains all important CLU information """

    # Query data: Conditions
    sl_cq_1 = np.where((data['z']<0.05) & (data['luminous_event']=='False') & (data['clu_d_kpc']<=30) & (data['clu_d_host']<100)) # z<0.05(200 Mpc), 30Kpc/100" seperation, and is subluminous
    sl_cq_2 = np.where((data['z']<=0.0331) & (data['luminous_event']=='False') & (data['clu_d_kpc']<=30) & (data['clu_d_host']<100)) # z<0.0331 (150 Mpc), 30Kpc/100" seperation, and is subluminous

    l_cq_1 = np.where((data['z']<0.05) & (data['luminous_event']=='True') & (data['clu_d_kpc']<=30) & (data['clu_d_host']<100)) # z<0.05(200 Mpc), 30Kpc/100" seperation, and is luminous
    l_cq_2 = np.where((data['z']<=0.0331) & (data['luminous_event']=='True') & (data['clu_d_kpc']<=30) & (data['clu_d_host']<100)) # z<0.0331 (150 Mpc), 30Kpc/100" seperation, and is luminous

    # Query data: CC
    sl_cc_1 = clu_cc(data, sl_cq_1) # subluminous CC 200 Mpc
    sl_cc_2 = clu_cc(data, sl_cq_2) # subluminous CC 150 Mpc

    l_cc_1 = clu_cc(data, l_cq_1) # luminous CC 200 Mpc
    l_cc_2 = clu_cc(data, l_cq_2) # luminous CC 150 Mpc

    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(10,5))
    fig.subplots_adjust(wspace=0)

    # Luminous CC
    ax[0].plot(l_cc_1[0], l_cc_1[1], color='k', ls='--', label="<200 Mpc, <30 kpc")
    ax[0].plot(l_cc_2[0], l_cc_2[1], color='navy', label="<150 Mpc, <30 kpc")

    # Subluminous CC
    ax[1].plot(sl_cc_1[0], sl_cc_1[1], color='m', ls='--', label="<200Mpc, <30 kpc")
    ax[1].plot(sl_cc_2[0], sl_cc_2[1], color='orange', label="<150 Mpc, <30 kpc")

    for i in range(2):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_xlim(14.5, 21)
        ax[i].set_ylim(0.1, 1.05)
        ax[i].set_xlabel(r"$m_{peak_{ZTF-r}}$", fontsize=20)
        ax[i].axhline(0.5, color='gray', alpha=0.5)
        ax[i].axhline(0.8, color='gray', alpha=0.5)
        ax[i].legend(fontsize=12)

    ax[0].set_ylabel(r"Completeness ($N_{classified} / N$)", fontsize=16)
    ax[0].text(16, 0.78, "80%", alpha=0.5)
    ax[0].text(16, 0.52, "50%", alpha=0.5)
    ax[0].set_title("Luminous", fontsize=16)
    ax[1].set_title("Subluminous", fontsize=16)
    plt.savefig(f"static/figstables/CC_fig_{today.year}-{today.month}-{today.day}.jpeg", format='jpeg', dpi=100, bbox_inches='tight') # sometimes this breaks and I don't know why
    plt.close()
    
def figure_3(data):
    """ Will generate & store figure for:

    Cumulative completeness for bins in apparent magnitude. Default is 0.5 mag bins in apparent magnitude.

    Input:
    data (astropy.Table): Astropy.Table that contains all important CLU information """

    # Query data: Conditions
    sl_cq_1 = np.where((data['z']<0.05) & (data['luminous_event']=='False') & (data['clu_d_kpc']<=30) & (data['peak_abs_mag']<=-10) & (data['peak_app_mag']<=21)) # z<0.05(200 Mpc), 30Kpc seperation, and is subluminous
    sl_cq_2 = np.where((data['z']<=0.0331) & (data['luminous_event']=='False') & (data['clu_d_kpc']<=30) & (data['peak_abs_mag']<=-10) & (data['peak_app_mag']<=21)) # z<0.0331 (150 Mpc), 30Kpc seperation, and is subluminous

    l_cq_1 = np.where((data['z']<0.05) & (data['luminous_event']=='True') & (data['clu_d_kpc']<=30) & (data['peak_abs_mag']<=-10) & (data['peak_app_mag']<=21)) # z<0.05(200 Mpc), 30Kpc seperation, and is subluminous
    l_cq_2 = np.where((data['z']<=0.0331) & (data['luminous_event']=='True') & (data['clu_d_kpc']<=30) & (data['peak_abs_mag']<=-10) & (data['peak_app_mag']<=21)) # z<0.0331 (150 Mpc), 30Kpc seperation, and is subluminous

    # Query data: CC
    sl_cc_1 = clu_cc(data, sl_cq_1) # subluminous CC 200 Mpc
    sl_cc_2 = clu_cc(data, sl_cq_2) # subluminous CC 150 Mpc

    l_cc_1 = clu_cc(data, l_cq_1) # luminous CC 200 Mpc
    l_cc_2 = clu_cc(data, l_cq_2) # luminous CC 150 Mpc

    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(10,5))
    fig.subplots_adjust(wspace=0)

    # Luminous CC
    ax[0].plot(l_cc_1[0], l_cc_1[1], color='k', ls='--', label="<200 Mpc, <30 kpc")
    ax[0].plot(l_cc_2[0], l_cc_2[1], color='navy', label="<150 Mpc, <30 kpc")

    # Subluminous CC
    ax[1].plot(sl_cc_1[0], sl_cc_1[1], color='m', ls='--', label="<200Mpc, <30 kpc")
    ax[1].plot(sl_cc_2[0], sl_cc_2[1], color='orange', label="<150 Mpc, <30 kpc")

    for i in range(2):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_xlim(14.5, 21)
        ax[i].set_ylim(0.1, 1.05)
        ax[i].set_xlabel(r"$m_{peak_{ZTF-r}}$", fontsize=20)
        ax[i].axhline(0.5, color='gray', alpha=0.5)
        ax[i].axhline(0.8, color='gray', alpha=0.5)
        ax[i].legend(fontsize=12)

    ax[0].set_ylabel(r"Completeness ($N_{classified} / N$)", fontsize=16)
    ax[0].text(16, 0.78, "80%", alpha=0.5)
    ax[0].text(16, 0.52, "50%", alpha=0.5)
    ax[0].set_title("Luminous", fontsize=16)
    ax[1].set_title("Subluminous", fontsize=16)
    plt.savefig(f"static/figstables/CC_fig_with_magcuts_{today.year}-{today.month}-{today.day}.jpeg", format='jpeg', dpi=100, bbox_inches='tight')
    plt.close()

def main(data):
    """Generate & render figures"""
    figure_1(data)
    figure_2(data)
    figure_3(data)

if __name__ == "__main__":
    main(data)
