import numpy as np
from datetime import date, datetime, timedelta
import time
from astropy.time import Time
from astropy.io import ascii
from astropy.table import Table, Column

global today
today = datetime.today()

def clf_counting(d0, cond, mode='app', bw=0.5):
    """Takes the data and a condition (numpy.where condition) and returns the cumulative number of candidates and total number of classified candidates.

    Input:
    -----
    d0 (astropy.table) original data
    cond (np.array): condition you want to query

    Supporting modes: app (apparent mag) OR abs (absolute magnitude)
    """
    
    D = d0[cond]

    if mode=='app':
        mag_bin = np.arange(15, 21+bw, step=bw)
    elif mode=='abs':
        mag_bin = np.arange(-25, -8, step=bw)
    else:
        logging.warning('Please add your mode type (i.e app or abs)')

    X, Y_all, Y_clf = [], [], []
    for i in range(len(mag_bin)):
        X.append(mag_bin[i])
        if mode=='app':
            w = np.where((D['peak_app_mag']<=mag_bin[i]))
            DD = D[w]
            if len(DD)>0:
                N_clf = len(np.where(DD['classification'].data.data!='0')[0])
                Y_all.append(len(DD))
                Y_clf.append(N_clf)
            else:
                Y_all.append(len(DD))
                Y_clf.append(0)

        elif mode=='abs':
            w = np.where((D['peak_abs_mag']<=mag_bin[i]))
            DD = D[w]
            if len(DD)>0:
                N_clf = len(np.where(DD['classification'].data.data!='0')[0])
                Y_all.append(len(DD))
                Y_clf.append(N_clf)
            else:
                Y_all.append(len(DD))
                Y_clf.append(0)

    return X, np.array(Y_all), np.array(Y_clf)

def get_cand_months(data):
    months_of_cand = []
    years_of_cand = []
    
    for T in data['saved_date']:
        mm = int(T.split('-')[1])
        months_of_cand.append(mm)
        yyyy = int(T.split('-')[0])
        years_of_cand.append(yyyy)
        
    return np.array(months_of_cand), np.array(years_of_cand)

def cc_by_month(data, month_num=today.month, year_num=today.year):
    cc = []
    cand_months, cand_years = get_cand_months(data)
    month_list = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    
    for year in np.unique(cand_years):
        if year<year_num:
            month_lst = month_list
        elif year==year_num:
            month_lst = month_list[:month_num]
        for m0 in enumerate(month_lst):
            ii, month_name = m0[0]+1, m0[1]
            CLU_in_month = data[np.where(cand_months==ii)[0]]
            mags, N_all, N_clf = clf_counting(CLU_in_month, np.where((CLU_in_month['in_clu']=='True') & (CLU_in_month['clu_d_kpc']<=30) & (CLU_in_month['peak_abs_mag']<=-10) & (CLU_in_month['peak_app_mag']<=21) & (CLU_in_month['clu_d_host']<100)))

            if ii==1:
                overall_N = N_all.copy()
                overall_clf = N_clf.copy()
            else:
                overall_N = overall_N + N_all
                overall_clf = overall_clf + N_clf

            N_frac = np.round(100*(N_clf/N_all), decimals=1)
            mags = [str(i) for i in mags]
            cc_in_month = [month_name+' '+str(year)]

            for j in range(len(mags)):
                if N_all[j]!=0:
                    cc_in_month.append(f'{N_clf[j]}/{N_all[j]} ({N_frac[j]}%)')
                elif N_all[j]==0:
                    cc_in_month.append(f'{N_clf[j]}/{N_all[j]}')

            cc.append(cc_in_month)
    
    cc_totals = ['Total']
    overall_fracs = np.round(100*(overall_clf/overall_N), decimals=1)
    for j in range(len(mags)):
        cc_totals.append(f'{overall_clf[j]}/{overall_N[j]} ({overall_fracs[j]}%)')
        
    cc.append(cc_totals)
    
    return np.array(cc, dtype=object), mags

def table1(data):
    cc_data, mags = cc_by_month(data)
    
    cols = []
    col_keys = ['Month']
    for m in mags:
        col_keys.append(m)
        
    for i, key in enumerate(col_keys):
        vals = [s for s in cc_data[:,i]]
        col = Column(vals, name=key)
        cols.append(col)
        
    data_table = Table(cols)
    data_table.write(f'static/figstables/CC_{today.year}-{today.month}-{today.day}.ascii', format='ascii')