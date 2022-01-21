"""
Author: Andy Tzanidakis (atzanida@caltech.edu//atzanida@uw.edu)
Edits: Tawny Sit (tsit@caltech.edu)

Generates interactive figures for CLU scanning web app.
"""

import plotly
import plotly.express as px
from astropy.io import ascii
from datetime import datetime, timedelta
from calendar import monthrange
from astropy.time import Time
import numpy as np

# Query date today
today = datetime.today()

# Fetch data from that date!
data = ascii.read(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii")

def data_prep(data, key_id='peak_abs_mag'):
    X1, Y1, N1 = [], [], []
    X2, Y2, N2 = [], [], []
    for yy in [2021,2022]:
        for ij in range(1, 13):
            max_days = monthrange(yy, ij)[1]
            T1 = Time(f"{yy}-{ij}-01")
            T2 = Time(f"{yy}-{ij}-{max_days}")

            w_n_clf = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2)& (data['classification'].data.data=="0"))
            w_clf = np.where((data['saved_date']>=T1) & (data['saved_date']<=T2) & (data['classification'].data.data!="0"))

            if len(w_n_clf[0])>0:
                for i in range(len(w_n_clf[0])):
                    #rand = np.random.normal(0, 0.09)
                    #X1.append(ij+rand)
                    X1.append(data[w_n_clf]['saved_date'][i])
                    Y1.append(data[w_n_clf][key_id][i])
                    N1.append(data[w_n_clf]['ZTF_id'][i])
                for i in range(len(w_clf[0])):
                    #rand = np.random.normal(0, 0.09)
                    #X2.append(ij+rand)
                    X2.append(data[w_clf]['saved_date'][i])
                    Y2.append(data[w_clf][key_id][i])
                    N2.append(data[w_clf]['ZTF_id'][i])
            else:
                continue
    X1, Y1, N1 = np.array(X1), np.array(Y1), np.array(N1)
    X2, Y2, N2 = np.array(X2), np.array(Y2), np.array(N2)

    return (X1, Y1, N1, X2, Y2, N2)


def galaxy_sep_int(data):
    """Returns the interactive pyplot figure of soruces per month for galaxy seperation"""
    X1, Y1, N1, X2, Y2, N2 = data_prep(data, key_id='clu_d_kpc')

    clf_tag = ["Not Classified" for _ in range(len(X1))]
    n_clf_tag = ["Classified" for _ in range(len(X2))]

    X0, Y0 = np.concatenate([X1, X2]), np.concatenate([Y1, Y2])
    N0 = np.concatenate([N1, N2])
    Nt = np.concatenate([clf_tag, n_clf_tag])

    fig = px.scatter(x=X0, y=Y0, color=Nt, range_y=(-3, 100), hover_name=N0, color_discrete_map={"Not Classified":"#ef553b", "Classified":"#00cc96"},
                    labels={
                        "x":r"Date Saved",
                        "y":r"CLU Sep. [kpc]",
                        "color":r"Classification Label"
                    })

    j = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

    return (j)

def peak_app_mag_int(data):
    """Returns the interactive pyplot figure of soruces per month for galaxy seperation"""
    X1, Y1, N1, X2, Y2, N2 = data_prep(data, key_id='peak_app_mag')
    clf_tag = ["Not Classified" for _ in range(len(X1))]
    n_clf_tag = ["Classified" for _ in range(len(X2))]

    X0, Y0 = np.concatenate([X1, X2]), np.concatenate([Y1, Y2])
    N0 = np.concatenate([N1, N2])
    Nt = np.concatenate([clf_tag, n_clf_tag])

    fig = px.scatter(x=X0, y=Y0, color=Nt, range_y=(22, 13), hover_name=N0, color_discrete_map={"Not Classified":"#ef553b", "Classified":"#00cc96"},
                    labels={
                        "x":r"Date Saved",
                        "y":r"Peak Apparent Magnitude",
                        "color":r"Classification Label"
                    })

    j = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

    return (j)

def peak_abs_mag_int(data):
    X1, Y1, N1, X2, Y2, N2 = data_prep(data, key_id='peak_abs_mag')
    clf_tag = ["Not Classified" for _ in range(len(X1))]
    n_clf_tag = ["Classified" for _ in range(len(X2))]

    X0, Y0 = np.concatenate([X1, X2]), np.concatenate([Y1, Y2])
    N0 = np.concatenate([N1, N2])
    Nt = np.concatenate([clf_tag, n_clf_tag])

    fig = px.scatter(x=X0, y=Y0, color=Nt, range_y=(-6, -24), hover_name=N0, color_discrete_map={"Not Classified":"#ef553b", "Classified":"#00cc96"},
                    labels={
                        "x":r"Date Saved",
                        "y":r"Peak Absolute Magnitude",
                        "color":r"Classification Label"
                    })

    j = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

    return (j)
