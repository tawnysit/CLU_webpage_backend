from flask import Flask, render_template, request, redirect
from flask_sqlalchemy import SQLAlchemy as sqa
from datetime import datetime, timedelta
from astropy.io import ascii
from astropy.time import Time
from flask import Markup
from calendar import monthrange
import numpy  as np
import os
from natsort import natsorted
import webbrowser
from threading import Timer

from cleanup import del_old_files
import query_clu_catalog as clu_cat
import generate_figures as gf
import generate_tables as gt

import argparse

global CLU_dat

# Query date today
today = datetime.today()
app = Flask(__name__) # initite flask app
app.config['UPLOAD_FOLDER'] = os.path.join('static', 'figstables') # flask app configuration - ideally move into a config.py for deployment

# command-line arguments for initiating the website - probably have to figure out how to automate these settings for eventual deployment...
parser = argparse.ArgumentParser(description="ZTF-II CLU Scanning and Statistics")
parser.add_argument('-d', '--date', type=str, help='add date from which to query argument -d YYYY-MM-DD or -d all; set to 2021-01-01 if no argument entered')
parser.add_argument('-r', '--refresh', type=str, help='add force refresh from date -r YYYY-MM-DD or -r all; does not force refresh by default')
args = parser.parse_args()

files = natsorted(os.listdir('query_data'))
if len(files)>0:
    file = files[-1]
    print (f"Most recent query: {file}")
else:
    file = None
    print("No recent queries.")

if file:
    query_key = False
    if int(file.split("_")[2])==today.year and int(file.split("_")[3])==today.month and int(file.split("_")[4].split(".")[0])==today.day:
        query_key = True

        print("Your CLU catalog is up to date!")
        
        # Query your CLU catalog or redownload data if force refresh & set as global parameter
        if args.refresh:
            print(f"Force refresh: {args.refresh}")
            
            if args.refresh=='all':
                if args.date:
                    _ = clu_cat.QUERYME(date_query=args.date) # force redownload entire catalog from date
                else:
                    _ = clu_cat.QUERYME(date_query='2021-01-01') # if no date provided, download from beginning of 2021
            
            else:
                if args.date:
                    _ = clu_cat.QUERYME(date_query=args.date, existing_query_file=file, force_refresh=args.refresh)
                else:
                    _ = clu_cat.QUERYME(date_query='2021-01-01', existing_query_file=file, force_refresh=args.refresh)
                    
            CLU_dat = ascii.read(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii")
            gf.main(CLU_dat) # generate updated figures
            gt.table1(CLU_dat)
        
        else:
            CLU_dat = ascii.read(f"query_data/{file}")

        from interactive_figs import *

    if query_key==False:
        print ("Looks like you don't have the lastest version the CLU catalog. Downloading now....")
        
        if args.refresh=='all':
            if args.date:
                _ = clu_cat.QUERYME(date_query=args.date)
            else:
                _ = clu_cat.QUERYME(date_query='2021-01-01')
                
        else:
            if args.date:
                _ = clu_cat.QUERYME(date_query=args.date, existing_query_file=file, force_refresh=args.refresh)
            else:
                _ = clu_cat.QUERYME(date_query='2021-01-01', existing_query_file=file, force_refresh=args.refresh)
        
        from interactive_figs import *
        CLU_dat = ascii.read(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii")
        gf.main(CLU_dat)
        gt.table1(CLU_dat)

# download data if no previous queries
else:
    print("Downloading latest CLU catalog...")
    
    if args.date:
        _ = clu_cat.QUERYME(date_query=args.date)
    else:
        _ = clu_cat.QUERYME(date_query='2021-01-01')
    
    from interactive_figs import *
    CLU_dat = ascii.read(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii")
    gf.main(CLU_dat)
    gt.table1(CLU_dat)

# clean up query and static data files older than 2 weeks
del_old_files('query_data', save_days=14)
del_old_files(app.config['UPLOAD_FOLDER'], save_days=14)

def open_browser():
      webbrowser.open_new('http://127.0.0.1:4000/')

@app.route("/")
def index():
    
    full_filename_1 = os.path.join(app.config['UPLOAD_FOLDER'], f'Ns_{today.year}-{today.month}-{today.day}.jpeg')
    if os.path.isfile(full_filename_1)==False: # if summary figure doesn't exist yet, render figure
        gf.figure_1(data)
        full_filename_1 = os.path.join(app.config['UPLOAD_FOLDER'], f'Ns_{today.year}-{today.month}-{today.day}.jpeg')
        
    return render_template('index.html', user_image=full_filename_1)

@app.route("/sources")
def assignment():
    # Sort all Id's by the date
    saved_date = Time(CLU_dat['saved_date'], format='isot')
    theta = CLU_dat[np.argsort(saved_date)]
    
    return render_template('posts.html', posts=theta[::-1])

# List of months of the year
month_lst = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
              'August', 'September', 'October', 'November', 'December']

@app.route("/sources/<any({}):segment>".format(str(month_lst)[1:-1]))
def monthy(segment):
    """ Fetch and push all transients saved by the month argument"""
    months_of_cand = []
    for T in CLU_dat['saved_date']:
        mm = int(T.split('-')[1])
        months_of_cand.append(mm)

    months_of_cand = np.array(months_of_cand)

    for m0 in enumerate(month_lst):
        ii, month_name = m0[0]+1, m0[1]
        if month_name==segment:
            CLU_at_that_month = CLU_dat[np.where(months_of_cand==ii)]

            # Sort them by the date
            CLU_at_that_month = CLU_at_that_month[np.argsort(CLU_at_that_month['saved_date'])]

            #return (render_template(f'{segment}.html', posts=CLU_at_that_month[::-1]))
            return (render_template('posts.html', month=segment, posts=CLU_at_that_month[::-1])) 


@app.route("/sources/unclassified")
def assignment_focus():
    # Sort all Id's by the date
    saved_date = Time(CLU_dat['saved_date'], format='isot')
    theta = CLU_dat[np.argsort(saved_date)] # sort by most recently saved
    
    w = np.where((theta['classification'].data.data=='0') & (theta['in_clu'].data=='True') & (theta['peak_app_mag'].data<21.5)) # get unclassified candidates within CLU cut
    theta_2 = theta[w]

    return render_template('notclf.html', posts=theta_2[::-1])

@app.route("/completeness")
def completeness_stats():
    
    data = ascii.read(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii")
    
    img_filename = os.path.join(app.config['UPLOAD_FOLDER'], f'CC_fig_with_magcuts_{today.year}-{today.month}-{today.day}.jpeg')
    if os.path.isfile(img_filename)==False: # if CC figure doesn't exist yet, render figure (for some reason plt.savefig breaks sometimes if I try to re-render every time...?)
        gf.figure_3(data)
        img_filename = os.path.join(app.config['UPLOAD_FOLDER'], f'CC_fig_with_magcuts_{today.year}-{today.month}-{today.day}.jpeg')
    
    try:
        table_filename = os.path.join(app.config['UPLOAD_FOLDER'], f'CC_{today.year}-{today.month}-{today.day}.ascii')
        table_dat = ascii.read(table_filename, header_start=0)
    except:
        gt.table1(data)
        table_filename = os.path.join(app.config['UPLOAD_FOLDER'], f'CC_{today.year}-{today.month}-{today.day}.ascii')
        table_dat = ascii.read(table_filename, header_start=0)
    
    return render_template('completeness.html', posts=table_dat, user_image=img_filename,
      int_fig_1=Markup(galaxy_sep_int(data)),
      int_fig_2=Markup(peak_app_mag_int(data)),
      int_fig_3=Markup(peak_abs_mag_int(data)))

@app.route("/search")
def search():
    search_query = request.args['search']
    results_inds = [i for i, name in enumerate(CLU_dat['ZTF_id']) if search_query.casefold() in name.casefold()]
    
    if len(results_inds)>0:
        search_results = CLU_dat[results_inds]
        return render_template('posts.html', query=search_query, posts=search_results)
    
    else: # return error page if search returned no results
        return render_template('noresult.html', query=search_query)

if __name__=='__main__':
    Timer(1, open_browser).start()
    #app.run(port=4000, debug=True)
    app.run(port=4000, debug=True, use_reloader=False)
