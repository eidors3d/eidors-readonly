import requests
import json
from datetime import datetime
import calendar
from os.path import getmtime
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

releases = [
            "EIDORS 3.0",
            "EIDORS 3.1",
            "EIDORS 3.2",
            "eidors-v3.3",
            "eidors-v3.4",
            "eidors-v3.5",
            "eidors-v3.6",
            "eidors_v3.7",
            "eidors-v3.7.1",
            "eidors-v3.8",
            "eidors-v3.9",
            "eidors-v3.9.1",
            "eidors-v3.10"
]

url = "https://sourceforge.net/projects/eidors3d/files/eidors-v3/{}/stats/json?start_date={}&end_date={}&period=monthly&os_by_country=false"
zipname = 'download-stats.zip'
ziptime = datetime.utcfromtimestamp(getmtime(zipname))


def plot_stats():
    dr = pd.date_range('2005-11-01',datetime.now(),freq='M')
    cntdf = pd.DataFrame(index=dr)
    with zipfile.ZipFile(zipname, 'r') as zipf:
        for filename in zipf.namelist():
            with zipf.open(filename) as file:
                js = json.load(file)
                if not js['countries']:
                    continue
                year = int(filename[-12:-8])
                month = int(filename[-7:-5])
                last_day = calendar.monthrange(year, month)[1]
                data = dict(js['countries'])
                for c in data.keys():
                    if c not in cntdf.columns:
                        cntdf[c] = 0
                    cntdf[c][datetime(year,month,last_day)] += data[c]
                pass

    # idx = cntdf.columns[cntdf.sum(axis=0) >= 250] # countries with 500 downloads
    idx = cntdf.sum(axis=0).sort_values(ascending=False)[0:19].index # top 19 countries
    small = cntdf[idx].copy()
    small['Other'] = cntdf.sum(axis=1) - small.sum(axis=1)
    smallyear = small.groupby(small.index.year).sum()
    frac = smallyear.divide(smallyear.sum(axis=1), axis='index')
    frac.plot.area(colormap='tab20c', stacked=True, linewidth=0,xlim=[2010,2022],ylim=[0,1])
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.gca().legend(handles[::-1], labels[::-1], loc='center left',  bbox_to_anchor=(1.03, 0.5),fontsize=9)
    plt.savefig("countries.jpg", bbox_inches="tight",dpi=150)


    reldf = pd.DataFrame(index=dr)
    with zipfile.ZipFile(zipname, 'r') as zipf:
        names = zipf.namelist()
        for release in releases:
            vec = np.zeros([len(dr)])
            for idx, date in enumerate(dr):
                fname = "{}-{}-{:02d}.json".format(release, date.year, date.month)
                if fname not in names:
                    continue
                with zipf.open(fname) as file:
                    js = json.load(file)
                    vec[idx] = js['total']
            reldf[release] = vec
    reldf.plot.area(colormap='cividis',stacked=True,linewidth=0,ylim=[0,500])
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.gca().legend(handles[::-1], labels[::-1])
    plt.savefig('releases.jpg', bbox_inches="tight",dpi=150)
    pass



def get_stats():
    for release in releases:
        for year in range(2004, datetime.now().year + 1):
            for month in range(1, 13):
                last_day = calendar.monthrange(year, month)[1]
                dt = datetime(year, month, last_day)
                if dt < ziptime or dt >= datetime.now():
                    continue
                fname = "{}-{}-{:02d}.json".format(release, year, month)

                start = "{}-{:02d}-{:02d}".format(year, month, 1)
                end = "{}-{:02d}-{:02d}".format(year, month, last_day)
                r = requests.get(url.format(release, start, end))
                js = r.json()
                if not js['messages'] or "no download statistics" not in js['messages'][0]:
                    with open(fname, 'w') as outfile:
                        json.dump(js, outfile)
                    with zipfile.ZipFile(zipname, 'a') as zipf:
                        zipf.write(fname)


if __name__ == '__main__':
    get_stats()
    plot_stats()




