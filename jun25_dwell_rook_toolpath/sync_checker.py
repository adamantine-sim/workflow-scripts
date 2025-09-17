import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import toolpath_writer
import os
from scipy.optimize import minimize

# XXX Note the hardcoded dwell values in rewrite_toolpath! Also this WILL overwrite an existing sim path

def load(tf_file=os.path.join(os.getcwd(),'tf_log.csv'),inp_file=os.path.join(os.getcwd(), "scan_path_test.inp")):
    df_e = pd.read_csv(tf_file,names=['te','xe','ye','ze','pe'])
    df_s = pd.read_csv(inp_file,names=['ts','xs','ys','zs','ps'])
    return df_e, df_s # experimental and simulated motion dataframes

def sync(df_e, df_s):
    '''function to stich two scans by XYZ location and calculate lag between the two'''
    # Set index as timedelta, but keep the time columns for reference
    df_e = df_e.set_index(pd.to_timedelta(df_e["te"], unit="s"), drop=False)
    df_s = df_s.set_index(pd.to_timedelta(df_s["ts"], unit="s"), drop=False)

    df_e = df_e[df_e.index.notnull()]
    df_s = df_s[df_s.index.notnull()]

    # Remove duplicate indices
    df_e = df_e[~df_e.index.duplicated(keep='first')]
    df_s = df_s[~df_s.index.duplicated(keep='first')]

    # Join by DISTANCE
    search_time = 200 # s
    search_distance = 0.001 # m
    match_distances = []
    ts_matches = []
    xs_matches = []
    ys_matches = []
    zs_matches = []
    ps_matches = []
    lags = []
    for index, row in tqdm(df_e.iterrows(),total=len(df_e)):
        mask = (df_s['ts'] > (row['te'] - search_time)) & (df_s['ts'] <= (row['te']+search_time))
        masked = df_s.loc[mask]
        if len(masked) == 0:
            ts_matches.append(np.nan)
            xs_matches.append(np.nan)
            ys_matches.append(np.nan)
            zs_matches.append(np.nan)
            ps_matches.append(np.nan)
            lags.append(np.nan)
            continue
        ds = np.sqrt((masked['xs']-row['xe'])**2+(masked['ys']-row['ye'])**2+(masked['zs']-row['ze'])**2)

        if ds.min() > search_distance:
            # try again with the retraction offset XXX
            ds = np.sqrt((masked['xs']-row['xe'])**2+(masked['ys']-row['ye'])**2+(masked['zs']-row['ze'])**2)
            if ds.min() > search_distance:
                ts_matches.append(np.nan)
                xs_matches.append(np.nan)
                ys_matches.append(np.nan)
                zs_matches.append(np.nan)
                ps_matches.append(np.nan)
                lags.append(np.nan)
                continue

        match_distances.append(ds.min())
        match_idx = ds.idxmin()

        ts_matches.append(masked.loc[match_idx,'ts'])
        xs_matches.append(masked.loc[match_idx,'xs'])
        ys_matches.append(masked.loc[match_idx,'ys'])
        zs_matches.append(masked.loc[match_idx,'zs'])
        ps_matches.append(masked.loc[match_idx,'ps'])
        lags.append(ts_matches[-1] - row['te'])

    df_e['ts_match'] = ts_matches
    df_e['xs_match'] = xs_matches
    df_e['ys_match'] = ys_matches
    df_e['zs_match'] = zs_matches
    df_e['ps_match'] = ps_matches
    df_e['sim_lag'] = lags
    return df_e

def plot(df_e):
    sns.lineplot(df_e,x='te',y='xe', color='blue', linewidth=2.0)
    sns.lineplot(df_e,x='te',y='xs_match', color='orange', linewidth=1.0)
    ax2 = plt.twinx()
    sns.lineplot(data=df_e, x='te', y='sim_lag', color="g", ax=ax2)
    plt.plot()
    input('press enter to continue')

def rewrite_toolpath(dwell_0_offset,dwell_1_offset,dwell_section_offset):
    cwd = os.getcwd()
    toolpath_info = {
        'print_path'    : os.path.join(cwd, "print_layers"),
        'reheat_path'   : os.path.join(cwd, "reheat_layers"),
        ### U6 / minimum dwells
        # 'dwell_0'       : [5,5,5,5],
        # 'reheat_power'  : [900,900,900,900],
        # 'dwell_1'       : [5],
        ### Preopt / O1
        'dwell_0'       : [111.95093, 46.80592, 98.10942, 48.97092],
        'reheat_power'  : [990.72487, 687.35519, 531.07650, 1340.78436],
        'dwell_1'       : [55.20886],
        'dwell_0_offset': dwell_0_offset,
        'dwell_1_offset': dwell_1_offset,
        'dwell_section_offset': dwell_section_offset,
        'scan_path_out' : "scan_path_test.inp",
        'includes_end_message': False,
        'set_dwell_every_n_layers': True
    }    
    toolpath_writer.write_toolpath(toolpath_info)

def score_sync(x):
    rewrite_toolpath(x[0],x[1],x[2])
    df_e, df_s = load()
    df_e = sync(df_e,df_s)
    df_e['sim_lag_med'] = df_e.rolling(window=7).median()['sim_lag']
    avg_lag = np.sum(df_e['sim_lag_med']) / np.sum(df_e['sim_lag_med'].notna())
    print(f'd0l:{x[0]},d1l:{x[1]},dsl{x[2]};avg_lag:{avg_lag}')
    return np.absolute(avg_lag)

if __name__ == "__main__":
    df_e, df_s = load() # experimental and simulation scans

    # plot initial sync
    df_e = sync(df_e,df_s)
    plot(df_e)

    # optimize sync
    x0 = np.array([0.586,0.7375,-1.94])
    res = minimize(score_sync, x0, method='nelder-mead',options={'xatol': 0.1, 'disp': True})
    print(res)

    # plot best fit
    rewrite_toolpath(res.x[0],res.x[1],res.x[2])
    df_e, df_s = load()
    df_e = sync(df_e,df_s)
    df_e.to_csv('tf_log_distance_joined.csv')
    plot(df_e)

