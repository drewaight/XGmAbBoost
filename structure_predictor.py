import os
import sys
import json
import re
import pickle
import natural_sort
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, ion, draw, show
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import seaborn as sns

def main(cwd, total_sc, cutoff, xgb_clf_file, xgb_reg_file):

    xgb_clf_file='/media/hdd1/roproQ3drew/drewdock_exec/final_xgb_class.dat'
    xgb_reg_file='/media/hdd1/roproQ3drew/drewdock_exec/final_xgb_regress.dat'

    df = total_sc.copy()
    df_total = pd.DataFrame()
    df.reset_index(drop=True, inplace=True)
    df2 = df.copy()
    df.drop(['model_name', 'name','H1_cluster', 'complex', 'id', 'index',
            'H1_distance', 'H1_sequence', 'H2_cluster', 'H2_distance', 'nres_all', 'ref',
            'H2_sequence', 'H3_cluster', 'H3_distance', 'H3_sequence', 'L1_cluster',
            'L1_distance', 'L1_sequence', 'L2_cluster', 'L2_distance',
            'L2_sequence', 'L3_cluster', 'L3_distance', 'L3_sequence'], axis=1, inplace=True)

    df_total = df_total.append(df, ignore_index=True)
    df_total.dropna(inplace=True)  
    print(len(df_total.columns))

    plt.figure(figsize=(12, 6), facecolor="white")
    sns.distplot(df_total['dG_separated'])
    draw()
    plt.savefig(os.path.join(cwd, 'dG_hist.png'))

    with open(xgb_clf_file, 'rb') as pickle_file:
        xgb_clf = pickle.load(pickle_file)

    df_predict = xgb_clf.predict_proba(df_total)[:,1]
    df_score = pd.DataFrame(df_predict, columns=['score'])
    df_results = df2[['model_name', 'name', 'complex']]
    df_results = df_results.join(df_score)
    df_results.sort_values(by=['score'], inplace=True, ascending=False)

    plt.figure(figsize=(20, 8), facecolor="white")
    plt.xticks(rotation=45, fontsize =12)
    sns.set(font_scale=2)
    sns.set_style('whitegrid')
    ax = sns.scatterplot(x="complex", y="score",
                        hue="model_name", data=df_results,  legend="brief",
                        s=200, alpha=1)
    xmin, xmax = ax.get_xlim()
    ax.set(xlim=(xmin+5,xmax-5))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.grid(True)
    # ax.autoscale(enable=True)
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(markerscale=3, bbox_to_anchor=(1.005, 1), loc=2, 
               borderaxespad=0., handles=handles[1:], labels=labels[1:])
    text = ax.text(0,1.01, "XGBoost Classifier", transform=ax.transAxes)
    plt.legend(markerscale=3, bbox_to_anchor=(1.005, 1), loc=2, 
               borderaxespad=0., handles=handles[1:], labels=labels[1:])

    draw()
    plt.savefig(os.path.join(cwd, 'Classifier_proba.png'), bbox_extra_artists=(legend,text), bbox_inches='tight')

    df_docked = df_results[df_results.score > cutoff]
    if df_docked.empty:
        print("No suitable docking poses found. Sorry")
        return(df_docked)
   
    else:
        xgb_reg = pickle.load(open(xgb_reg_file, "rb"))
        dftd = df_total.loc[df_docked.complex]
        rmsd_rank = xgb_reg.predict(dftd)

        rr = pd.DataFrame(rmsd_rank, columns=['rmsd'])

        df_docked.reset_index(drop=True, inplace=True)
        df_docked = df_docked.join(rr)
        df_docked.sort_values(by='rmsd', inplace=True)
      
        df_docked['name'] = df_docked['name'].str.extract(r"(^complex_(\d{3}|\d{2}|\d{1}))")
        df_docked.to_json(os.path.join(cwd, 'results.json'))
        df_docked.to_csv(os.path.join(cwd, 'results.csv'), sep ='\t')

        plt.figure(figsize=(12, 8), facecolor="white")
        plt.xticks(rotation=45, fontsize =18)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.2)
        sns.set(font_scale=2)
        sns.set_style('whitegrid')
        ax = sns.barplot(x="name", y="rmsd", data=df_docked)
        ax.set_xlabel('')
        ax.set_ylabel('rmsd score')
        draw()
        plt.savefig(os.path.join(cwd, 'Regressor_rank.png'))
        show()
        return(df_docked)

if __name__ == '__main__':
    cwd = '/media/hdd1/XGmAbBoost/dara_gromacs/complex14_LH_complex14_A'
    df = pd.read_json(os.path.join(cwd,'complex14_LH_complex14_A_TOTAL.json')) 
    xgb_clf_file='/media/hdd1/roproQ3drew/drewdock_exec/final_xgb_class.dat'
    xgb_reg_file='/media/hdd1/roproQ3drew/drewdock_exec/final_xgb_regress.dat'
    main(cwd=cwd, total_sc=df, cutoff=0.01, xgb_clf_file=xgb_clf_file, xgb_reg_file=xgb_reg_file)