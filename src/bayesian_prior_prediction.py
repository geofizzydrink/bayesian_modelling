"""
.. note:: 

    **Module Details**
    
    **Bayesian Prior Prediction Module**
    
    This Module contains the collection of operations to perform Bayesian Statistical Modelling of Biomedical Clinical Trial Data.
    
    Created: 6 October 2023
    
    @author: Dr. Matthew B.J. Purss

"""
#===============================================================================
# represents a production grade implementation of `Patent WO2018/136996 <https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2018136996>`_ :cite:t:`Purss_Patent_2018` 
#===============================================================================
           

import argparse
import scipy.integrate as sp_integrate
import scipy.interpolate as sp_interpolate
import scipy.stats as sp_stats
import scipy.special as sp_special
import scipy.optimize as sp_optimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colours
import json
import logging
import os
import sys
import datetime
import subprocess
import traceback

from bioinfokit.analys import get_data, stat
from statistics import NormalDist
from statsmodels.stats.power import TTestIndPower

from pprint import pprint

import colored
from colored import stylize

process_step = colored.fg("chartreuse_2b") + colored.attr("bold") + colored.attr("underlined")
cleanup_step = colored.fg("red_3b") + colored.attr("bold") + colored.attr("underlined")
script_parameters = colored.fg("dark_orange") + colored.attr("bold")
step_separator = colored.fg("dark_orange") + colored.attr("bold")

colour_dict = mpl_colours.get_named_colors_mapping()

 
logging.basicConfig(filename='output_file.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s - %(funcName)s %(lineno)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)
logger = logging.getLogger('Bayesian Model Simulation Logs')  
logfilename = 'output_file.log'


def fit_curve(x_val,
              y_val):
    
    def func1(x, a, b, c):
        return a*x**2+b*x+c
    
    def func2(x, a, b, c, d):
        return a*x**3+b*x**2+c*x + d
    
    def func3(x, a, b, c):
        return a*x**3+b*x**2+c
    
    def func4(x, a, b):
        return a*np.exp(-b*x)
     
    
    #===========================================================================
    # params, covs = sp_optimize.curve_fit(func1, x_val, y_val)
    # print('===========func1================')
    # print("params: ", params) 
    # print("covariance: ", covs) 
    # SE = np.sqrt(np.diag(covs))
    # print("SE: ", SE) 
    #===========================================================================
    
    #===========================================================================
    # params, covs = sp_optimize.curve_fit(func2, x_val, y_val)
    # print('===========func2================')
    # print("params: ", params) 
    # print("covariance: ", covs) 
    # SE = np.sqrt(np.diag(covs))
    # print("SE: ", SE) 
    #===========================================================================
    
    #===========================================================================
    # params, covs = sp_optimize.curve_fit(func3, x_val, y_val)
    # print('===========func3================')
    # print("params: ", params) 
    # print("covariance: ", covs) 
    # SE = np.sqrt(np.diag(covs))
    # print("SE: ", SE) 
    #===========================================================================
    xdata = np.hstack([[0], x_val])
    ydata = np.hstack([[100], y_val])
    print('xdata',xdata)
    print('ydata',ydata)
    params, covs = sp_optimize.curve_fit(func1, xdata, ydata)#, p0=[1,1,100])
    print('===========func4================')
    print("params: ", params) 
    print("covariance: ", covs) 
    SE = np.sqrt(np.diag(covs))
    print("SE: ", SE) 
    
    fit_test = func1(x_val, *params)
    print('fit_test',fit_test)
    
    lower_params = params-SE 
    lower_fit_test = func1(x_val, *lower_params)
    print('lower_fit_test',lower_fit_test)
    
    
    upper_params = params+SE 
    upper_fit_test = func1(x_val, *upper_params)
    print('upper_fit_test',upper_fit_test)
    
    return fit_test, lower_fit_test, upper_fit_test
    

def fit_data_bspline(x_val,
                    y_val,
                    degree=3):
    
    spl = sp_interpolate.BSpline(x_val,
                                 y_val,
                                 degree)
    return spl 


def fit_data_values(x_val,
                    y_val,
                    degree=3):
    xdata = np.hstack([np.linspace(-10,0,num=100), x_val])
    ydata = np.hstack([np.ones(100)*100, y_val])
    fit_output = np.polynomial.Polynomial.fit(x=xdata,
                               y=ydata,
                               deg=degree,
                               full=False)

    return fit_output
    
    

def load_data(model_input_params, model_run_flag=False): 
    
    def get_study_data(datafile,data_label): 
        if os.path.isfile(datafile):        
            study_data = json.load(open(datafile))
            study_data_df = pd.DataFrame()
            for k in study_data.keys():
                print(study_data[k])
                print('================')
                for study in study_data[k]:
                    study_dict = {'author': k,
                                  'study_id': study['study_id'],
                                  'gender': [],
                                  'total_mice_start': [],
                                  'total_mice_end':[],
                                  'study_duration_weeks': [],
                                  'tumours': []
                                  }
                    df_index = 0
                    for j in study.keys():
                        print(j, study[j])
                        if (j == 'male') and (len(study['male']) > 0):
                            study_dict['gender'].append('male')
                            study_dict['total_mice_start'].append(study['male']['total_mice_start'])
                            study_dict['total_mice_end'].append(study['male']['total_mice_end'])
                            study_dict['study_duration_weeks'].append(study['male']['study_duration_weeks'])
                            study_dict['tumours'].append(study['male']['tumours'])
                            df_index += 1
                                                        
                        elif (j == 'female') and (len(study['female']) > 0):    
                            study_dict['gender'].append('female')
                            study_dict['total_mice_start'].append(study['female']['total_mice_start'])
                            study_dict['total_mice_end'].append(study['female']['total_mice_end'])
                            study_dict['study_duration_weeks'].append(study['female']['study_duration_weeks'])
                            study_dict['tumours'].append(study['female']['tumours'])
                            df_index += 1
                        else:        
                            study_dict[j] = study[j]
                    pprint(study_dict)
                    print(study_data_df)
                    print('...')
                    study_data_df = pd.concat([study_data_df,
                                         pd.DataFrame(study_dict,
                                                      index=np.arange(df_index))])
                        
                        
            if model_input_params['mouse_gender'].lower() == 'both':
                if model_input_params['mouse_line'].lower() == 'all':
                    q = ('(total_mice_end != "NR") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                   
                    study_data_df = pd.DataFrame(study_data_df.query(q.format(mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                               'total_mice_end': 'int',
                                                               'study_duration_weeks': np.float64})
                    
                elif model_input_params['mouse_line'].lower() in ['fvb/n', 'cd-1']:
                    q = ('(total_mice_end != "NR") & (mouse_line == "{mouse_line}") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
           
                    study_data_df = pd.DataFrame(study_data_df.query(q.format(mouse_line=model_input_params["mouse_line"],
                                                                                    mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                               'total_mice_end': 'int',
                                                               'study_duration_weeks': np.float64})
            else:
                if model_input_params['mouse_line'].lower() == 'all':
                    q = ('(total_mice_end != "NR") & (gender == "{mouse_gender}") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                  
                    study_data_df = pd.DataFrame(study_data_df.query(q.format(mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    mouse_gender=model_input_params["mouse_gender"],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                               'total_mice_end': 'int',
                                                               'study_duration_weeks': np.float64})
                    
                elif model_input_params['mouse_line'].lower() in ['fvb/n', 'cd-1']:
                    q = ('(total_mice_end != "NR") & (gender == "{mouse_gender}") & (mouse_line == "{mouse_line}") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                   
                    study_data_df = pd.DataFrame(study_data_df.query(q.format(mouse_line=model_input_params["mouse_line"],
                                                                                    mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    mouse_gender=model_input_params["mouse_gender"],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                               'total_mice_end': 'int',
                                                               'study_duration_weeks': np.float64})
                                                                                    
                                                            
            study_data_df['percent_survived'] = (study_data_df['total_mice_end']/study_data_df['total_mice_start'])*100   
            study_data_df.attrs = {'name': data_label}
            
                
        else:
            raise Exception('model run MUST include a valid set of "current study" data!')    
        
        return study_data_df
    
    def get_historic_data(datafile, data_label):  
        if os.path.isfile(datafile):    
            historic_data = json.load(open(datafile))
            historic_data_df = pd.DataFrame()
            for k in historic_data.keys():
                print(historic_data[k])
                print('================')
                for study in historic_data[k]:
                    study_dict = {'author': k,
                                  'study_id': study['study_id'],
                                  'gender': [],
                                  'total_mice_start': [],
                                  'total_mice_end':[],
                                  'study_duration_weeks': [],
                                  'tumours': []
                                  }
                    df_index = 0
                    for j in study.keys():
                        if (j == 'male') and (len(study['male']) > 0):
                            study_dict['gender'].append('male')
                            study_dict['total_mice_start'].append(study['male']['total_mice_start'])
                            study_dict['total_mice_end'].append(study['male']['total_mice_end'])
                            study_dict['study_duration_weeks'].append(study['male']['study_duration_weeks'])
                            study_dict['tumours'].append(study['male']['tumours'])
                            df_index += 1
                                                        
                            
                        elif (j == 'female') and (len(study['female']) > 0):
                            study_dict['gender'].append('female')
                            study_dict['total_mice_start'].append(study['female']['total_mice_start'])
                            study_dict['total_mice_end'].append(study['female']['total_mice_end'])
                            study_dict['study_duration_weeks'].append(study['female']['study_duration_weeks'])
                            study_dict['tumours'].append(study['female']['tumours'])
                            df_index += 1
                        else:        
                            study_dict[j] = study[j]
                        
                    historic_data_df = pd.concat([historic_data_df,
                                         pd.DataFrame(study_dict,
                                                      index=np.arange(df_index))])
            
            #===========================================================================
            # trim the historic_data_df using the model input parameters
            #===========================================================================
            
            if model_input_params['mouse_gender'].lower() == 'both':
                if model_input_params['mouse_line'].lower() == 'all':
                    q = ('(total_mice_end != "NR") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                   
                    historic_data_df = pd.DataFrame(historic_data_df.query(q.format(mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                                   'total_mice_end': 'int',
                                                                   'study_duration_weeks': np.float64})
                    
                elif model_input_params['mouse_line'].lower() in ['fvb/n', 'cd-1']:
                    q = ('(total_mice_end != "NR") & (mouse_line == "{mouse_line}") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                   
                    historic_data_df = pd.DataFrame(historic_data_df.query(q.format(mouse_line=model_input_params["mouse_line"],
                                                                                    mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                                   'total_mice_end': 'int',
                                                                   'study_duration_weeks': np.float64})
                    
            else:
                if model_input_params['mouse_line'].lower() == 'all':
                    q = ('(total_mice_end != "NR") & (gender == "{mouse_gender}") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                   
                    historic_data_df = pd.DataFrame(historic_data_df.query(q.format(mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    mouse_gender=model_input_params["mouse_gender"],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                                   'total_mice_end': 'int',
                                                                   'study_duration_weeks': np.float64})
                                                                                    
                    
                    
                elif model_input_params['mouse_line'].lower() in ['fvb/n', 'cd-1']:
                    q = ('(total_mice_end != "NR") & (gender == "{mouse_gender}") & (mouse_line == "{mouse_line}") & (study_duration_weeks >= {mice_age_min}) & (study_duration_weeks <= {mice_age_max})')
                   
                    historic_data_df = pd.DataFrame(historic_data_df.query(q.format(mouse_line=model_input_params["mouse_line"],
                                                                                    mice_age_min=model_input_params["mouse_age_weeks"][0],
                                                                                    mice_age_max=model_input_params["mouse_age_weeks"][1],
                                                                                    mouse_gender=model_input_params["mouse_gender"],
                                                                                    ))).sort_values(by='study_duration_weeks').astype({'total_mice_start': 'int',
                                                                   'total_mice_end': 'int',
                                                                   'study_duration_weeks': np.float64})
                                                                                    
                    
                                                                                    
                                                            
            
            historic_data_df['percent_survived'] = (historic_data_df['total_mice_end']/historic_data_df['total_mice_start'])*100
            historic_data_df.attrs = {'name': data_label}
              
    
        else:
            historic_data_df = pd.DataFrame()
        
        return historic_data_df
    
    if model_run_flag: 
        study_test_data_df = get_study_data(model_input_params['test_datasets']['study_datafile'],
                                            model_input_params['test_datasets']["current_data_label"])
        historic_test_data_df = get_historic_data(model_input_params['test_datasets']['historic_datafile'],
                                            model_input_params['test_datasets']["historic_data_label"])
                                                
                                                
        study_control_data_df = get_study_data(model_input_params['control_datasets']['study_datafile'],
                                               model_input_params['control_datasets']["current_data_label"])
        historic_control_data_df = get_historic_data(model_input_params['control_datasets']['historic_datafile'],
                                            model_input_params['control_datasets']["historic_data_label"])
        return study_test_data_df, study_control_data_df, historic_test_data_df, historic_control_data_df
    else:
        study_data_df = get_study_data(model_input_params['study_datafile'],
                                       model_input_params['current_data_label'])        
        historic_data_df = get_historic_data(model_input_params['historic_datafile'],
                                             model_input_params['historic_data_label'])
    
    
    
    
        return historic_data_df, study_data_df


def compare_historic_and_current_data(historic_data_df, 
                                      current_data_df,
                                      datafields = ['total_mice_end'],
                                      alpha_levels= [0.01,0.025,0.05,0.075,0.1],
                                      k_values=[3,5,8],
                                      mean_values=[None],
                                      plot_results=True):
    if len(mean_values) == 1 and mean_values[0] is None:
        mean_values = []
        for i in datafields:
            mean_values.append(None)
    test_inputs = []
    for data_to_test in range(len(datafields)):
        #=======================================================================
        # for current_data in [current_data_df,pd.DataFrame()]:
        #=======================================================================
        #=======================================================================
        # for mean_value in mean_values:
        #=======================================================================
       
        for alpha in alpha_levels:
            for k in k_values:          
                test_inputs.append({'model_run': {'historic_data': historic_data_df.attrs.get('name'),
                                                  'current_data': current_data_df.attrs.get('name'),
                                                  'data_field': datafields[data_to_test],
                                                  'mean_value': mean_values[data_to_test],
                                                  'alpha': alpha,
                                                  'k': k
                                                  },
                                    'historic_data': historic_data_df,
                                    'current_data': current_data_df,
                                    'data_to_test': datafields[data_to_test],
                                    'mean_value': mean_values[data_to_test],
                                    'alpha': alpha,
                                    'k': k})
                        
                    
    power_prior_df = pd.DataFrame()                    
    for pp_input in test_inputs:
        power_prior = get_power_prior(historic_df=pp_input['historic_data'], 
                                      current_df= pp_input['current_data'], 
                                      data_to_test = pp_input['data_to_test'], 
                                      mean_value = pp_input['mean_value'], 
                                      alpha = pp_input['alpha'], 
                                      k=pp_input['k'],
                                      model_run=pp_input['model_run'])
        power_prior_df = pd.concat([power_prior_df,
                                    pd.DataFrame(power_prior,
                                                 index= [0])])
    
    
            #===================================================================
            # q = '(data_to_test in {datafield_to_test}) & (p_value > alpha) & (k == {k_val})'
            #===================================================================
            
            #===================================================================
            # q = '(data_to_test in {datafield_to_test})  & (k == {k_val})'
            #===================================================================
            
            #===================================================================
            # q = '(k == {k_val})'
            # 
            # passed_df = pd.DataFrame(power_prior_df.query(q.format(datafield_to_test = [datafields[d]],
            #                                                        k_val=k_val)))
            #===================================================================
    
    power_prior_df['model_run_id'] = np.arange(len(power_prior_df))        
    
    q = '(p_value >= alpha)'
    passed_df = pd.DataFrame(power_prior_df.query(q))
    #===========================================================================
    # print(passed_df.describe().T)
    # print(passed_df.to_markdown())
    #===========================================================================
    
    if plot_results:
        '''
        #===========================================================================
        # create the Figure as a set of subplots
        #===========================================================================
        '''
        mosaic = np.zeros((len(datafields), 3), dtype=object)
        
        for p in range(len(datafields)):
            mosaic[p,0] = datafields[p]+'_hist'
            mosaic[p,1] = datafields[p]+'_pp'
        mosaic[:,2] = 'Kaplan-Meier'
        
        fig, axes = plt.subplot_mosaic(
            mosaic,
        )
        fig.suptitle("Histogram Distributions of Historic Control Data", fontsize=16)
        fig.set_layout_engine('constrained')
        
        
        secondary_axes = {}
        for d in range(len(datafields)):
            axes[datafields[d]+'_hist'].set_xlabel(datafields[d]) 
            axes[datafields[d]+'_pp'].set_xlabel('Power Prior Model Run') 
            secondary_axes[datafields[d]+'_hist'] = axes[datafields[d]+'_hist'].twinx()
            
            if not current_data_df.empty:
                axes[datafields[d]+'_hist'].set_title('Histogram Comparison between\n'+\
                                        current_data_df.attrs.get('name')+' (current data) & \n'+\
                                        historic_data_df.attrs.get('name')+' (historic data)') 
                
                mu = np.linspace(np.min([historic_data_df[datafields[d]].min(),
                                         current_data_df[datafields[d]].min()]), 
                                 np.max([historic_data_df[datafields[d]].max(),
                                         current_data_df[datafields[d]].max()]),  num = 10000)
               
                historic_normal_dist = sp_stats.norm.pdf(mu, historic_data_df[datafields[d]].mean(), historic_data_df[datafields[d]].std())
                current_normal_dist = sp_stats.norm.pdf(mu, current_data_df[datafields[d]].mean(), current_data_df[datafields[d]].std())
                
                
                
                current_data_df[datafields[d]].plot(kind='hist',
                                                ax=axes[datafields[d]+'_hist'], 
                                                label=current_data_df.attrs['name'],
                                                color='g',
                                               legend=True,
                                               alpha=0.5, 
                                               bins=20)
                axes[datafields[d]+'_hist'].axvline(x = current_data_df[datafields[d]].mean(), 
                                color = 'darkgreen', 
                                linestyle='dashed',
                                linewidth=1.0,
                                label = current_data_df.attrs['name']+' - mean')
                axes[datafields[d]+'_hist'].axvline(x = current_data_df[datafields[d]].median(), 
                                color = 'lightgreen', 
                                linestyle='solid',
                                linewidth=2,
                                label = current_data_df.attrs['name']+' - median')       
                secondary_axes[datafields[d]+'_hist'].plot(mu, 
                             current_normal_dist,
                             color = 'darkgreen',
                             linestyle='dashdot',
                             label = current_data_df.attrs['name']+' Normal Dist.') 
                
                
                historic_data_df[datafields[d]].plot(kind='hist',
                                                ax=axes[datafields[d]+'_hist'], 
                                                label=historic_data_df.attrs['name'],
                                               legend=True,
                                               alpha=0.5, 
                                               bins=20)
                axes[datafields[d]+'_hist'].axvline(x = historic_data_df[datafields[d]].mean(), 
                                color = 'royalblue', 
                                linestyle='dashed',
                                linewidth=1.0,
                                label = historic_data_df.attrs['name']+' - mean')
                axes[datafields[d]+'_hist'].axvline(x = historic_data_df[datafields[d]].median(), 
                                color = 'steelblue', 
                                linestyle='solid',
                                linewidth=2,
                                label = historic_data_df.attrs['name']+' - median')
                
                
                secondary_axes[datafields[d]+'_hist'].plot(mu, 
                             historic_normal_dist,
                             color = 'royalblue',
                             linestyle='dashdot',
                             label = historic_data_df.attrs['name']+' Normal Dist.') 
                print('historic_normal_dist',historic_normal_dist)
                print(np.max([historic_normal_dist])+0.05*np.max([historic_normal_dist]))
                if not np.isnan(np.max([historic_normal_dist])+0.05*np.max([historic_normal_dist])):
                    secondary_axes[datafields[d]+'_hist'].set_ylim(0,np.nanmax([historic_normal_dist])+0.05*np.nanmax([historic_normal_dist]))
                
                
                secondary_axes[datafields[d]+'_hist'].fill_between(mu, 
                                             0, 
                                             historic_normal_dist,
                                            where=(historic_normal_dist < current_normal_dist), 
                                             interpolate=True,
                                             facecolor=(0,0,0,.4), 
                                             edgecolor=(0,0,0,.2),
                                             step='mid', 
                                             linewidth=1.0)
                secondary_axes[datafields[d]+'_hist'].fill_between(mu, 
                                             0, 
                                             current_normal_dist,
                                            where=(current_normal_dist < historic_normal_dist), 
                                             facecolor=(0,0,0,.4), 
                                             edgecolor=(0,0,0,.2),
                                             step='mid', 
                                             linewidth=0.0,
                                             interpolate=False)
                
                
            else:
                axes[datafields[d]+'_hist'].set_title('Histogram Comparison between\n'+\
                                    ' (current data) & \n'+\
                                    historic_data_df.attrs.get('name')+' (historic data)') 
            
                mu = np.linspace(np.min([historic_data_df[datafields[d]].min(),
                                         mean_values[d]]), 
                                 np.max([historic_data_df[datafields[d]].max(),
                                         mean_values[d]]),  num = 10000)
               
                historic_normal_dist = sp_stats.norm.pdf(mu, historic_data_df[datafields[d]].mean(), historic_data_df[datafields[d]].std())
                
                axes[datafields[d]+'_hist'].axvline(x = mean_values[d], 
                                color = 'darkgreen', 
                                linestyle='dashed',
                                linewidth=1.0,
                                label = 'Current Data - mean')
                
                
                historic_data_df[datafields[d]].plot(kind='hist',
                                                ax=axes[datafields[d]+'_hist'], 
                                                label=historic_data_df.attrs['name'],
                                               legend=True,
                                               alpha=0.5, 
                                               bins=20)
                axes[datafields[d]+'_hist'].axvline(x = historic_data_df[datafields[d]].mean(), 
                                color = 'royalblue', 
                                linestyle='dashed',
                                linewidth=1.0,
                                label = historic_data_df.attrs['name']+' - mean')
                axes[datafields[d]+'_hist'].axvline(x = historic_data_df[datafields[d]].median(), 
                                color = 'steelblue', 
                                linestyle='solid',
                                linewidth=2,
                                label = historic_data_df.attrs['name']+' - median')
                
                
                secondary_axes[datafields[d]+'_hist'].plot(mu, 
                             historic_normal_dist,
                             color = 'royalblue',
                             linestyle='dashdot',
                             label = historic_data_df.attrs['name']+' Normal Dist.') 
                
                secondary_axes[datafields[d]+'_hist'].set_ylim(0,np.max([historic_normal_dist])+0.05*np.max([historic_normal_dist]))
                
            
            axes[datafields[d]+'_hist'].legend()
            
            for k_val in k_values:
                
                axes[datafields[d]+'_pp'].set_title('Power Prior Analysis') 
            
                
                axes[datafields[d]+'_pp'].plot(passed_df['model_run_id'].values, 
                             passed_df['p_value'].values,
                             marker='o',
                            markevery=1,
                            linestyle='solid',
                            color = colour_dict[list(colour_dict.keys())[k_val**2]],
                             label = 'p_value - k='+str(k_val)) 
                axes[datafields[d]+'_pp'].plot(passed_df['model_run_id'].values, 
                             passed_df['prior_lambda_ak'].values,
                             marker='v',
                            color = colour_dict[list(colour_dict.keys())[k_val**2]],
                            linestyle='dotted',
                            markevery=1,
                             label = 'prior_lambda_ak - k='+str(k_val)) 
                axes[datafields[d]+'_pp'].plot(passed_df['model_run_id'].values, 
                             passed_df['prior_lambda_apk'].values,
                             marker='*',
                            color = colour_dict[list(colour_dict.keys())[k_val**2]],
                            linestyle='dashdot',
                            markevery=1,
                             label = 'prior_lambda_apk - k='+str(k_val)) 
            axes[datafields[d]+'_pp'].set_yscale('log')
            axes[datafields[d]+'_pp'].grid(axis = 'y', which='both')
            if not passed_df.empty:
                axes[datafields[d]+'_pp'].set_ylim(passed_df['prior_lambda_apk'].min()/2,1.05)
            else:            
                axes[datafields[d]+'_pp'].set_ylim(power_prior_df['prior_lambda_apk'].min()/2,1.05)
                
            
            axes[datafields[d]+'_pp'].legend(ncol=2).set_draggable(state=True)
            
        '''
        #===========================================================================
        # create a Kaplan-Meier survival curve plot to get the trend curve for survival vs age
        #===========================================================================
        '''
        axes['Kaplan-Meier'].set_title('Kaplan-Meier survival curves') 
        historic_data_df.plot(kind='scatter',
                      x='study_duration_weeks',
                      y='percent_survived',
                        ax=axes['Kaplan-Meier'], 
                        label= historic_data_df.attrs.get('name')+' (Historic Dataset) - percent_survived',
                       legend=True)
        if not current_data_df.empty:
            current_data_df.plot(kind='scatter',
                          x='study_duration_weeks',
                          y='percent_survived',
                          marker='v',
                            ax=axes['Kaplan-Meier'], 
                            label=current_data_df.attrs.get('name')+' (Current Dataset) - percent_survived',
                            color='r',
                           legend=True)
        
        tmp_df = pd.concat([historic_data_df,
                            current_data_df])
       
        tmp_df.sort_values(by='study_duration_weeks', inplace=True, ignore_index=True)
        y_fit = fit_data_values(tmp_df['study_duration_weeks'].values,
                        tmp_df['percent_survived'].values,
                        3)
        
        axes['Kaplan-Meier'].plot(y_fit.linspace(n=100,domain=[0,tmp_df['study_duration_weeks'].max()])[0], 
                     y_fit.linspace(n=100,domain=[0,tmp_df['study_duration_weeks'].max()])[1],
                    color = 'darkgreen',
                    linestyle='dashdot',
                    label = 'Polyfit: '+str(y_fit))
        axes['Kaplan-Meier'].set_ylim(0,101)
        axes['Kaplan-Meier'].set_xlim(0,tmp_df['study_duration_weeks'].max()+0.1*tmp_df['study_duration_weeks'].max())
        
        axes['Kaplan-Meier'].legend()

    
    return power_prior_df, passed_df
    
    


def main(input_parameters='',
         build_documentation=False,
         compare_historic_data=False,
         model_data=False,
         log_filename=None):
    '''
    .. note:: 

        **Function Details**
        
        Main Function to drive the Bayesian Statistical Analysis (with Power Prior consideration) for
        Medical Clinical trial data.
    
    Parameters
    ----------
    
    input_parameters: JSON (or file input), required        
        JSON representation of Model Input Parameters. This can be either a JSON string or a JSON formatted file.            
        
    build_documentation: Boolean, optional, default=False
        Flag to build documentation for this project.
    
    compare_historic_data: Boolean, optional, default=False
        Flag to perform a Historic vs Current Study data comparison analysis.    
    
    log_filename: str, optional, default="bayesian_model_output" 
        Filename to use for storing logging information for this model run.
        
    Returns
    -------
    
    output_dict: JSON
        JSON compatible dictionary object is returned with details of the `tested` DGGS Zone.
     
       
    ''' 
    #===========================================================================
    # 
    # 
    #     .. note:: 
    #     
    #         **For Historical Data Comparison analyises**
    #         
    #         (used to determine the best Power Prior parameters to use for the main model run)
    #         the JSON structure should be similar to the following:
    #     
    #         .. code-block:: python
    #                 
    #                 {"historic_data": "all_historic_data", 
    #                  "current_data": "study_data", 
    #                  "datafields": ["percent_survived"], 
    #                  "alpha_levels": [0.025], 
    #                  "k_values": [3, 5, 8], 
    #                  "mean_values": [67.5], 
    #                  "mouse_line": "FVB/N", 
    #                  "mouse_age_weeks": [52, 150]}
    #              
    #                          
    #     For Main Bayesian analyises 
    #     the JSON structure should be similar to the following:
    # 
    #     .. code-block:: python
    #             
    #             {"historic_data": <name of historic data to use>, 
    #              "current_data": <name of current study data to use>, 
    #              "datafields": <list of datafields to use for data comparison>, 
    #              "alpha_levels": <list of alpha threshold values to test for statistical significance>, 
    #              "k_values": <list of "k" calibration parameters to scale the "Power Prior" parameter as defined by Shi et. al. 2023>, 
    #              "mean_values": <list of mean study data values to test for each datafield specified -- only used if no "current_data" dataset is providedf (i.e. not enough data to produce a normal distribution)>, 
    #              "mouse_line": <the "mouse line" to use (e.g. FVB/N, CD-1, All, etc...), 
    #              "mouse_age_weeks": <minimum and maximum age of mice to compare - as a list (i.e. [<min_age>, <max_age>])>}
    #              
    #             e.g.
    #             
    #             {"historic_data": "all_historic_data", 
    #              "current_data": "study_data", 
    #              "datafields": ["percent_survived"], 
    #              "alpha_levels": [0.025], 
    #              "k_values": [3, 5, 8], 
    #              "mean_values": [67.5], 
    #              "mouse_line": "FVB/N", 
    #              "mouse_age_weeks": [52, 150]}
    #              
    # 
    #===========================================================================
    
    global logfilename
    timestamp = datetime.datetime.utcnow().strftime(format='%Y%m%d-%H%M%S.%f')
    if log_filename is None:
        logfilename = 'bayesian_model_output.'+timestamp+'.log'
        logging.basicConfig(filename=logfilename,
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s - Function: %(funcName)s Lineno: %(lineno)s %(levelname)s >> %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG,
                            force=True)
    else:
        logfilename = log_filename+'.'+timestamp+'.log'
        logging.basicConfig(filename=logfilename,
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s - Function: %(funcName)s Lineno: %(lineno)s %(levelname)s >> %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG,
                            force=True)
    
    try:
        if build_documentation:
            logger.info('Build Software Documentation') 
            print(stylize('\n>>> Build Software Documentation',process_step))
            
            base_path = os.path.dirname(os.path.abspath(__file__))
            doc_path = os.path.join(os.path.dirname(base_path),'docs')
            print(base_path)
            print(doc_path)
            os.chdir(doc_path)
            cmd = 'make clean'
            os.system(cmd)
            cmd = 'make pdf'
            os.system(cmd)
            
            os.chdir(base_path)
        else:
            if input_parameters != '':
                if os.path.isfile(input_parameters):
                    logger.info('load JSON from file: '+str(input_parameters))
                    with open(input_parameters) as fp:
                        input_parameters_JSON = json.load(fp)
                else:
                    logger.info('parse JSON from string')  
                    input_parameters_JSON = json.loads(input_parameters)
                logger.info('input_parameters JSON Document: \n'+json.dumps(input_parameters_JSON,indent=4))
                
                print(stylize('input_parameters_JSON:\n\n'+ json.dumps(input_parameters_JSON,indent=4),script_parameters))
            
            print(stylize('\n>>> load data...',process_step))
            
            if compare_historic_data:
                if input_parameters != '':
                    if isinstance(input_parameters_JSON, list):
                        for model_run_params in input_parameters_JSON:
                            print(stylize('model_run_params\n\n'+json.dumps(model_run_params, indent=4), script_parameters))
                            
                            historic_data_df, study_data_df = load_data(model_run_params)
                            
                            print('historic_data_df',historic_data_df)
                            print('study_data_df',study_data_df)
                            if 'mean_values' in model_run_params.keys():
                                power_prior_df, passed_df = compare_historic_and_current_data(historic_data_df, 
                                                                  study_data_df,
                                                                  datafields = model_run_params['datafields'],
                                                                  alpha_levels= model_run_params['alpha_levels'],
                                                                  k_values=model_run_params['k_values'],
                                                                  mean_values=model_run_params['mean_values'])
                            else:
                                power_prior_df, passed_df = compare_historic_and_current_data(historic_data_df, 
                                                                  study_data_df,
                                                                  datafields = model_run_params['datafields'],
                                                                  alpha_levels= model_run_params['alpha_levels'],
                                                                  k_values=model_run_params['k_values'])
                    
                            print('power_prior_df')
                            pd.set_option("expand_frame_repr", False)
                            print(power_prior_df.describe().T)
                            print(power_prior_df.to_markdown())
                            print('--------------------')
                            
                            print('passed_df')
                            print(passed_df.describe().T)
                            print(passed_df.to_markdown())
                            print('--------------------')
                            comparison_run_output_file = 'output_comparison_run.out'
                            
                            comparison_run_output_data = {'model_run_params':model_run_params,
                                                         'comparison_run_output': passed_df[['area_under_curve',
                                                                                                'distribution_intercept',
                                                                                                'data_to_test',
                                                                                                'p_value',
                                                                                                'alpha',
                                                                                                'k',
                                                                                                'prior_lambda_ak',
                                                                                                'prior_lambda_apk']].to_dict(orient='records')}                                        
                                                                    
                            with open(comparison_run_output_file,'w') as fd:
                                fd.write(json.dumps(comparison_run_output_data, indent=4))
                                         
                            
                            
            
            elif model_data:    
                rng = np.random.default_rng()            
                if input_parameters != '':
                    if isinstance(input_parameters_JSON, list):
                        for model_run_params in input_parameters_JSON:
                            print(stylize('model_run_params\n\n'+json.dumps(model_run_params, indent=4), script_parameters))
                            
                            study_test_data_df, study_control_data_df, historic_test_data_df, historic_control_data_df = load_data(model_run_params,
                                                                                                                                   model_run_flag=True)
                            
                            
                            
                            active_test_arr = np.zeros(study_test_data_df[model_run_params['data_x_field']].values[0])
                            for i in range(study_test_data_df[model_run_params['data_y_field']].sum()):
                                active_test_arr[i] = 1
                            rng.shuffle(active_test_arr)
                            print(active_test_arr)
                            
                            
                            active_control_arr = np.zeros(study_control_data_df[model_run_params['data_x_field']].values[0])
                            for i in range(study_control_data_df[model_run_params['data_y_field']].sum()):
                                active_control_arr[i] = 1
                            rng.shuffle(active_control_arr)
                            print(active_control_arr)
                            
                            
                            historic_test_arr = []
                            for rec in historic_test_data_df.iterrows():
                                print(rec[1])
                                print(rec[1].to_json())
                                
                                tmp_arr = np.zeros(rec[1][model_run_params['data_x_field']])
                                for i in range(rec[1][model_run_params['data_y_field']]):
                                    tmp_arr[i] = 1
                                rng.shuffle(tmp_arr)
                                print(tmp_arr)
                                historic_test_arr.append(tmp_arr)
                                print(np.mean(tmp_arr))
                            
                            historic_control_arr = []
                            for rec in historic_control_data_df.iterrows():
                                print(rec[1])
                                print(rec[1].to_json())
                                
                                tmp_arr = np.zeros(rec[1][model_run_params['data_x_field']])
                                for i in range(rec[1][model_run_params['data_y_field']]):
                                    tmp_arr[i] = 1
                                rng.shuffle(tmp_arr)
                                print(tmp_arr)
                                historic_control_arr.append(tmp_arr)
                                print(np.mean(tmp_arr))
                            
                            
                            
                            '''   
                            #===================================================
                            # Compute Control Arm Posterior Probability Distribution
                            #===================================================
                            '''
                            print('historic_control_data_df',historic_control_data_df)
                            print('historic_test_data_df',historic_test_data_df)
                            if len(historic_control_data_df) > 0:
                                control_posteriour_dist = find_posterior(input_data=active_control_arr,
                                                                historic_data=historic_control_arr,
                                                                 prior_dist_type=model_run_params['prior_dist_type'],
                                                                power_prior_lambda=model_run_params['power_prior_lambda'],
                                                                multiple_historic_datasets_flag = True
                                                                 )
                                model_output_dict = {'control_posterior_mean': control_posteriour_dist[0],
                                                     'control_posterior_std': control_posteriour_dist[1],
                                                     'control_num_samples': control_posteriour_dist[2]}
                                if 'number_of_regressions' in model_run_params.keys():
                                    if (model_run_params['power_prior_lambda'] > 0.0) and (model_run_params['number_of_regressions'] > 1):
                                    
                                        tmp_mu_0 = []
                                        tmp_sigma_0 = []
                                        tmp_num_0 = []
                                        for a in range(model_run_params['number_of_regressions']):
                                            tmp_posteriour_dist = find_posterior(input_data=active_control_arr,
                                                                        historic_data=historic_control_arr,
                                                                         prior_dist_type=model_run_params['prior_dist_type'],
                                                                        power_prior_lambda=model_run_params['power_prior_lambda'],
                                                                        multiple_historic_datasets_flag = True,
                                                                        plot_data=False
                                                                         )
                                            tmp_mu_0.append(tmp_posteriour_dist[0])                           
                                            tmp_sigma_0.append(tmp_posteriour_dist[1])  
                                            tmp_num_0.append(tmp_posteriour_dist[2]) 
                                        tmp_df = pd.DataFrame({'mu_0': tmp_mu_0,
                                                               'sigma_0': tmp_sigma_0,
                                                               'tmp_num_0': tmp_num_0})
                                        print(tmp_df)
                                        print(tmp_df.describe())
                                        
                                        
                                        control_posteriour_dist = tmp_df['mu_0'].mean(), tmp_df['sigma_0'].mean(), int(tmp_df['tmp_num_0'].mean())
                                        print('average control_posteriour_dist', control_posteriour_dist)
                                        model_output_dict = {'control_all_regressions_posterior_mean': tmp_df['mu_0'].values.tolist(),
                                                             'control_all_regressions_posterior_std': tmp_df['sigma_0'].values.tolist(),
                                                             'control_all_regressions_num_samples': tmp_df['tmp_num_0'].values.tolist(),
                                                             'control_posterior_mean': control_posteriour_dist[0],
                                                             'control_posterior_std': control_posteriour_dist[1],
                                                             'control_num_samples': control_posteriour_dist[2]}
                            else:
                                control_posteriour_dist = find_posterior(input_data=active_control_arr,
                                                                historic_data=historic_control_arr,
                                                                 prior_dist_type=model_run_params['prior_dist_type'],
                                                                power_prior_lambda=model_run_params['power_prior_lambda'],
                                                                multiple_historic_datasets_flag = False
                                                                 )
                                model_output_dict = {'control_posterior_mean': control_posteriour_dist[0],
                                                     'control_posterior_std': control_posteriour_dist[1],
                                                     'control_num_samples': control_posteriour_dist[2]}
                            
                            
                            mu_0 = control_posteriour_dist[0]                            
                            sigma_0 = control_posteriour_dist[1]  
                            num_samples_0 =  control_posteriour_dist[2] 
                            control_distribution = sp_stats.norm(loc = mu_0, scale = sigma_0) 
                            control_pdf = control_distribution.pdf(np.linspace(mu_0-10*sigma_0,mu_0+10*sigma_0,num_samples_0))
                             
                            if np.any(np.isnan(control_pdf)):
                                null_control_data = True 
                               
                            else:
                                null_control_data = False
                                
                            '''   
                            #===================================================
                            # Compute Active Arm Posterior Probability Distribution
                            #===================================================
                            '''
                            if len(historic_test_data_df) > 0:
                                test_posteriour_dist = find_posterior(input_data=active_test_arr,
                                                                historic_data=historic_test_arr,
                                                                 prior_dist_type=model_run_params['prior_dist_type'],
                                                                power_prior_lambda=model_run_params['power_prior_lambda'],
                                                                multiple_historic_datasets_flag = True
                                                                 )
                            
                            else:
                                test_posteriour_dist = find_posterior(input_data=active_test_arr,
                                                                historic_data=historic_test_arr,
                                                                 prior_dist_type=model_run_params['prior_dist_type'],
                                                                power_prior_lambda=model_run_params['power_prior_lambda'],
                                                                multiple_historic_datasets_flag = False
                                                                 )
                            model_output_dict['test_posterior_mean']= test_posteriour_dist[0]
                            model_output_dict['test_posterior_std']= test_posteriour_dist[1]
                            model_output_dict['test_num_samples']= test_posteriour_dist[2]
                            
                            mu_1 = test_posteriour_dist[0]
                            sigma_1 = test_posteriour_dist[1]
                            num_samples_1 =  test_posteriour_dist[2] 
                            test_distribution = sp_stats.norm(loc = mu_1, scale = sigma_1) 
                            test_pdf = test_distribution.pdf(np.linspace(mu_1-10*sigma_1,mu_1+10*sigma_1,num_samples_1))
                            print('num_samples_0',num_samples_0)
                            print('num_samples_1',num_samples_1)
                            print('sigma_0',sigma_0)
                            print('sigma_1',sigma_1)
                            print('mu_0',mu_0)
                            print('mu_1',mu_1)
                            print('mu_1 - mu_0', mu_1 - mu_0)
                            
                            
                            #===================================================
                            # plt.figure(figsize = (10, 8))
                            #===================================================
                            fig, ax1 = plt.subplots()
                            ax1.plot(np.linspace(mu_1-10*sigma_1,mu_1+10*sigma_1,num_samples_1), 
                                     test_pdf, label = 'test_posteriour_dist Distribution')    
                            ax1.axvline(x=mu_0, ls = '--', color = 'k', label = f'mu_0 = {mu_0}')
                            ax1.axvline(x=mu_1, ls = '--', color = 'r', label = f'mu_1 = {mu_1}')
                            if null_control_data:
                                print('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++')
                                print('Assuming a proxy Standard Deviation for the NULL Control Arm based on a sampling uncertanty of 0.5 (i.e. std = sqrt((0.5 - _mu0)^2) = 1.58113883×10⁻¹')
                                sigma_0 = np.sqrt((0.5 - mu_0)**2)
                                control_distribution = sp_stats.norm(loc=mu_0, scale=sigma_0)
                                control_pdf = control_distribution.pdf(np.linspace(mu_0-10*sigma_0,mu_0+10*sigma_0,num_samples_0))
                             
                                
                            ax1.plot(np.linspace(mu_0-10*sigma_0,mu_0+10*sigma_0,num_samples_0),  
                                     control_pdf, label = 'control_posteriour_dist Distribution') 
                            ax1.axis([np.min([mu_1-3*sigma_1,
                                              mu_0-3*sigma_0]),
                                      np.max([mu_1+3*sigma_1,
                                              mu_0+3*sigma_0]),
                                      0,np.max([np.nanmax(test_pdf)+0.05*np.nanmax(test_pdf),
                                               np.nanmax(control_pdf)+0.05*np.nanmax(control_pdf)])])
                             
                            res = stat()
                            analysis = TTestIndPower()
                            if null_control_data:
                                print('-----------------------------------------------------')
                                print('=====================================================')
                                print('Bayesian Regression Results')
                                print('=====================================================')
                                print('null_control_data',null_control_data, '(i.e. insufficient variation in the data to create a proper normal distribution based on the data')
                                
                                test_df = pd.DataFrame({'test_data': test_distribution.rvs(size=num_samples_1, random_state=rng)})
                                test_df['diff'] = test_df['test_data'] - 0.0
                                
                                res.ztest(df=test_df, 
                                         x = 'test_data', 
                                         x_std = sigma_1, 
                                         mu=mu_0,
                                         alpha = 0.05, 
                                         test_type = 1)
                                print('One Sided Z Score Statistics between Test Arm Distribution and the Control Arm mean')                                
                                print(res.summary)
                                print('Z Score = ',res.result[1],'p value = ',res.result[2])
                                if (np.abs(res.result[1]) < 1.96):
                                    print('Test Arm Distribution and Control Arm Mean is Statistically Similar')
                                else:
                                    print('Test Arm Distribution and Control Arm Mean is NOT Statistically Similar')
                                    
                                res.ztest(df=test_df, 
                                                 x = 'diff', 
                                                x_std = np.std(test_df['diff'].values), 
                                                mu=0.0,
                                                 alpha = 0.05, 
                                                 test_type = 1)
                                print('One Sided Z-score Statistics between the Test Arm vs Control Arm Difference Distribution and the Null Hypothesis H0:_mu1 - _mu0 = 0 (i.e. H0:_mu1 = _mu0)')
                                print(res.summary)
                                print('Z Score = ',res.result[1],'p value = ',res.result[2])
                                if (np.abs(res.result[1]) > 1.96) and (res.result[2] < 0.05):
                                    print('Difference between Test Arm and Control Arm is Statistically Significant')
                                else:
                                    print('Difference between Test Arm and Control Arm is NOT Statistically Significant')
                                cohens_d = (mu_1 - mu_0) / (np.sqrt((sigma_1 ** 2 + sigma_0 ** 2) / 2))                                
                                power = analysis.solve_power(cohens_d, power=None, nobs1=num_samples_1, ratio=num_samples_0/num_samples_1, alpha=0.05)
                                print('Statistical Power: %.3f' % power)    
                                
                                print('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++')
                                print('Assuming a proxy Standard Deviation for the NULL Control Arm based on a sampling uncertanty of 0.5 (i.e. std = sqrt((0.5 - _mu0)^2) = 1.58113883×10⁻¹')
                                sigma_0 = np.sqrt(0.5**2)
                                control_distribution = sp_stats.norm(loc=mu_0, scale=sigma_0)
                                
                                if num_samples_0 == num_samples_1:
                                    test_df['control_data'] = control_distribution.rvs(size=num_samples_0, random_state=rng)
                                
                                else:
                                    if num_samples_0 > num_samples_1:
                                        '''
                                        #=======================================
                                        # number of samples in control arm is greater than in the test arm
                                        #=======================================
                                        '''
                                        test_df = pd.DataFrame({'control_data': control_distribution.rvs(size=num_samples_0, random_state=rng)})
                                        rvs_test = test_distribution.rvs(size=num_samples_1, random_state=rng)
                                        tmp_arr = np.ones(num_samples_0)
                                        tmp_arr = np.where(tmp_arr==1, np.nan, tmp_arr)
                                        for i in range(num_samples_0):
                                            if i < num_samples_1:
                                                tmp_arr[i] = rvs_test[i] 
                                      
                                        test_df['test_data'] = tmp_arr
                                    else:
                                        '''
                                        #=======================================
                                        # number of samples in test arm is greater than in the control arm
                                        #=======================================
                                        '''
                                        rvs_control = control_distribution.rvs(size=num_samples_0, random_state=rng)
                                        tmp_arr = np.ones(num_samples_1)
                                        tmp_arr = np.where(tmp_arr==1, np.nan, tmp_arr)
                                        for i in range(num_samples_1):
                                            if i < num_samples_0:
                                                tmp_arr[i] = rvs_control[i] 
                                      
                                        test_df['control_data'] = tmp_arr
                                
                                test_df['diff'] = test_df['test_data'] - test_df['control_data']  
                                
                                res.ztest(df=test_df, 
                                         x = 'test_data', 
                                         y = 'control_data', 
                                        x_std = sigma_1, 
                                        y_std = sigma_0, 
                                         alpha = 0.05, 
                                         test_type = 2)
                                print('res',res.summary)
                                print('Two Sided Z Score Statistics between Test Arm Distribution and the Control Arm Distribution')                                
                                print(res.summary,res.result)
                                print('Z Score = ',res.result[2],'p value = ',res.result[4])
                                if (np.abs(res.result[2]) < 1.96):
                                    print('Test Arm Distribution and Control Arm Distribution are Statistically Similar')
                                    significant_flag = True
                                else:
                                    print('Test Arm Distribution and Control Arm Distribution are NOT Statistically Similar')
                                    significant_flag = False
                                
                                model_output_dict['stats_distribution_comparision'] = {'zscore': res.result[2],
                                                                                       'p_value': res.result[4],
                                                                                       'Test_and_Control_Distributions_Statistically_Similar': significant_flag}
                              
                                res.ztest(df=test_df, 
                                                 x = 'diff', 
                                                x_std = np.nanstd(test_df['diff'].values), 
                                                mu=0.0,
                                                 alpha = 0.05, 
                                                 test_type = 1)
                                print('Z-score Statistics between the Test Arm vs Control Arm Difference Distribution and the Null Hypothesis H0:_mu1 - _mu0 = 0 (i.e. H0:_mu1 = _mu0)')
                                print(res.summary)
                                print('Z Score = ',res.result[1],'p value = ',res.result[3])
                                if (np.abs(res.result[1]) > 1.96) and (res.result[3] < 0.05):
                                    print('Difference between Test Arm and Control Arm is Statistically Significant')
                                    significant_flag = True
                                else:
                                    print('Difference between Test Arm and Control Arm is NOT Statistically Significant')
                                    significant_flag = False
                                cohens_d = (mu_1 - mu_0) / (np.sqrt((sigma_1 ** 2 + sigma_0 ** 2) / 2))                                
                                power = analysis.solve_power(cohens_d, power=None, nobs1=num_samples_1, ratio=num_samples_0/num_samples_1, alpha=0.05)
                                print('Statistical Power: %.3f' % power)    
                                
                                model_output_dict['stats_difference_comparision'] = {'zscore': res.result[1],
                                                                                       'p_value': res.result[3],
                                                                                       'Difference_between_test_and_control_Statistically_Significant': significant_flag,
                                                                                       'statistical_power': power}
                              
                                    
                                print('+++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')
                                
                                print('=====================================================')
                                
                                
                                
                            else:  
                                print('-----------------------------------------------------')
                                print('=====================================================')
                                print('Bayesian Regression Results')
                                print('=====================================================')  
                                if num_samples_0 == num_samples_1:
                                    test_df = pd.DataFrame({'test_data': test_distribution.rvs(size=num_samples_1, random_state=rng),
                                                            'control_data': control_distribution.rvs(size=num_samples_0, random_state=rng)})
                                
                                else:
                                    if num_samples_0 > num_samples_1:
                                        '''
                                        #=======================================
                                        # number of samples in control arm is greater than in the test arm
                                        #=======================================
                                        '''
                                        test_df = pd.DataFrame({'control_data': control_distribution.rvs(size=num_samples_0, random_state=rng)})
                                        rvs_test = test_distribution.rvs(size=num_samples_1, random_state=rng)
                                        tmp_arr = np.ones(num_samples_0)
                                        tmp_arr = np.where(tmp_arr==1, np.nan, tmp_arr)
                                        for i in range(num_samples_0):
                                            if i < num_samples_1:
                                                tmp_arr[i] = rvs_test[i] 
                                      
                                        test_df['test_data'] = tmp_arr
                                    else:
                                        '''
                                        #=======================================
                                        # number of samples in test arm is greater than in the control arm
                                        #=======================================
                                        '''
                                        test_df = pd.DataFrame({'test_data': test_distribution.rvs(size=num_samples_1, random_state=rng)})
                                        rvs_control = control_distribution.rvs(size=num_samples_0, random_state=rng)
                                        tmp_arr = np.ones(num_samples_1)
                                        tmp_arr = np.where(tmp_arr==1, np.nan, tmp_arr)
                                        for i in range(num_samples_1):
                                            if i < num_samples_0:
                                                tmp_arr[i] = rvs_control[i] 
                                      
                                        test_df['control_data'] = tmp_arr
                                
                                test_df['diff'] = test_df['test_data'] - test_df['control_data']  
                                
                                
                                res.ztest(df=test_df, 
                                         x = 'test_data', 
                                         y = 'control_data', 
                                        x_std = sigma_1, 
                                        y_std = sigma_0, 
                                         alpha = 0.05, 
                                         test_type = 2)
                                print('Two Sided Z Score Statistics between Test Arm Distribution and the Control Arm Distribution')                                
                                print(res.summary)
                                print('Z Score = ',res.result[2],'p value = ',res.result[4])
                                if (np.abs(res.result[2]) < 1.96):
                                    print('Test Arm Distribution and Control Arm Distribution are Statistically Similar')
                                    significant_flag = True
                                else:
                                    print('Test Arm Distribution and Control Arm Distribution are NOT Statistically Similar')
                                    significant_flag = False
                                model_output_dict['stats_distribution_comparision'] = {'zscore': res.result[2],
                                                                                       'p_value': res.result[4],
                                                                                       'Test_and_Control_Distributions_Statistically_Similar': significant_flag}
                                  
                                res.ztest(df=test_df, 
                                                 x = 'diff', 
                                                x_std = np.nanstd(test_df['diff'].values), 
                                                mu=0.0,
                                                 alpha = 0.05, 
                                                 test_type = 1)
                                print('Z-score Statistics between the Test Arm vs Control Arm Difference Distribution and the Null Hypothesis H0:_mu1 - _mu0 = 0 (i.e. H0:_mu1 = _mu0)')
                                print(res.summary)
                                print('Z Score = ',res.result[1],'p value = ',res.result[2])
                                if (np.abs(res.result[1]) > 1.96) and (res.result[2] < 0.05):
                                    print('Difference between Test Arm and Control Arm is Statistically Significant')
                                    significant_flag = True
                                else:
                                    print('Difference between Test Arm and Control Arm is NOT Statistically Significant')
                                    significant_flag = False
                                cohens_d = (mu_1 - mu_0) / (np.sqrt((sigma_1 ** 2 + sigma_0 ** 2) / 2))                                
                                power = analysis.solve_power(cohens_d, power=None, nobs1=num_samples_1, ratio=num_samples_0/num_samples_1, alpha=0.05)
                                print('Statistical Power: %.3f' % power)    
                                
                                
                            
                                print('=====================================================')  
                                model_output_dict['stats_difference_comparision'] = {'zscore': res.result[1],
                                                                                       'p_value': res.result[2],
                                                                                       'Difference_between_test_and_control_Statistically_Significant': significant_flag,
                                                                                       'statistical_power': power}
                              

                            
                            ax1.legend()  
                            json.dump(model_output_dict,open('output_model_run.json','w'))  
                            
                      
        plt.show()          
    except Exception as e:
        exc_traceback = sys.exc_info()
        err_message = ''
        for exc_line in traceback.format_exception(exc_traceback[0],
                                                   exc_traceback[1],
                                                   exc_traceback[2]):
            err_message = err_message+'\n\n'+exc_line
        logger.error(err_message)
        print(subprocess.call(["grep", "-r", 'ERROR', logfilename]))
        raise Exception(e)     
                
    
def get_power_prior(historic_df, current_df=None, data_to_test='', mean_value=None, alpha=0.05, k=5, model_run=''):
    '''
    .. note:: 
    
        **Function Details**
        
        Main Function to drive the Bayesian Statistical Analysis (with Power Prior consideration) for
        Medical Clinical trial data.
    
    Parameters
    ----------
    
    input_parameters: JSON (or file input), required        
        JSON representation of Model Input Parameters. This can be either a JSON string or a JSON formatted file.            
        
    build_documentation: Boolean, optional, default=False
        Flag to build documentation for this project.
    
    compare_historic_data: Boolean, optional, default=False
        Flag to perform a Historic vs Current Study data comparison analysis.    
    
    log_filename: str, optional, default="bayesian_model_output" 
        Filename to use for storing logging information for this model run.
        
    Returns
    -------
    
    output_dict: JSON
        JSON compatible dictionary object 
       
    
    .. note::
    
        find the Power Prior "λ" value to determine the proportion of data from historic studies to be included in the 
        statistical analyisis of the current study.
        
        Using the Power Prior proposed by :cite:t:`doi:10.1177/09622802221146309`
        
        λ = (A*P)^k 
        
        Where: 
            A = area under the overlap between the historic and current normal distributions 
            
            P = is the 'p-value' of the two-sided "z-test" with H0:_muh = _mu0 versus H1:_muh ≠ _mu0 
            
                (NB: if sample number is < 30 a "t_test" with H0:_muh = _mu0 versus H1:_muh ≠ _mu0 should be used) 
                
            k = a calibration parameter to control the shape and scale of λ
        
    '''
    print(historic_df.describe().T)
    print('==================')
    print(historic_df[data_to_test].to_markdown())
    print('==================')
    mu_1 = historic_df[data_to_test].mean()
    sigma_1 = historic_df[data_to_test].std()
    num_obs_1 = len(historic_df[data_to_test])
    
    res = stat()
    if current_df.empty:
        current_df = pd.DataFrame()
        mu_2 = mean_value
        sigma_2 = 1.0
        num_obs_2 = 1
    else:
        mu_2 = current_df[data_to_test].mean()
        sigma_2 = current_df[data_to_test].std()
        if sigma_2 == 0.0:
            sigma_2 = 1e-15
        num_obs_2 = len(current_df[data_to_test])
    
        
    if (len(historic_df) > 30) and (len(current_df) > 30):
        '''
        #=======================================================================
        # use "z-test" 
        #=======================================================================
        '''
        if len(historic_df[data_to_test].values) == len(current_df[data_to_test].values):
            test_df = pd.DataFrame({'historic_data': historic_df[data_to_test].values,
                                    'current_data': current_df[data_to_test].values})
        else:
            test_df = pd.DataFrame({'historic_data': historic_df[data_to_test].values})
            tmp_arr = np.ones_like(historic_df[data_to_test].values)
            tmp_arr = np.where(tmp_arr==1, np.nan, tmp_arr)
            for i in range(len(historic_df[data_to_test].values)):
                if i < len(current_df[data_to_test].values):
                    tmp_arr[i] = current_df[data_to_test].values[i] 
            
            test_df['current_data'] = tmp_arr
            
        historic_std = test_df['historic_data'].std()
        current_std = test_df['current_data'].std()
        print(historic_df[data_to_test])
        print(current_df[data_to_test])
        print(test_df)
        if len(historic_df[data_to_test].values) != len(current_df[data_to_test].values):
            evarFlag = False 
        else:
            evarFlag = True
        res.ztest(df=test_df, 
                         x = 'historic_data', 
                         y = 'current_data', 
                         x_std = historic_std, 
                         y_std = current_std, 
                         alpha = alpha, 
                         test_type = 2)
        p_value = res.result[4]
    else:
        '''
        #=======================================================================
        # use "t-test"
        #=======================================================================
        '''
        print('current_df',current_df)
        if current_df.empty:
            '''
            #===================================================================
            # perform one-sided t-test against the provided mean value (_mu)
            #===================================================================
            '''
            test_df = pd.DataFrame({data_to_test: historic_df[data_to_test].values,
                                    'data_type': 'historic'})
            print('t-test',test_df)
            res.ttest(df=test_df, 
                     res=data_to_test, 
                     alpha=alpha, 
                     test_type=1,
                     mu=mean_value)
            p_value = res.result[-2]
        else:
            '''
            #===================================================================
            # perform two-sided t-test
            #===================================================================
            '''
            test_df = pd.DataFrame({data_to_test: historic_df[data_to_test].values,
                                    'data_type': 'historic'})
            test_df = pd.concat([test_df, pd.DataFrame({data_to_test: current_df[data_to_test].values,
                                                        'data_type': 'current'})])
        
            
            print('t-test',test_df)
            if len(historic_df[data_to_test].values) != len(current_df[data_to_test].values):
                evarFlag = False 
            else:
                evarFlag = True
            res.ttest(df=test_df, 
                             xfac='data_type', 
                             res=data_to_test, 
                             evar=evarFlag, 
                             alpha=alpha, 
                             test_type=2)
            print('res',res.summary,res.result)
            p_value = res.result[-1]
    print('res',res.summary,res.result)
    print('p_value',p_value)
    

    '''
    #===========================================================================
    # compute Area under the overlapping normal distributions
    c=_mu2σ21−σ2(_mu1σ2+σ1(_mu1−_mu2)2+2(σ21−σ22)log(σ1σ2)−−−−−−−−−−−−−−−−−−−−−−−−−−√)σ21−σ22
    #===========================================================================
    '''
    
    #===========================================================================
    # mu_1 =current_df[data_to_test].mean()
    # mu_2 =  historic_df[data_to_test].mean()
    # sigma_1 = current_df[data_to_test].std()
    # sigma_2 = historic_df[data_to_test].std()
    #===========================================================================
    
    
    #===========================================================================
    # mu_1 = mu_2
    #===========================================================================
    
    #===========================================================================
    # std_1 = std_2
    #===========================================================================
    print('mu_1',mu_1,'mu_2',mu_2)
    print('sigma_1',sigma_1,'sigma_2',sigma_2)
    
    if sigma_1 != sigma_2:
        if (sigma_1 > 0) and (sigma_2 > 0):
            c = (mu_2*sigma_1**2 - sigma_2*(mu_1*sigma_2 + sigma_1*np.sqrt((mu_1 - mu_2)**2 + 2*(sigma_1**2 - sigma_2**2)*np.log(sigma_1/sigma_2)))) / (sigma_1**2 - sigma_2**2)
        else:
            c = (mu_1 + mu_2)/2
            
    else:
        c = (mu_1 + mu_2)/2
        
    print('c',c)
    #===========================================================================
    # if mu_1 <= mu_2:
    #     area_under_curve = 1 - 0.5*sp_special.erf((c-mu_1)/(np.sqrt(2)*sigma_1)) + 0.5*sp_special.erf((c-mu_2)/(np.sqrt(2)*sigma_2))
    # else:
    #     area_under_curve = 1 - 0.5*sp_special.erf((c-mu_2)/(np.sqrt(2)*sigma_2)) + 0.5*sp_special.erf((c-mu_1)/(np.sqrt(2)*sigma_1))
    #===========================================================================
    if (sigma_1 > 0) and (sigma_2 > 0):
        area_under_curve = 1 - 0.5*sp_special.erf((c-mu_1)/(np.sqrt(2)*sigma_1)) + 0.5*sp_special.erf((c-mu_2)/(np.sqrt(2)*sigma_2))
    else:
        area_under_curve = 0.0     
    print('area_under_curve 2: ',area_under_curve)
    k_value = k
    prior_lambda_ak = (area_under_curve)**k_value 
    prior_lambda_apk = (area_under_curve*p_value)**k_value 
    print('prior_lambda (A)^k',prior_lambda_ak)
    print('prior_lambda (AP)^k',prior_lambda_apk)
    
    
    cohens_d = (mu_1 - mu_2) / (np.sqrt((sigma_1 ** 2 + sigma_2 ** 2) / 2))
    print('cohens_d',cohens_d)
    analysis = TTestIndPower()
    power = analysis.solve_power(cohens_d, power=None, nobs1=num_obs_1, ratio=num_obs_2/num_obs_1, alpha=alpha)
    print('historic power: %.3f' % power)
    #===========================================================================
    # result = analysis.solve_power(cohens_d, power=None, nobs1=num_obs_2, ratio=1.0, alpha=alpha)
    # print('current power: %.3f' % result)
    #===========================================================================
    
    result_dict = {'model_run': str(model_run),
                   'area_under_curve': area_under_curve,
                   'distribution_intercept': c,
                   'data_to_test': data_to_test,
                   'p_value': p_value, 
                   'alpha': alpha,
                   'k': k_value,
                   'power': power,
                   'prior_lambda_ak': prior_lambda_ak, 
                   'prior_lambda_apk':prior_lambda_apk}
    return result_dict
   

def get_prior(mu, 
              distribution_type='uniform',
              norm_mean=0,
              norm_sigma=1):
    '''
    
    P(θ)
    
    '''
    if distribution_type == 'uniform':
        prior_dist = sp_stats.uniform.pdf(mu, loc=norm_mean) + 1 
    elif distribution_type == 'normal':
        prior_dist = sp_stats.norm.pdf(mu, loc=norm_mean, scale= norm_sigma) + 1 

    prior_dist = prior_dist/prior_dist.sum() 
        
    return prior_dist
 

def get_likelihood(datum,
                      data_values,
                      sigma):
    '''
    
        P(Data|θ)
        
    '''
    likelihood_out = sp_stats.norm.pdf(datum, data_values, scale = sigma) 
    return likelihood_out/likelihood_out.sum()

    #===========================================================================
    # data_mean = np.mean(data_values)
    # 
    # likelihood = (2*np.pi*sigma**2)^((-1*num_samples)/2) * np.exp(-np.sum((data_values - data_mean)**2)/(2*sigma**2))
    # return likelihood
    #===========================================================================

def get_P_data(unnormalized_posterior,
               mu):
    '''
    
       P(Data) = ∫P(Data|θ)∗P(θ)dθ
    
    '''
    p_data = sp_integrate.trapz(unnormalized_posterior, mu)
    
    return p_data 


    
def compute_percentile(parameter_values, distribution_values, percentile):
    cumulative_distribution = sp_integrate.cumtrapz(distribution_values, 
                                                    parameter_values)
    percentile_index = np.searchsorted(cumulative_distribution, 
                                       percentile)
    #===========================================================================
    # print(cumulative_distribution)
    # print(parameter_values)
    # print(percentile)
    # print(percentile_index)
    # print(parameter_values[percentile_index])
    # #===========================================================================
    # # plt.figure(figsize = (10, 8))
    # # plt.plot(parameter_values, cumulative_distribution)
    # # plt.show()
    # #===========================================================================
    # print(asdf)
    #===========================================================================
   
    
    return parameter_values[percentile_index]
    

def find_posterior(
                   input_data,
                   prior_dist_type='uniform',
                   historic_data=[],
                   power_prior_lambda=0.0,
                   multiple_historic_datasets_flag=False,
                   plot_data=True):
    '''
    .. note:: 
    
        **Function Details**
        
        Function to compute the Posterior Distribution for the following scenarios:
        
        1) Standard Posterior Distribution Computation (i.e. Current Study Data only);
        2) Power Prior Posterior Distribution Computation for a single historic dataset; and,
        3) Power Prior Posterior Distribution Computation for multiple historic datasets.
        
        
        The standard Posterior Relation is: 
        
        
        P(θ|Data) ∝ L(θ|Data) * P(θ)
        
        
        Where:
            P(θ|Data) = The Posterior Distribution θ, factoring in the data from the current study
            
            L(θ|Data) = Likelihood Function of the Data for the Current Study = P(Data|θ) / P(Data)
            
            P(θ) = The initial Prior Distribution Function
            
            
            
        When including the influence of the Power Prior factor (λ) to include historical data this becomes:
        
        
    
        P(θ|Data, Data0, λ) ∝ L(θ|Data) * (L(θ|Data0)^λ * P(θ)).
        
        
        Where:  
            P(θ|Data, Data0, λ) = The Posterior Distribution θ, factoring in both the current and historical data and the Power Prior Factor.
            
            L(θ|Data) = Likelihood Function of the Data for the Current Study = P(Data|θ) / P(Data)
            
            L(θ|Data0) = Likelihood Function for the Data of a Historic Study = P(Data0|θ) / P(Data0)
            
            P(θ) = The initial Prior Distribution Function
            
            λ = Power Prior Factor
          
          
          
        For multiple historical datasets the above relation becomes:
                      
        P(θ|Data, Data0, λ) ∝ L(θ|Data) *  ∏ (L(θ|Data0)^λ * P(θ))
                                            
        Where: 
            ∏ = the "Prod" function from k=1 to k=k0 ("Prod" is the multiplication version of "Sum")
            
            k0 = the number of historical datasets

    
    Parameters
    ----------
    
    input_data: Array-like, required        
        An array of input data values representing the "Current Study" dataset            
    
    historic_data: Array-like, optional, default=[]
        An array of input data values representing the "Historic Study" dataset. If the historic data object
        represents multiple historical datasets this must be a 2-dimensional array.
        
    prior_dist_type: str, optional, default='uniform'
        Prior Distribution Type - either "uniform" or "normal"
    
    power_prior_lambda: number, optional, default=0.0
        Power Prior Parameter to use in the computation. Must be between 0 and 1.
    
    multiple_historic_datasets_flag: Boolean, optional, default=False
        Flag to perform the computation for multiple historic datasets.    
    
    plot_data: Boolean, optional, default=True
        Flag to plot the data for the Bayesian Regression.    
    
        
    Returns
    -------
    
    output: tuple
        A tuple object including the "posterior mean", "posterior standard deviation", and the number of historical data samples borrowed 
       
        
    '''
    
    if plot_data:
        mosaic = [["prior", "posterior_models"],
                  ["regression", "final_model"]]
        
        fig, axes = plt.subplot_mosaic(
            mosaic,
        )
        fig.suptitle("Bayesian Posterior Distribution Regression", fontsize=16)
        fig.set_layout_engine('constrained')
    
    if (np.min(input_data) == 0.0) and (np.max(input_data) == 0.0):
        mu = np.linspace(0, 1, num = len(input_data))
        
    else:
        mu = np.linspace(np.min(input_data), np.max(input_data), num = len(input_data))
    
    '''
    #===========================================================================
    # Code Validation Step - running a baseline Bayesian Regression without including the power prior
    #===========================================================================
    '''
    if prior_dist_type == 'uniform':
        prior = get_prior(mu)
    elif prior_dist_type == 'normal':
        if np.std(input_data) == 0.0:
            '''
            #===================================================================
            # for the case of a binary flagged dataset with values of 0's and 1's where the 
            # observed standard deviation is 0.0, we will assume a 
            # nominal standard deviation of 0.5 - which is half of the sampling interval
            #===================================================================
            '''
            prior = get_prior(mu, 
                              distribution_type=prior_dist_type, 
                              norm_mean=np.mean(input_data), 
                              norm_sigma=0.5)
        
        else:
            prior = get_prior(mu, 
                              distribution_type=prior_dist_type, 
                              norm_mean=np.mean(input_data), 
                              norm_sigma=np.std(input_data))
        
        prior = get_prior(mu, 
                          distribution_type=prior_dist_type, 
                          norm_mean=np.mean(input_data), 
                          norm_sigma=np.std(input_data))
        power_prior = prior
       
    
    max_historic_data = 0
    min_historic_data = 0
    historic_data_borrowed = []
    if multiple_historic_datasets_flag:
        historic_priors = []
        data_id = 0
        for historic_dataset in historic_data: 
            if np.min(historic_dataset) < min_historic_data:
                min_historic_data = np.min(historic_dataset)
            if np.max(historic_dataset) > max_historic_data:
                max_historic_data = np.max(historic_dataset)
                    
            if len(historic_dataset)>0:
                historic_mu = np.linspace(min_historic_data, max_historic_data, num = len(historic_dataset))
            else:
                historic_mu = np.array([])
            
            historic_likelihood_dict = {}
            #===================================================================
            # historic_prior = get_prior(historic_mu)
            #===================================================================
            
            if prior_dist_type == 'uniform':
                historic_prior = get_prior(historic_mu)
            elif prior_dist_type == 'normal':
                if np.std(historic_dataset) == 0.0:
                    '''
                    #===================================================================
                    # for the case of a binary flagged dataset with values of 0's and 1's where the 
                    # observed standard deviation is 0.0, we will assume a 
                    # nominal standard deviation of sqrt(0.5^2) - which is the standard deviation of the sampling uncertainty (i.e. half of the sampling interval)
                    #===================================================================
                    '''
                    historic_prior = get_prior(historic_mu, 
                                      distribution_type=prior_dist_type, 
                                      norm_mean=np.mean(historic_dataset), 
                                      norm_sigma=np.sqrt(0.5**2))
                
                else:
                    historic_prior = get_prior(historic_mu, 
                                      distribution_type=prior_dist_type, 
                                      norm_mean=np.mean(historic_dataset), 
                                      norm_sigma=np.std(historic_dataset))
                
            for h_ind, h_datum in enumerate(historic_dataset):
                historic_likelihood = get_likelihood(h_datum, historic_mu, np.std(historic_mu))  
                historic_likelihood = historic_likelihood**power_prior_lambda     
                historic_likelihood_dict[h_ind] = historic_likelihood                
                unnormalized_posterior = historic_likelihood * historic_prior
                p_data = get_P_data(unnormalized_posterior, historic_mu)
                if p_data == 0.0:
                    normalized_posterior = unnormalized_posterior*p_data
                else:
                    normalized_posterior = unnormalized_posterior/p_data
                historic_prior = normalized_posterior
                historic_likelihood_dict[h_ind] = normalized_posterior
                
            num_samples_to_borrow = int(power_prior_lambda*len(historic_dataset))
            borrowed_data = np.random.choice(historic_dataset, size=num_samples_to_borrow, replace=False)
            historic_data_borrowed.append(borrowed_data)
            
            
            if len(historic_dataset) > 0:   
                if plot_data: 
                    axes['posterior_models'].plot(historic_mu, 
                             historic_likelihood_dict[len(historic_dataset)-1], 
                             color = 'r',
                             linestyle='dashed',
                             linewidth=2,
                             label = f'Power Prior Model final  - '+str(data_id))
            
            #===================================================================
            # 
            #     print('posterior_mean',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.5))
            #     print('posterior_std',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.5-0.341),compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.5+0.341))
            #     print('posterior_99.5',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.995))
            #     print('posterior_0.5',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.005))
            #     
            #     print('first model posterior_mean',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.5))
            #     print('first model posterior_std',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.5-0.341),compute_percentile(historic_mu, historic_likelihood_dict[0], 0.5+0.341))
            #     print('first model posterior_99.5',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.995))
            #     print('first model posterior_0.5',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.005))
            #     
            #===================================================================
                historic_mean = compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.5)
                historic_std = (compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.5+0.341) - compute_percentile(historic_mu, historic_likelihood_dict[len(historic_dataset)-1], 0.5-0.341))/2
                #===============================================================
                # print('historic_mean',historic_mean,'historic_std',historic_std)
                #===============================================================
                
                power_prior =  get_prior(mu, 
                                      distribution_type=prior_dist_type, 
                                      norm_mean=historic_mean, 
                                      norm_sigma=historic_std)
                if not np.any(np.isnan(power_prior)):
                    historic_priors.append(power_prior)
                if plot_data:
                    axes['prior'].plot(mu, power_prior, label = 'power_prior Distribution - '+str(data_id))
            data_id += 1
        power_prior = np.prod(historic_priors, axis=0)
        print('historic_data_borrowed',historic_data_borrowed)
        if plot_data:
            axes['prior'].plot(mu, power_prior, label = 'Master power_prior Distribution')
    else:
        historic_likelihood_dict = {}
        if len(historic_data)>0:
            historic_mu = np.linspace(np.min(historic_data), np.max(historic_data), num = len(historic_data))
        else:
            historic_mu = np.array([])
        
        if prior_dist_type == 'uniform':
            historic_prior = get_prior(historic_mu)
        elif prior_dist_type == 'normal':
            if np.std(historic_data) == 0.0:
                '''
                #===================================================================
                # for the case of a binary flagged dataset with values of 0's and 1's where the 
                # observed standard deviation is 0.0, we will assume a 
                # nominal standard deviation of 0.5 - which is half of the sampling interval
                #===================================================================
                '''
                historic_prior = get_prior(historic_mu, 
                                  distribution_type=prior_dist_type, 
                                  norm_mean=np.mean(historic_data), 
                                  norm_sigma=np.sqrt(0.5**2))
            
            else:
                historic_prior = get_prior(historic_mu, 
                                  distribution_type=prior_dist_type, 
                                  norm_mean=np.mean(historic_data), 
                                  norm_sigma=np.std(historic_data))
                
            historic_prior = get_prior(historic_mu, 
                              distribution_type=prior_dist_type, 
                              norm_mean=np.mean(historic_data), 
                              norm_sigma=np.std(historic_data))
            
        
        for h_ind, h_datum in enumerate(historic_data):
            historic_likelihood = get_likelihood(h_datum, historic_mu, np.std(historic_mu))  
            historic_likelihood = historic_likelihood**power_prior_lambda     
            unnormalized_posterior = historic_likelihood * historic_prior
            p_data = get_P_data(unnormalized_posterior, historic_mu)
            normalized_posterior = unnormalized_posterior/p_data
            historic_prior = normalized_posterior
            historic_likelihood_dict[h_ind] = normalized_posterior
            
            if plot_data:
                if h_ind%(len(input_data)//10) == 0:
                    axes['posterior_models'].plot(historic_mu, normalized_posterior, label = f'Model after observing {h_ind} samples')
                
            num_samples_to_borrow = int(power_prior_lambda*len(historic_data))
            borrowed_data = np.random.choice(historic_data, size=num_samples_to_borrow, replace=False)
            historic_data_borrowed.append(borrowed_data)
            
        if len(historic_data) > 0:   
            if plot_data: 
                axes['posterior_models'].plot(historic_mu, 
                         historic_likelihood_dict[len(historic_data)-1], 
                         color = 'r',
                         linestyle='dashed',
                         linewidth=2,
                         label = f'Power Prior Model final')
        
        
            #===================================================================
            # print('posterior_mean',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.5))
            # print('posterior_std',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.5-0.341),compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.5+0.341))
            # print('posterior_99.5',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.995))
            # print('posterior_0.5',compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.005))
            # 
            # print('first model posterior_mean',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.5))
            # print('first model posterior_std',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.5-0.341),compute_percentile(historic_mu, historic_likelihood_dict[0], 0.5+0.341))
            # print('first model posterior_99.5',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.995))
            # print('first model posterior_0.5',compute_percentile(historic_mu, historic_likelihood_dict[0], 0.005))
            #===================================================================
            historic_mean = compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.5)
            historic_std = (compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.5+0.341) - compute_percentile(historic_mu, historic_likelihood_dict[len(historic_data)-1], 0.5-0.341))/2
            print('historic_mean',historic_mean,'historic_std',historic_std)
            
            power_prior =  get_prior(historic_mu, 
                                  distribution_type=prior_dist_type, 
                                  norm_mean=historic_mean, 
                                  norm_sigma=historic_std)
            
            if plot_data:
                axes['prior'].plot(historic_mu, power_prior, label = 'power_prior Distribution')
            
    posterior_dict = {}    
    for ind, datum in enumerate(input_data):
        likelihood = get_likelihood(datum, mu, np.std(mu))
        unnormalized_posterior = likelihood * prior
        p_data = get_P_data(unnormalized_posterior, mu)
        normalized_posterior = unnormalized_posterior/p_data
        prior = normalized_posterior
        posterior_dict[ind] = normalized_posterior
        
        if plot_data:
            if ind%(len(input_data)//5) == 0:
                axes['posterior_models'].plot(mu, normalized_posterior, label = f'Baseline Model after observing {ind} samples')
    if plot_data:    
        axes['posterior_models'].plot(mu, 
                 posterior_dict[len(input_data)-1], 
                 color = 'k',
                 linestyle='dashed',
                 linewidth=2,
                 label = f'Baseline Model final')
    
    
        axes['posterior_models'].legend()
        axes['posterior_models'].set_title('Bayesian Posterior Models')
        axes['posterior_models'].set_xlabel("Modelled Data Values")
        axes['posterior_models'].set_ylabel("Modelled Probability Density")
     
#===============================================================================
#     print('mu',mu)
#     print('mu_mean',np.mean(mu))
#     print('mu_std',np.std(mu))
#           
#     print('prior_norm_mean',compute_percentile(mu, prior, 0.5))
#     print('prior_norm_std_lower',compute_percentile(mu, prior, 0.5-0.341))
#     print('prior_norm_std_upper',compute_percentile(mu, prior, 0.5+0.341))
#     print('posterior_mean',compute_percentile(mu, posterior_dict[len(input_data)-1], 0.5))
#     print('posterior_std',compute_percentile(mu, posterior_dict[len(input_data)-1], 0.5-0.341),compute_percentile(mu, posterior_dict[len(input_data)-1], 0.5+0.341))
#     print('posterior_99.5',compute_percentile(mu, posterior_dict[len(input_data)-1], 0.995))
#     print('posterior_0.5',compute_percentile(mu, posterior_dict[len(input_data)-1], 0.005))
#     
#     print('first model posterior_mean',compute_percentile(mu, posterior_dict[0], 0.5))
#     print('first model posterior_std',compute_percentile(mu, posterior_dict[0], 0.5-0.341),compute_percentile(mu, posterior_dict[0], 0.5+0.341))
#     print('first model posterior_99.5',compute_percentile(mu, posterior_dict[0], 0.995))
#     print('first model posterior_0.5',compute_percentile(mu, posterior_dict[0], 0.005))
# 
#     print('test_arr_mean',np.mean(input_data))
#     print('test_arr_std',np.std(input_data))
#     print('test_arr_var',np.var(input_data))
#===============================================================================

    if len(historic_data_borrowed) > 0:
        master_borrowed_data = np.hstack(historic_data_borrowed)
        borrowed_input_data = np.hstack([master_borrowed_data,input_data])
    else:
        borrowed_input_data = input_data
    
    if len(borrowed_input_data)>0:
        print('min',np.min(borrowed_input_data))
        print('max',np.max(borrowed_input_data))
        print('len(borrowed_input_data)',len(borrowed_input_data))
        borrowed_mu = np.linspace(np.min(borrowed_input_data), np.max(borrowed_input_data), num = len(borrowed_input_data))
    else:
        borrowed_mu = np.array([])
   
    if prior_dist_type == 'uniform':
        borrowed_prior = get_prior(borrowed_mu)
    elif prior_dist_type == 'normal':
        if np.std(input_data) == 0.0:
            '''
            #===================================================================
            # for the case of a binary flagged dataset with values of 0's and 1's where the 
            # observed standard deviation is 0.0, we will assume a 
            # nominal standard deviation of np.sqrt(0.5**2) - which is half of the sampling interval
            #===================================================================
            '''
            borrowed_prior = get_prior(borrowed_mu, 
                              distribution_type=prior_dist_type, 
                              norm_mean=np.mean(borrowed_input_data), 
                              norm_sigma=np.sqrt(0.5**2))
         
        else:
            borrowed_prior = get_prior(borrowed_mu, 
                              distribution_type=prior_dist_type, 
                              norm_mean=np.mean(borrowed_input_data), 
                              norm_sigma=np.std(borrowed_input_data))
            
        borrowed_prior = get_prior(borrowed_mu, 
                          distribution_type=prior_dist_type, 
                          norm_mean=np.mean(borrowed_input_data), 
                          norm_sigma=np.std(borrowed_input_data))
    
    if plot_data:
        axes['prior'].plot(mu, power_prior, label = 'power_prior Distribution')
                
        axes['prior'].legend()      
        axes['prior'].set_title('Bayesian Prior Model')
        axes['prior'].set_xlabel("Prior Data Values")
        axes['prior'].set_ylabel("Prior Probability Density")
        axes['prior'].set_yscale("log")
    
    #===========================================================================
    # borrowed_prior = power_prior
    #===========================================================================
    
    prior_posterior_dict = {}    
    for ind, datum in enumerate(borrowed_input_data):
        likelihood = get_likelihood(datum, borrowed_mu, np.std(borrowed_mu))
        unnormalized_posterior = likelihood * borrowed_prior
        p_data = get_P_data(unnormalized_posterior, borrowed_mu)
        normalized_posterior = unnormalized_posterior/p_data
        borrowed_prior = normalized_posterior
        prior_posterior_dict[ind] = normalized_posterior
        if plot_data:
            if ind%(len(borrowed_input_data)//5) == 0:
                axes['posterior_models'].plot(borrowed_mu, normalized_posterior, label = f'Main Model after observing {ind} samples')
            
    
    if plot_data: 
        axes['posterior_models'].plot(borrowed_mu, 
                 prior_posterior_dict[len(borrowed_input_data)-1], 
                 color = 'b',
                 linestyle='dashed',
                 linewidth=2,
                 label = f'Bayesian Model final')
        
        axes['posterior_models'].legend()
      
    num_x_samples= len(borrowed_input_data)
    
    top_bound = means = [compute_percentile(borrowed_mu, prior_posterior_dict[i], 0.995) for i in range(num_x_samples)]
    bottom_bound = means = [compute_percentile(borrowed_mu, prior_posterior_dict[i], 0.005) for i in range(num_x_samples)]
    means = [compute_percentile(borrowed_mu, prior_posterior_dict[i], 0.5) for i in range(num_x_samples)]
    
    if plot_data:
        axes['regression'].plot(range(num_x_samples), top_bound, label = '99.5% Confidence - top bound')
        axes['regression'].plot(range(num_x_samples), bottom_bound, label = '99.5% Confidence - bottom bound')
        axes['regression'].plot(range(num_x_samples), means, label = 'regression mean')
        axes['regression'].legend()
        axes['regression'].set_title('Model Convergence')
        axes['regression'].set_xlabel("Number of samples used")
        axes['regression'].set_ylabel("Predicted Mean $\mu$")
        
           
        axes['final_model'].plot(mu, posterior_dict[len(input_data)-1], label = 'Final Posterior Model')
        axes['final_model'].plot(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], label = 'Final Power Posterior Model')
        axes['final_model'].axvline(x = compute_percentile(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], 0.5), ls = '--', color = 'r', label = 'Prior Model Mean')
        axes['final_model'].axvline(x = compute_percentile(mu, posterior_dict[len(input_data)-1], 0.5), ls = '--', color = 'k', label = 'Model Mean')
        axes['final_model'].axvline(x = compute_percentile(mu, posterior_dict[len(input_data)-1], 0.005), ls = '--', color = 'y', label = '99% Confidence Interval')
        axes['final_model'].axvline(x = compute_percentile(mu, posterior_dict[len(input_data)-1], 0.995), ls = '--', color = 'y')
        axes['final_model'].axvline(x = compute_percentile(mu, posterior_dict[len(input_data)-1], 0.5-0.341), ls = '--', color = 'g', label = 'Standard Deviation')
        axes['final_model'].axvline(x = compute_percentile(mu, posterior_dict[len(input_data)-1], 0.5+0.341), ls = '--', color = 'g')
        axes['final_model'].legend()
        axes['final_model'].set_title('Final Posterior Model')
        axes['final_model'].set_xlabel("Modelled Data Values")
        axes['final_model'].set_ylabel("Modelled Probability Density")
        axes['final_model'].set_xlim(compute_percentile(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], 0.00005),
                                     compute_percentile(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], 0.99995))
        
        
        plt.show()
    
    
    posterior_mean = compute_percentile(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], 0.5)
    posterior_std = (compute_percentile(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], 0.5+0.341) - compute_percentile(borrowed_mu, prior_posterior_dict[len(borrowed_input_data)-1], 0.5-0.341))/2
    
    print('posterior_mean',posterior_mean)
    print('posterior_std',posterior_std)
    return posterior_mean, posterior_std, len(borrowed_input_data)
    
    
def identify_axes(ax_dict, fontsize=48):
    """
    Helper to identify the Axes in the examples below.

    Draws the label in a large font in the center of the Axes.

    Parameters
    ----------
    ax_dict : dict[str, Axes]
        Mapping between the title / label and the Axes.
    fontsize : int, optional
        How big the label should be.
    """
    kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
    for k, ax in ax_dict.items():
        ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)
        
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Bayesian Analysis of Control Studies with Power Prior')    
    parser.add_argument('-m','--model_data',
                        help='Flag to perform the Bayesian Regression Model on the data. [default = False]',
                        default=False,
                        action='store_true')
    parser.add_argument('-c','--compare_historic_data',
                        help='Flag to analyse Historic vs Current data and determine the Power Prior parameter. [default = False]',
                        default=False,
                        action='store_true')
    parser.add_argument('-d','--build_documentation',
                        help='Flag to build documenation for this application/project. [default = False]',
                        default=False,
                        action='store_true')
    parser.add_argument('-i','--input_parameters',
                        help='Model Parameters for the selected Bayesian Analysis run in JSON data format - (e.g. {"historic_data": "all_historic_df", "current_data": "study_data", "datafields": ["percent_survived"], "alpha_levels": [0.025], "k_values": [3, 5, 8], "mean_values": [67.5], "mouse_line": "FVB/N", "mouse_age_weeks": [52, 150]})',
                        required=False)
    parser.add_argument('--log_filename',
                        '-l',
                        help='The output log filename from the model run.',
                        required=False)
    args = parser.parse_args()
    print(args)
    main(model_data=args.model_data,
         compare_historic_data=args.compare_historic_data,
         build_documentation=args.build_documentation,
         input_parameters=args.input_parameters,
         log_filename=args.log_filename)
    
    
    #===========================================================================
    # (all_historic_df, 
    #                                   all_historic_df,
    #                                   datafields = ['percent_survived','total_mice_end','study_duration_weeks'],
    #                                   alpha_levels= [0.025,0.05],
    #                                   k_values=[1,2,3,4,5,6,7,8])
    #===========================================================================
    
    
    
    
