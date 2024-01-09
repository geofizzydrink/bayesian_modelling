
import os
import json


def main():
    '''
    #===========================================================================
    # Run two pass bayesian regression: 
    # 
    # Pass 1: Comparison of Current vs Historic control data 
    # Pass 2: Model Regression Run using the Power Prior Lambda value obtained from Pass 1.
    #===========================================================================
    '''
    
    '''
    #===========================================================================
    # Pass 1
    #===========================================================================
    '''
    
    os.system('python bayesian_prior_prediction.py -c -i default_historic_data_comparison_input.json')
    '''
    #===========================================================================
    # get power prior lambda values
    #===========================================================================
    '''
    with open('output_comparison_run.out','r') as fd:
        comparison_output_data = json.load(fd)
        print(comparison_output_data)
        print(comparison_output_data['comparison_run_output'])
        
        prior_lambda_ak = comparison_output_data['comparison_run_output'][0]['prior_lambda_ak']
        prior_lambda_apk = comparison_output_data['comparison_run_output'][0]['prior_lambda_apk']
        print('prior_lambda_ak',prior_lambda_ak)
        print('prior_lambda_apk',prior_lambda_apk)
       
       
    '''   
    #===========================================================================
    # Pass 2 Model Regression Run    
    #===========================================================================
    '''
    model_input_filename = 'default_model_input.json'
    with open(model_input_filename) as fp:
        model_input_json = json.load(fp)
    final_model_run_output_json = {}
    '''    
    #===========================================================================
    # (a) using A^k Prior Lambda value (area under curve only)
    #===========================================================================
    '''        
    '''    
    #===========================================================================
    # update "default_model_input.json" model input parameter file with "prior_lambda_ak" prior lambda value
    #===========================================================================
    '''    
    model_input_json[0]['power_prior_lambda'] = prior_lambda_ak
    json.dump(model_input_json,open('tmp_model_run_a.json','w'))
    os.system('python bayesian_prior_prediction.py -m -i tmp_model_run_a.json')  
    model_run_output_json = json.load(open('output_model_run.json','r'))   
    final_model_run_output_json['prior_lambda_ak-output'] = model_run_output_json
    '''    
    #===========================================================================
    # (b) using A.p^k Moderated Prior Lambda value (area under curve moderated by the statistical p-value)
    #     - this is a more conservative estimate of the Prior Lambda Value 
    #===========================================================================
    '''    
  
    '''    
    #===========================================================================
    # update "default_model_input.json" model input parameter file with "prior_lambda_apk" prior lambda value
    #===========================================================================
    '''      
    model_input_json[0]['power_prior_lambda'] = prior_lambda_apk
    json.dump(model_input_json,open('tmp_model_run_b.json','w'))
    os.system('python bayesian_prior_prediction.py -m -i tmp_model_run_b.json') 
    model_run_output_json = json.load(open('output_model_run.json','r'))   
    final_model_run_output_json['prior_lambda_apk-output'] = model_run_output_json
    '''    
    #===========================================================================
    # (c) using Prior Lambda Value of 0.0 (i.e. no borrowing of historic data) 
    #===========================================================================
    '''    
  
    '''    
    #===========================================================================
    # update "default_model_input.json" model input parameter file with prior lambda value = 0.0
    #===========================================================================
    '''      
    model_input_json[0]['power_prior_lambda'] = 0.0
    json.dump(model_input_json,open('tmp_model_run_c.json','w'))
    os.system('python bayesian_prior_prediction.py -m -i tmp_model_run_c.json') 
    model_run_output_json = json.load(open('output_model_run.json','r'))   
    final_model_run_output_json['no_borrow-output'] = model_run_output_json
    
    
    
    json.dump(final_model_run_output_json,open('model_run_for_paper_output.out','w'),indent=4)
    print(json.dumps(final_model_run_output_json,indent=4))

if __name__ == '__main__':
    
    main()