def active_learning(regressor, gold_regressor, allowed_conc, test_size = 100, steps = 10, verbose=0):
    ## first step
    if verbose:
        print('step:  1')    
    # make first dataset
    X_train_1 = random_input(allowed_conc, test_size)
    
    # first fit
    regressor.fit(X_train_1, gold_regressor.predict(X_train_1))
    
    # save results
    result = pd.DataFrame(X_train_1)
    result['gold_yield'] = gold_regressor.predict(X_train_1)
    result['pred_yield'] = 0.0 # not available but choose 0.0 to avoid further error
    result['step'] = 'step_1'
    
    ## next steps loop
    for step in range(steps-1):
        if verbose>=2:
            print('step: ',step+2)
        # make i th dataset
        X_train_1_1 = random_input(allowed_conc, 100000)
        df_1 = pd.DataFrame(X_train_1_1)
        df_1['pred_yield'] = regressor.predict(X_train_1_1)
        df_1 = df_1.sort_values(['pred_yield'], ascending=False)
        X_train_2 = df_1.iloc[0:test_size,0:11].values
        
        # save and add results
        temp_result = pd.DataFrame(X_train_2)
        temp_result['gold_yield'] = gold_regressor.predict(X_train_2)
        temp_result['pred_yield'] = df_1.iloc[0:test_size,11:12].values
        temp_result['step'] = 'step_{}'.format(step+2)
        result = pd.concat([result, temp_result], ignore_index=True)
        
        # update and refit regressor
        regressor.fit(result.iloc[:,0:11].values, result.iloc[:,11].values)

    return result, regressor


def bayesian_optimization(regressors_list,
                          gold_regressor,
                          allowed_conc,
                          exploitation=1, exploration=1, test_size=100, steps=10, verbose=0):
    ## first step
    if verbose:
        print('step:  1')    
    # make first dataset
    X_train_1 = random_input(allowed_conc, test_size)
    
    # first fit
    for regressor in regressors_list:
        regressor.fit(X_train_1, gold_regressor.predict(X_train_1))
    
    # save results
    result = pd.DataFrame(X_train_1)
    result['gold_yield'] = gold_regressor.predict(X_train_1)
    result['pred_yield'] = 0.0 # not available but choose 0.0 to avoid further error
    result['step'] = 'step_1'
    
    ## next steps loop
    for step in range(steps-1):
        if verbose>=2:
            print('step: ',step+2)
        # make i th dataset
        X_train_1_1 = random_input(allowed_conc, 100000)
        df_1 = pd.DataFrame(X_train_1_1)
        
        #upper Confidence Bound
        for index, regressor in enumerate(regressors_list):
            df_1['pred_yield_{}'.format(index)] = regressor.predict(X_train_1_1)
        
        df_1['regressors_std'] = df_1[[str(i) for i in df_1.columns if 'pred_yield' in str(i)]].std(axis=1)
        df_1['mean_vote'] = df_1[[str(i) for i in df_1.columns if 'pred_yield' in str(i)]].mean(axis=1)
        df_1['UCB'] = exploitation * df_1['mean_vote']+ exploration * df_1['regressors_std'] 
        df_1 = df_1.sort_values(['UCB'], ascending=False)
        X_train_2 = df_1.iloc[0:test_size,0:11].values
        
        # save and add results
        temp_result = pd.DataFrame(X_train_2)
        temp_result['gold_yield'] = gold_regressor.predict(X_train_2)
        #temp_result['pred_yield'] = df_1.iloc[0:test_size,11:12].values
        temp_result['pred_yield'] = df_1.mean_vote[0:test_size].values
        temp_result['step'] = 'step_{}'.format(step+2)
        result = pd.concat([result, temp_result], ignore_index=True)
        
        # update and refit regressor
        regressor.fit(result.iloc[:,0:11].values, result.iloc[:,11].values)

    return result, regressor