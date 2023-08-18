# based on scripts from limix_qtl

def define_correction_function(top_pvalues_perm, cis_mode):
    #Always try to use the MLE estimator, new default to 10 permutations.
    #If the MLE estimator fails we go back to the cruder estimation of the beta distribution.
    offset = (np.finfo(np.double).tiny*100)
    ##Replace zero's value with smallest number not 0.
    top_pvalues_perm[top_pvalues_perm == 0] = offset
    ##Replace highest value with highest number not 1.
    top_pvalues_perm[top_pvalues_perm == 1] = 1-offset
    try :
        alpha_para,beta_para,loc,fscale =  beta.fit(top_pvalues_perm,floc=0,fscale=1)
    except (scipy.stats._continuous_distns.FitSolverError):
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    except (scipy.stats._continuous_distns.FitDataError):
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    if(cis_mode):
        if(alpha_para<BETA_SHAPE1_MIN or alpha_para>BETA_SHAPE1_MAX or alpha_para<BETA_SHAPE2_MIN_CIS or alpha_para>BETA_SHAPE2_MAX_CIS):
            alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
            ### If pvalues become more significant after multiple testing correction we put them back to the orignal test Pvalue in a seperate step.
    else :
        if(alpha_para<BETA_SHAPE1_MIN or alpha_para>BETA_SHAPE1_MAX or alpha_para<BETA_SHAPE2_MIN_TRANS or alpha_para>BETA_SHAPE2_MAX_TRANS):
            alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    
    beta_dist = scipy.stats.beta(alpha_para,beta_para)
    correction_function = lambda x: beta_dist.cdf(x)
    #Would be good to replace 0 with minimal double value of python.
    return [correction_function, alpha_para, beta_para]

def apply_pval_correction(self,feature_id,top_pvalues_perm,cis_mode):
    '''Function to correct p values based on nominal p values and the top
    hits from permutation runs for the given feature.'''
    table = self.h5file.get_node('/'+feature_id)
    if(np.mean(top_pvalues_perm)>=0.999999999 and np.var(top_pvalues_perm)==0):
        for row in table:
            row['empirical_feature_p_value'] = row['p_value']
            row.update()
        alpha_para=-9
        beta_para=-9
    else:
        correction_function, alpha_para, beta_para = qtl_fdr_utilities.define_correction_function(top_pvalues_perm,cis_mode)
        for row in table:
            row['empirical_feature_p_value'] = correction_function(row['p_value'])
            row.update()
    table.flush()
    return [alpha_para, beta_para]

if(bestPermutationPval[perm] > min(relevantOutput)):
  bestPermutationPval[perm] = min(relevantOutput)

alpha_para, beta_para = output_writer.apply_pval_correction(feature_id.replace("/","-"),bestPermutationPval, cis_mode)
