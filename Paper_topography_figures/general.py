# HELPERS

from scipy.stats import mannwhitneyu, wilcoxon, kruskal, ttest_1samp
# Stats
def get_star(p):
    '''
    Helper that returns significance stars or "ns" (non-significant) string
    depending on p value 
    '''
    # Create string representation of significance for figure legend
    if   p < 0.001:
        sig = '***'
    elif p < 0.01:
        sig = '**'
    elif p < 0.05:
        sig = '*'
    else:
        sig = 'ns'
    return sig 


def print_mannwhitneyu(values_A, 
                       values_B, 
                       label_A='A',
                       label_B='B', 
                       alternative='two-sided'
                       ):
    ''' Convenience function to calculate and print results of Mann Whitney U test'''
    
    u, p = mannwhitneyu(values_A, values_B, alternative=alternative)
    if p < 0.0001:
        pstr = f'{p:.2e}'
    else:
        pstr = f'{p:.5f}'

    sig = get_star(p)

    if len(values_A) == len(values_B):
        print(f'Mann-Whitney U={u:.3f}, n{label_A}=n{label_B}={len(values_B)}, p={pstr}{sig} {alternative}')
    else: 
        print(f'Mann-Whitney U={u:.3f}, n{label_A}={len(values_A)}, n{label_B}={len(values_B)}, p={pstr}{sig} {alternative}')
    return u,p

def print_wilcoxon(values, 
                   label='A',
                   alternative='two-sided'
                   ):
    ''' 
    Convenience function to calculate and print results of 
    Wilcoxon signed-rank test
    '''
    
    z, p = wilcoxon(values, alternative=alternative)

    if p < 0.0001:
        pstr = f'{p:.2e}'
    else:
        pstr = f'{p:.5f}'

    sig = get_star(p)
    print(f'Wilcoxon signed-rank test Z={z:.3f}, n{label}={len(values)}, p={pstr}{sig} {alternative}')
    return z,p


def print_ttest_1samp(values, 
                      popmean,
                      label='A',
                      ):
    ''' 
    Convenience function to calculate and print results of 
    (two-sided) 1 sample ttest 
    '''
    
    t, p = ttest_1samp(values, popmean)

    if p < 0.0001:
        pstr = f'{p:.2e}'
    else:
        pstr = f'{p:.5f}'

    sig = get_star(p)
    print(f'T-test against population mean of {popmean:.1f} t={t:.3f}, n{label}={len(values)}, p={pstr}{sig} two-sided')
    return t,p


def print_kruskal(values):
    ''' 
    Convenience function to calculate and print results of
    Kruskal-Wallis H-test for independent samples
    ''' 

    args_kruskal = [l for l in values]
    h, p = kruskal(*args_kruskal)

    if p < 0.0001:
        pstr = f'{p:.2e}'
    else:
        pstr = f'{p:.5f}'

    sig = get_star(p)

    print(f'Kruskal-Wallis H={h:.3f}, n={len(values)}, p={pstr}{sig}')

    return h, p 