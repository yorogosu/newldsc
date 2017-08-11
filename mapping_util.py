import pandas as pd
import scipy as sp
import scipy.interpolate

def bp_to_cM(target_snp_positions, chromosome_str, genet_map_dir='/home/yorgos/HeritPartition/faststorage/genetic_maps'):
    """
    This function takes a reference genetic map and calculates genetic positions (in cM) for a set of target SNP positions
    (in bp) through linear interpolation. IMPORTANT: interpolation only works for scipy v0.17.0 or later!!!
    """    
    genet_map_file = genet_map_dir+'/genetic_map_'+str(chromosome_str)+'_combined_b37.txt' #these file names are pretty standard
    g_map = pd.read_table(genet_map_file, names=['position', 'rate', 'cM'], skiprows=1, delim_whitespace=True)
    
    try:
        interpolation = sp.interpolate.interp1d(g_map.position, g_map.cM, bounds_error=False, fill_value=(0, max(g_map.cM)))
    except ValueError:
            print 'Interpolation only works for scipy v0.17.0 or later!!!'
    
    interpolated = interpolation(target_snp_positions)
    return interpolated

