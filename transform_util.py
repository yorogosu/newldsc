import scipy as sp
from scipy.linalg import pinv, cholesky

def transform_genotypes(norm_genotypes, grm, normalize=True):
    """
    This function performs genotype whitening by use of the Cholesky root of the inverse of a genetic relationship matrix.
    """    
    if normalize:
        norm_genotypes = norm_genotypes.T
        norm_genotypes = (norm_genotypes - sp.mean(norm_genotypes, 0)) / sp.std(norm_genotypes, 0, ddof=1)
        norm_genotypes = norm_genotypes.T
    
    transformation_matrix = cholesky(pinv(grm))
    transformed_genotypes = sp.dot(norm_genotypes, transformation_matrix)
    return transformed_genotypes

