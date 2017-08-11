import scipy as sp
from scipy.linalg import eigh

def calc_pca_score(norm_genotypes, eigenvecs, eigenvals, n_pcs, normalize=True, summing=True):
    """
    This function calculates PCA scores. Summing over k PC's is optional.
    """    
    if normalize:
        norm_genotypes = norm_genotypes.T
        norm_genotypes = (norm_genotypes - sp.mean(norm_genotypes, 0)) / sp.std(norm_genotypes, 0, ddof=1)
        norm_genotypes = norm_genotypes.T
    
    M,N = sp.shape(norm_genotypes)
    evecs = eigenvecs[:,:n_pcs]
    evecs = (evecs - sp.mean(evecs, 0)) / sp.std(evecs, 0, ddof=1)
    scaled_evals = eigenvals[:n_pcs] / sp.sum(eigenvals)
    lambdas = sp.diag(scaled_evals)
    
    snp_weights = sp.dot(norm_genotypes, evecs) / (N-1)
    pc_scores = snp_weights**2
    pc_scores = pc_scores - ((1-pc_scores)/(N-2))
    if summing:
        pc_scores = sp.sum(sp.dot(pc_scores, lambdas), 1)
    return pc_scores
    

def pca(grm):
    """
    Runs a normal PCA on a genetic relationship matrix but also sorts eigenvalues (and accordingly eigenvectors) from larger to smaller.
    """
    eigenval, eigenvec = eigh(sp.array(grm))
    sort_indices = sp.argsort(eigenval,)[::-1]
    eigenval = eigenval[sort_indices]
    eigenvec = eigenvec[sort_indices]
    return eigenval, eigenvec

