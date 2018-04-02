from cobra import __version__ as cobra_version

def pyMADE(cobra_model, fold_change, pvals, gene_names, obj_frac, weighting, bounds, objs, p_thresh, p_eps,
           transition_matrix, remove_rev, theoretical_match, log_fold_change, return_models, verbose, set_IntFeasTol,
           weight_thresh, verify, round_states):
    """MADE: Metabolic Adjustment by Differential Expression.



    :param cobra_model: An object of cobra.Model class. pyMADE also accepts a list of models for each condition
    :param fold_change: Measured fold change from expression data.  Columns correspond to conditions, rows correspond to genes.
    :param pvals: P-values for changes.  Format is the same as for fold_change.
    :param gene_names:
    :param obj_frac:
    :param weighting:
    :param bounds:
    :param objs:
    :param p_thresh:
    :param p_eps:
    :param transition_matrix:
    :param remove_rev:
    :param theoretical_match:
    :param log_fold_change:
    :param return_models:
    :param verbose:
    :param set_IntFeasTol:
    :param weight_thresh:
    :param verify:
    :param round_states:
    :return:
    """
