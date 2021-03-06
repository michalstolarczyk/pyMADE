# reading and generating the data for testing
import cobra
import numpy as np
import pandas as pd

cobra_model = cobra.io.read_sbml_model(
    "/home/mstolarczyk/Uczelnia/UVA/iPfal17_curation/output/curated_model_2018_02_27_13:04:25.xml")
df = pd.read_csv('/home/mstolarczyk/Uczelnia/UVA/pyMADE/testData.csv', header=None, index_col=0,
                 names=['gene_name', 'logFC', 'p_value'], skiprows=1)
logFC = df[['logFC']]
p_value = df[['p_value']]

fold_change = logFC
pvals = p_value
transition_matrix = None
obj_frac = 0.5
lb_list = []
ub_list = []
bounds = []
for r in cobra_model.reactions:
    lb_list.append(r.bounds[0])
    ub_list.append(r.bounds[1])
bounds.append(pd.DataFrame({'lb': lb_list, 'ub': ub_list}))
bounds.append(pd.DataFrame({'lb': lb_list, 'ub': ub_list}))
del lb_list, ub_list

objs = []
objs_inner = [rxn.objective_coefficient for rxn in cobra_model.reactions]
objs.append(objs_inner)
objs.append(objs_inner)

# The actual function
import cobra
import numpy as np
import pandas as pd


def pyMADE(cobra_model, fold_change, gene_names, pvals=None, obj_frac=0.3, weighting="log", bounds=None, objs=None,
           p_thresh=0.5, p_eps=0.0000000001, transition_matrix=None, remove_rev=False, theoretical_match=True,
           log_fold_change=False, return_models=True, verbose=True, set_IntFeasTol=0.0000000001,
           weight_thresh=0.00000001, verify=True, round_states=True):
    """MADE: Metabolic Adjustment by Differential Expression.


    :param cobra_model: An object of cobra.Model class. pyMADE also accepts a list of cobra_models for each condition
    :param fold_change: An object of pandas.DataFrame class or object that can be converted to one. Measured fold change from expression data.  Columns correspond to conditions, rows correspond to genes.
    :param pvals: An object of pandas.DataFrame class or object that can be converted to one. P-values for changes.  Format is the same as for fold_change.
    :param gene_names: Cell array of names for genes in expression dataset. These correspond to the rows in fold_change and pvals. If none is given, the rows correspond to cobra.Model.genes.
    :param obj_frac: Fraction of metabolic objective required in the resulting model (v_obj >= frac*v_obj_max). Default is 0.3.  Input can also be a list giving a separate fraction for each condition.
    :param weighting:  Method to convert PVALS to weights.  Options include:
                 'log'     w(p) = -log(p)    (default)
                 'linear'  w(p) = 1 - p
                 'unit'    w(p) = 1
                 'none'    No transformation -- PVALS are weights and
                             FOLD_CHANGE is a direction matrix
    :param bounds: List of objects of pandas.DataFrame class with two columns ('lb','up') and length and order corresponding to length and order of reactions in cobra_model. The number of objects must be equal the number of conditions. For example, bounds=[pd.DataFrame({'lb': lbi_list, 'ub': ubi_list})
,pd.DataFrame({'lb': lbj_list, 'ub': ubj_list})] are the lower and upper bounds for the ith and jth condition, respectively. If not specified, the bounds from cobra_model are copied to each condition
    :param objs: List of lists with condition-specific objectives. Each inner list must be of length equal to number of reactions and consist of 0 or 1 only, representing non-obective or objective reaction, respectively. The number of outer list elements (inner lists) must be equal the number of conditions. If not specified, the objectives coefficients from cobra.Reaction.obejctive_vcoefficient are copied.
 is copied to each condition.
    :param p_thresh: Threshold above which a P-value is not considered significant.  If a gene increases with a P-value p > P_THRESH, then it will be held constant with P-value 1 - p.  The default is 0.5.
    :param p_eps: P-values below P_EPS are considered equal to P_EPS.  This is used with log weighting to avoid taking the logarithm of very small P-values.  The default is 1e-10.
    :param transition_matrix:  A matrix (T) describing the interaction between conditions.  If T(i,j) = k, then the kth column in FOLD_CHANGE describes the change in expression between condition i to condition j.  If not given, FOLD_CHANGE assumes 1 -> 2, 2 -> 3, ..., n-1 -> n.
    :param remove_rev: Remove reversibility constraints from the model before converting to a MIP.  May improve performance if activity cycles are not a concern.  (Default = false)
    :param theoretical_match: Calculate the theoretical matches possible and adjust the match statistics.  (Default = true
    :param log_fold_change: If true, consider FOLD_CHANGE to be log fold change. (Default = false)
    :param return_models: If true, the models for each condition are returned. (Default = true)
    :param verbose: If true (default), a results table is printed.  Otherwise, pyMADE is silent
    :param set_IntFeasTol: A reduced value for the IntFeasTol solver parameter to avoid "integrality leaks" that return models with zero objective flux.  The default value is 1e-10.  If false, the parameter is not adjusted.  All solver parameters are reset to previous values before the function returns.
    :param weight_thresh: Threshold for weights on variables to be held constant.  MADE only attempts to keep variables constant if the weight from the P-value is above this value. (Default = 1e-8.)
    :param verify: If true (default), MADE checks that the generated models are feasible for the objective flux fraction.
    :param round_states: If true (default), binary variables in GENE_STATES are rounded
    :return:m
    """

    # Arguments processing
    if isinstance(transition_matrix, type(None)):  # check if transition_matrix was provided
        if isinstance(fold_change,
                      pd.DataFrame):  # check if the fold_change is an object of pandas.DataFrame class and convert it to on if needed
            ntrans = fold_change.shape[1]  # get the number of transitions
        else:
            try:
                fold_change = pd.DataFrame(fold_change)
            except ValueError:
                raise ValueError(
                    "The provided fold_change object cannot be converted to the object of pandas.DataFrame class.")
            ntrans = fold_change.shape

        ncond = ntrans + 1  # get the number of conditions
    else:
        try:
            transition_matrix = pd.DataFrame(transition_matrix)
            ntrans = transition_matrix.max
            ncond = transition_matrix.shape[0]
        except ValueError:
            raise ValueError(
                "The provided transition_matrix object cannot be converted to the object of pandas.DataFrame class.")
    if isinstance(pvals,
                  pd.DataFrame):  # check if the fold_change is an object of pandas.DataFrame class and convert it to on if needed
        assert pvals.shape == fold_change.shape, "fold_change and pavals must have the same dimensions."
    else:
        if isinstance(pvals, type(None)):
            pvals = pd.DataFrame(np.ones((fold_change.shape[0], fold_change.shape[1])))
            pvals.index = fold_change.index
            weighting = 'none'
            print("No P-values given; applying unit weighting.\n")
        try:
            pvals = pd.DataFrame(pvals)
        except ValueError:
            raise ValueError(
                "The provided pvals object cannot be converted to the object of pandas.DataFrame class.")

    print("The dataset includes " + str(ntrans) + " transitions between " + str(ncond) + " conditions.")

    if isinstance(cobra_model, list):  # check if multiple models were provided
        assert len(cobra_model) == ncond, "Number of models does not match number of conditions."
        models = cobra_model
    else:  # if just one - replicate it
        models = [cobra_model] * ncond

    for i in models:  # check if each input model is an object of cobra.Model class
        assert isinstance(i, cobra.Model), "All modes have to be objects of cobra.Model class."

    if isinstance(bounds, type(None)):
        pass  # keep the original models bounds
    else:
        assert len(bounds) == ncond, "Number of bounds holding objects must match the number of conditions."
        for i in range(len(models)):
            assert len(models[i].reactions) == bounds[i].shape[
                0], "Number of bounds within bounds holding object must match number of reactions in the model."
            for j in range(len(models[i].reactions)):
                models[i].reactions[j].bounds = tuple(bounds[i].iloc[j, :])

    if isinstance(objs, type(None)):
        pass  # keep the original models objectives
    else:
        assert len(objs) == ncond, "Number of objectives holding objects must match the number of conditions."
        for i in range(len(models)):
            assert len(models[i].reactions) == len(objs[i]), "Number of bounds within bounds holding object must match number of reactions in the model."
            for j in range(len(models[i].reactions)):
                models[i].reactions[j].objective_coefficient = objs[i][j]

    frac = obj_frac
    if len(obj_frac) == 1:  # set growth fraction for each condition if it was not specified
        frac = [obj_frac] * ncond
    assert len(frac) == ncond, "Number of objective fraction elements has to be equal the number of conditions."
