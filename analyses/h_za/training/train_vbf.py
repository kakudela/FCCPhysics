import uproot
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
import ROOT
import pickle
import os


ROOT.gROOT.SetBatch(True)
# e.g. https://root.cern/doc/master/tmva101__Training_8py.html

def load_process(proc, variables, target=0):
    weight = 1
    if ":" in proc:
        tmp = proc.split(":")
        proc = tmp[0]
        weight = float(tmp[1])

    fIn = f"{inputDir}/{proc}.root"

    f = uproot.open(fIn)
    print(f.keys())
    tree = f["events"]
    xsec = f["crossSection"].value
    ev_proc = f["eventsProcessed"].value
    weight *= xsec*3e6/ev_proc
    weight = 1.0
    #weight = 0.2*10.8e6 / ev_proc # 200 fb-1 for each Higgs decay == equal

    
    print(f"Load {fIn.replace('.root', '')} with {tree.num_entries} events and weight {weight} and cross-section {xsec}")

    df = tree.arrays(variables, library="pd") # convert the signal and background data to pandas DataFrames
    #df = df.head(int(0.1*ev_proc))
    df['target'] = target # add a target column to indicate signal (1) and background (0)
    df['weight'] = weight
    return df





################################################################################

tag = "vbf" # hqqa hvva

inputDir = "/ceph/submit/data/group/fcc/ee/analyses/h_za/treemaker/ecm365/vbf/"

variables = ["qq_m", "qq_p", "qqa_p", "qq_trans", "vv_trans", "vv_long", "qq_long", "dr_vv_qqa", "dr_vv_qq", "dr_a_qq", "dr_a_vv", "cos_qqa", "cos_qq", "jet1_p", "jet2_p", "jet1_theta", "jet2_theta", "photon_p", "photon_theta"]

sig_procs = ['wzp6_ee_nuenueH_HZa_ecm365']
bkg_procs = ['p8_ee_ZZ_ecm365']

parms = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'n_estimators': 2000,
    'max_depth': 6,
    #'eval_metric': ["merror", "mlogloss"],
    'early_stopping_rounds': 5,
}


################################################################################


print("Parse inputs")

dfs = []





for sig in sig_procs:
    df = load_process(sig, variables, target=1)
    dfs.append(df)



for bkg in bkg_procs:
    df = load_process(bkg, variables)
    dfs.append(df)



# Concatenate the dataframes into a single dataframe
data = pd.concat(dfs, ignore_index=True)


# split data in train/test events
train_data, test_data, train_labels, test_labels, train_weights, test_weights  = train_test_split(
    data[variables], data['target'], data['weight'], test_size=0.3, random_state=42
)



# conversion to numpy needed to have default feature_names (fN), needed for conversion to TMVA
train_data = train_data.to_numpy()
test_data = test_data.to_numpy()
train_labels = train_labels.to_numpy()
test_labels = test_labels.to_numpy()
train_weights = train_weights.to_numpy()
test_weights = test_weights.to_numpy()



params = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'eta': 0.1,
    'max_depth': 3,
    'subsample': 0.5,
    'colsample_bytree': 0.5,
    'seed': 42,
    'n_estimators': 350, # low number for testing purposes (default 350)
    'early_stopping_rounds': 1,
    'num_rounds': 20,
    'learning_rate': 0.20,
    'gamma': 3,
    'min_child_weight': 10,
    'max_delta_step': 0,
}

parms_opt = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'eta': 0.1,
    'max_depth': 5,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'seed': 42,
    'n_estimators': 500, # low number for testing purposes (default 350)
    'early_stopping_rounds': 10,
    'num_rounds': 20,
    'learning_rate': 0.05,
    'gamma': 0.1,
    'min_child_weight': 10,
    'max_delta_step': 0,
    'reg_alpha': 0.1,
    'reg_lambda': 1.0,
    'use_label_encoder': False,
}

parms_softprob = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'n_estimators': 2000,
    'max_depth': 8,
    #'eval_metric': ["merror", "mlogloss"],
    'early_stopping_rounds': 10,
}



# train the XGBoost model
print("Start training")
#eval_set = [(train_data, train_labels), (test_data, test_labels)]
eval_set = [(train_data, train_labels), (test_data, test_labels)]
bdt = xgb.XGBClassifier(**parms_opt)
#bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight=train_weights)
bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight_eval_set=[train_weights, test_weights])



# export model (to ROOT and pkl)
print("Export model")
fOutName = f"/ceph/submit/data/group/fcc/ee/analyses/za/training/{tag}.root"
ROOT.TMVA.Experimental.SaveXGBoost(bdt, tag, fOutName, num_inputs=len(variables))

# append the variables
variables_ = ROOT.TList()
for var in variables:
     variables_.Add(ROOT.TObjString(var))
fOut = ROOT.TFile(fOutName, "UPDATE")
fOut.WriteObject(variables_, "variables")


save = {}
save['model'] = bdt
save['train_data'] = train_data
save['test_data'] = test_data
save['train_labels'] = train_labels
save['test_labels'] = test_labels
save['variables'] = variables
pickle.dump(save, open(f"{fOutName.replace('.root', '.pkl')}", "wb"))

print(tag)
print(f"Save to {fOutName}")

os.system(f"python FCCPhysics/analyses/h_za/training/evaluate.py --tag {tag}")