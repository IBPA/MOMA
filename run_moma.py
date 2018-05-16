#!/usr/bin/python

import sys
import pandas as pd
import numpy  as np

import keras.backend  as K
from keras.constraints import Constraint
from keras.models      import Sequential
from keras.layers      import SimpleRNN
from keras.callbacks   import EarlyStopping
from keras.optimizers  import RMSprop, SGD

### CONSTANTS ###

OUTPUT_COLUMN_NAME_PREFIX    = "m"       # e.g. m.b0001
STRAIN_COLUMN_NAME_PREFIX    = "strain"  # e.g. strain.b3943
MEDIUM_COLUMN_NAME_PREFIX    = "medium"  # e.g. medium.Glu
STRESS_COLUMN_NAME_PREFIX    = "stress"  # e.g. stress.O2-starvation
GP_COLUMN_NAME_PREFIX        = "gp"      # e.g. gp.b1039_KO
PREFIX_DELIMETER             = "."
WT_CONDITIONS                = ['MG1655.MD001.none.na_WT','MG1655.MD018.none.na_WT'] # WILDTYPE CONDITION WITH M9+Glu or LB

### check if each of {strain,medium,stress,genetic_perturbations} of test_condition is in train_conditions. 
### return True only when there are all of four information in train_conditions.
def is_condition_validatable(test_condition, train_conditions):
    list_test_condition = test_condition.split(PREFIX_DELIMETER)
    np_train_conditions = np.matrix(list(map(lambda x: x.split(PREFIX_DELIMETER), train_conditions)))
    for i in range(len(list_test_condition)):
        if [list_test_condition[i]] not in np_train_conditions[:,i]:
            return False
    return True

def get_wildtype_exp(conditions, outputs):
    wt_conditions = list(map(lambda x: x in WT_CONDITIONS, conditions))
    wt_profiles   = outputs[wt_conditions].mean()
    return wt_profiles

def load_gene_regulatory_network(grn_file, column_names):
    pd_grn        = pd.read_csv(grn_file,sep="\t")
    np_grn_matrix = np.zeros((len(column_names), len(column_names)))
    for idx, row in pd_grn.iterrows():
        source_gene, target_gene = map(lambda x: OUTPUT_COLUMN_NAME_PREFIX + PREFIX_DELIMETER + x, row)
        if source_gene in column_names and target_gene in column_names:
            np_grn_matrix[column_names.index(source_gene), column_names.index(target_gene)] = 1
            np_grn_matrix[column_names.index(target_gene), column_names.index(source_gene)] = 1
    return np_grn_matrix

class GeneRegulatoryConstraint(Constraint):
    def __init__(self, gene_regulatory_network):
        self._gene_regulatory_network = gene_regulatory_network

    def __call__(self, w):
        return K.tf.multiply(w, self._gene_regulatory_network)

### DATA READING ###
pd_data                    = pd.read_csv(sys.argv[1],sep="\t")
column_names               = pd_data.columns.tolist()
output_column_names        = list(filter(lambda x: x.startswith(OUTPUT_COLUMN_NAME_PREFIX + PREFIX_DELIMETER), column_names))
input_strain_column_names  = list(filter(lambda x: x.startswith(STRAIN_COLUMN_NAME_PREFIX + PREFIX_DELIMETER), column_names))
input_medium_column_names  = list(filter(lambda x: x.startswith(MEDIUM_COLUMN_NAME_PREFIX + PREFIX_DELIMETER), column_names))
input_stress_column_names  = list(filter(lambda x: x.startswith(STRESS_COLUMN_NAME_PREFIX + PREFIX_DELIMETER), column_names))
input_gp_column_names      = list(filter(lambda x: x.startswith(GP_COLUMN_NAME_PREFIX     + PREFIX_DELIMETER), column_names))
input_column_names         = input_strain_column_names + input_stress_column_names + input_medium_column_names + input_gp_column_names
n_conditions               = pd_data.shape[0]
n_features                 = len(input_column_names)
n_timesteps                = 4
X                          = np.zeros((n_conditions, n_timesteps, n_features))
Y                          = pd_data[output_column_names].as_matrix()
gene_regulatory_network    = load_gene_regulatory_network(sys.argv[2], output_column_names)

earlystop_on               = True
n_hidden_nodes             = len(output_column_names)
n_epochs                   = 1000
hidden_node_activation     = 'sigmoid'
recurrent_dropout          = 0.7
optimization_method        = eval(sys.argv[3]) # RMSProp or SGD
condition_idx_to_test      = int(sys.argv[4])
note                       = sys.argv[5]

print("[DATA READ] {} samples, {} timesteps, {} features".format(n_conditions, n_timesteps, n_features))

#### DATA TRANSFORMATION (2D to 3D tensor) ###
for idx, row in pd_data.iterrows():
    for timestep in range(n_timesteps):
        X[idx, timestep, :] = row[input_column_names]
    
### CREATE MODEL ###
model = Sequential()
model.add(SimpleRNN(n_hidden_nodes,
            return_sequences=False,
            activation=hidden_node_activation,
            recurrent_constraint=GeneRegulatoryConstraint(gene_regulatory_network),
            recurrent_dropout=recurrent_dropout,
            input_shape=(n_timesteps, n_features)))
print("[MODEL CREATED]")
print(model.summary())

### MODEL COMPILE
earlystop  = EarlyStopping(monitor='loss', min_delta=0.0001, patience=4, verbose=0, mode='auto')
callbacks  = [earlystop] if earlystop_on else []
rmsprop    = optimization_method()
loss       = 'mse'
model.compile(loss=loss, optimizer=rmsprop, metrics=["mse"])

### LEAVE-ONE-CONDITION-OUT CROSS-VALIDATION ###
pd_train_data = pd_data.drop(condition_idx_to_test)
if not is_condition_validatable(row['Cond'], pd_train_data['Cond']):
    print("[MODEL EVALUATION {}] {} is skipped because it is not validatable".format(condition_idx_to_test, row['Cond']))
    exit()

train_X = np.delete(X, condition_idx_to_test, 0)
train_Y = np.delete(Y, condition_idx_to_test, 0)
test_X  = X[condition_idx_to_test,:,:]
test_Y  = Y[condition_idx_to_test,:]

### set initial state of recurrent hidden nodes to mean expression level of wildtype profiles
wildtype_state            = get_wildtype_exp(pd_train_data['Cond'], pd_train_data[output_column_names])
model.layers[0].states[0] = K.variable(value=np.array(wildtype_state))
history                   = model.fit(train_X, train_Y, epochs=n_epochs, verbose=1, callbacks=callbacks, batch_size=32)   

### EVALUATE MODEL PERFORMANCE ###
test_predicted_Y  = model.predict(test_X[np.newaxis,:,:])
test_pcc          = np.corrcoef(test_predicted_Y, test_Y)[0,1]
wt_baseline_pcc   = np.corrcoef(wildtype_state, test_Y)[0,1]

print("[MODEL EVALUATION {}] {} PCC: {}, WT Baseline: {} (note: {})".format(condition_idx_to_test, pd_data['Cond'][condition_idx_to_test], test_pcc, wt_baseline_pcc, note))


