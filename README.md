# MOMA

## Steps to run MOMA
### Step 1: Install software requirements

* R 3.4.0 or above
* Python 3.6.3 or above
* Python TensorFlow package
* Python numpy package
* Python pandas package

Tip: once you install python TensorFlow, you can simply install all other required python packages by
<code>pip install -r requirements.txt</code>.

### Step 2: Download the dataset
Then you download the Ecomcis transcriptome data from [here](https://www.dropbox.com/sh/t3zs3jbmq1efj3q/AAATQNlJimWT1bnTI9uK81S9a?dl=0) and place it in Dataset folder.

### Step 3: Preprocess the dataset
This step will preprocess the original dataset in the format that MOMA can train. For this, type
```Rscript preprocess_dataset.R Dataset/Ecomics.transcriptome.no_avg.v8.txt Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt```

Note that the code reads information from <code>Dataset/Meta.txt</code>, <code>Dataset/Meta.Medium.txt, <code>Dataset/Meta.Strain.txt</code>. This will save the preprocessed dataset in the file <code>Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt</code>

### Step 4: Run MOMA
Then you can run MOMA (prediction of transcriptomic response from characteristics of experimental condition) by

```python3 run_moma.py Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt Dataset/GRN.txt OPTIMIZATION_METHOD CONDITION_INDEX_TO_TEST```

* <code>Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt</code> is the dataset to be used for leave-one-condition-out cross-validation.
* <code>Dataset/GRN.txt</code> is the list of gene-regulatory relations (gene-regulatory network). This information is used to regulate the recurrent weight matrix. That is, we constrain the weight matrix not to have nonzero weights on any connections between genes that are not in the gene-regulatory network.
* <code>OPTIMIZATION_METHOD</code> can be SGD or RMSProp. To speed up the model training, RMSprop is recommended for <code>OPTIMIZATION_METHOD</code>.
* <code>CONDITION_INDEX_TO_TEST</code> is the index of condition (that is, a row index, ranging from 0 to 492 as there are 493 conditions or 493 rows in the <code>Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt</code>) to test its prediction from the model that is built on the rest of conditions (Leave-one-condition-out cross-validation; refer to [Kim et al. Nature commms 2016](https://www.nature.com/articles/ncomms13090) for more information).

Please note the following before running the model:
* Note that some of test conditions will not produce prediction results if the conditions are cross-validatable (for example, strain of the test condition is JM109 but this strain is not in the training data).
* The prediction results will display in the console in PCC metric (that is, PCC between predicted expression levels and known expression levels for the test condition) in comparison to the wildtype baseline (that is, PCC between mean expression levels of wildtype profiles and known expression levels for the test condition).
