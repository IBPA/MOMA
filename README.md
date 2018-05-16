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

This will save the preprocessed dataset in the file <code>Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt</code>

### Step 4: Run MOMA
Then you can run MOMA (prediction of transcriptomic response from characteristics of experimental condition) by

```python3 run_moma.py Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt Dataset/GRN.txt OPTIMIZATION_METHOD CONDITION_INDEX_TO_TEST```

where OPTIMIZATION_METHOD can be SGD or RMSProp and CONDITION_INDEX is the index of condition (ranging from 0 to 492 as there are 493 conditions) to test its prediction from the model that is built on the rest of conditions (Leave-one-condition-out cross-validation; refer to Kim et al. Nature commms 2016 for more information). To speed up the model training, RMSprop is recommended for OPTIMIZATION_METHOD. Note that some of test conditions will not produce prediction results if the conditions are cross-validatable (for example, strain of the test condition is JM109 but this strain is not in the training data). 
