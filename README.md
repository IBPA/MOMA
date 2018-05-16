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
Then you download the Ecomcis transcriptome data from [here](https://www.dropbox.com/s/85iifwwc0eevfti/ecomics.transcriptome.no_avg.v8.txt?dl=0) and place it in Dataset folder.

### Step 3: Preprocess the dataset
This step will preprocess the original dataset in the format that MOMA can train. For this, WILDTYPE
<code>Rscript preprocess_dataset.R Dataset/Ecomics.transcriptome.no_avg.v8.txt Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt</code>

This will save the preprocessed dataset in the file <code>Dataset/Ecomics.transcriptome_with_meta.avg.v8.txt</code>

### Step 4: Run MOMA
Then you can 
