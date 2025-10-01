# Chebifier

Chebifier is a tool for automated classification of chemicals in the [ChEBI](https://www.ebi.ac.uk/chebi/) ontology.

## News
- 25/10/01: Fixed issue where server crashed if when running predict without adding a SMILES string.

## Installation

### Setup Backend

Some dependencies require that pytorch is already installed:

`pip install torch`

After that, you can install the prediction system and web framework:

`pip install -r backend/requirements.txt`

*Chebifier* comes with a number of mandatory configuration files. `config.template.json` contains a template for a *Chebifier* configuration. Copy the contents of this file 

`cp backend/config.template.json backend/config.json`

and change the path for each setting according to your setup.

 * ELECTRA_CHECKPOINT : Path to a chebai-electra checkpoint,
 * BATCH_SIZE: Number of molecules that are passed to the model at once,
 * CLASS_HEADERS: Mapping of prediction labels to ChEBI classes,
 * CHEBI_JSON: Hierarchy of labels defined in `CLASS_HEADERS` (this file can be generate with [robot export](http://robot.obolibrary.org/export) using the options `--header "ID|LABEL|SubClasses" --entity-format ID`)



### Setup Frontend

Change to the respective directory and build the node.js files
```
cd react-app
npm run build
```

### Run in development

You can now start the development server with 

```
cd backend
flask run
```

The server should now run at [localhost:5000](localhost:5000)

## What's New?

## Citation

If you found Chebifier useful, please cite: 
[Martin Glauer, Fabian Neuhaus, Simon Flügel, Marie Wosny, Till Mossakowski, Adel Memariani, Johannes Schwerdt and Janna Hastings "Chebifier: Automating Semantic Classification in ChEBI to Accelerate Data-driven Discovery."Digital Discovery, 2024, 3, 896.](https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00238a)

