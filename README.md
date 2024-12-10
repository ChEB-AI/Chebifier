# Installation

## Setup Backend

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
 * CHEBI_JSON: Hierarchy of labels defined in `CLASS_HEADERS` (this file can be generate with [robot export](http://robot.obolibrary.org/export) using `--input chebi.obo` and ``--header ID|LABEL|SubClasses|``)



## Setup Frontend

Change to the respective directory and build the node.js files
```
cd react-app
npm run build
```

## Run in development

You can now start the development server with 

```
cd backend
flask run
```

The server should now run at [localhost:5000](localhost:5000)
