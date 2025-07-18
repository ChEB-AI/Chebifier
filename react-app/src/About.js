import Paper from '@mui/material/Paper';
import Text from '@mui/material/Typography';
import Link from '@mui/material/Link';
import Box from '@mui/material/Box';
import axios from "axios";
import * as React from "react";
import Tooltip from "@mui/material/Tooltip";
import FormControlLabel from "@mui/material/FormControlLabel";
import Switch from "@mui/material/Switch";

const About = () => {
    const [availableModels, setAvailableModels] = React.useState([]);
    const [availableModelsInfoTexts, setAvailableModelsInfoTexts] = React.useState([]);

    if (availableModels.length === 0) {
        axios.get('/api/hierarchy').then(response => {
            setAvailableModels(response.data.available_models);
            setAvailableModelsInfoTexts(response.data.available_models_info_texts);
        });
    }

    const modelList = availableModels.map((model, index) => (
        <Text>
			<h3>{model}</h3>
            <span dangerouslySetInnerHTML={{ __html: availableModelsInfoTexts[index]} } />
        </Text>

    ));


    return (
        <div className="App">
            <header className="App-header">
                <Paper sx={{width: "100%"}}>
                    <Box
                        sx={{
                            height: '100%',
                            width: '100%',
                            '& .actions': {
                                color: 'text.secondary',
                            },
                            '& .textPrimary': {
                                color: 'text.primary',
                            },
                            paddingBottom: 5,
                        }}
                    >
                        <h2>About</h2>
                        <Text>
						Chebifier is a tool for automated classification of chemicals in the <Link href="https://www.ebi.ac.uk/chebi/">ChEBI</Link> ontology.
                        <br/>
						To run a prediction, add a new row using "ADD SMILES" and enter a SMILES string or upload a file. Then,
						start the prediction with "PREDICT CLASSES". You can get additional insights by clicking on the light bulb icon on the right-hand side.
                        <br/>
                        </Text>
                        <h3>The model</h3>
                        <Text>
                            Chebifier uses an ensemble of machine learning models and rule-based methods to classify molecules into ChEBI classes. Details about the ensemble and the implementation can be found <Link href="https://github.com/ChEB-AI/python-chebifier">here</Link>.
						    At the moment, the following prediction models are supported by Chebifier. You can select specific models using the buttons at the bottom of the classification page.
						{modelList}
                        </Text>
                            <h3>Main Publication for Chebifier</h3>
                        <Text>
                            Glauer, Martin, et al.: Chebifier: Automating Semantic Classification in ChEBI to Accelerate
                            Data-driven Discovery; Digital Discovery 3.5 (2024), <Link
                            href="https://doi.org/10.1039/D3DD00238A">Link</Link>
                            <br/>
                        </Text>
                    </Box>
                </Paper>
            </header>
        </div>

    );
};

export default About;
