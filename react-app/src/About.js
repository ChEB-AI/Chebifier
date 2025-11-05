import Paper from '@mui/material/Paper';
import Typography from '@mui/material/Typography';
import Link from '@mui/material/Link';
import Box from '@mui/material/Box';
import axios from "axios";
import * as React from "react";
import {SlChemistry} from "react-icons/sl";

const About = () => {
    const [availableModels, setAvailableModels] = React.useState([]);
    const [availableModelsInfoTexts, setAvailableModelsInfoTexts] = React.useState([]);

    // Load once on mount so About content is fetched when the site loads
    React.useEffect(() => {
        axios.get('/api/modelinfo').then(response => {
            setAvailableModels(response.data.available_models || []);
            setAvailableModelsInfoTexts(response.data.available_models_info_texts || []);
        }).catch(() => {
            // silently ignore, page content still renders
        });
    }, []);

    const modelList = availableModels.map((model, index) => (
        <Box key={`model-${model}-${index}`} sx={{mb: 2}}>
            <Typography variant="h6" component="h3" gutterBottom>{model}</Typography>
            <Typography variant="body1" component="div" sx={{"& p": {margin: 0}}} dangerouslySetInnerHTML={{ __html: availableModelsInfoTexts[index] }} />
        </Box>
    ));

    return (
        <div className="App">
            <header className="App-header">
                <Box sx={{
                    width: '100%',
                    minHeight: '100vh',
                    backgroundColor: '#f7f2e7',
                    display: 'flex',
                    flexDirection: 'column',
                }}>
                    <Paper sx={{
                        width: '90%',
                        maxWidth: '900px',
                        mx: 'auto',
                        my: 2,
                        p: 2,
                        backgroundColor: '#ffffff',
                        borderRadius: 2,
                        boxShadow: 3
                    }}>
                        <Typography variant="h4" component="h2" gutterBottom>About</Typography>
                        <Typography variant="body1" paragraph>
                            Chebifier is a tool for automated classification of chemicals in the <Link href="https://www.ebi.ac.uk/chebi/">ChEBI</Link> ontology.
                            Currently, it can predict 1,700+ ChEBI classes.
                        </Typography>
                        <Typography variant="body1" paragraph>
                            To run a prediction, enter a SMILES string (or multiple ones, line-separated) or upload a file. Then,
                            hit the predict button (running the model might take a few seconds).
                            You can get more information about a result by clicking on it.
                        </Typography>

                        <Typography variant="h5" component="h3" gutterBottom><SlChemistry /> News</Typography>
                        <Typography variant="body1" paragraph>
                            <b>11/2025</b>: Added new models (Graph Attention Networks and augmented Graph Neural Networks).
                            Improved the Ensemble weighting mechanism. Redesigned the user interface.
                        </Typography>
                        <Typography variant="body1" paragraph>
                            <b>08/2025</b>: Added the ensemble. Added ChemLog, C3P and Graph Convolutional Networks.
                        </Typography>

                        <Typography variant="h5" component="h3" gutterBottom><SlChemistry /> The Ensemble</Typography>
                        <Typography variant="body1" paragraph>
                            Chebifier uses an ensemble of machine learning models and rule-based methods to classify molecules into ChEBI classes.
                            A weighting mechanism and inconsistency resolution are applied to ensure you get the best our models can over.
                            Details about the ensemble and the implementation can be found <Link href="https://github.com/ChEB-AI/python-chebifier">here</Link>.
                        </Typography>

                        <Typography variant="h5" component="h3" gutterBottom><SlChemistry /> Models</Typography>
                        <Typography variant="body1" paragraph>
                            At the moment, the following prediction models are supported by Chebifier.
                            You can either use them in the ensemble or select a specific model.
                        </Typography>

                        {/* Available models inside a Paper for emphasis */}
                        {/*<Paper sx={{p: 2, backgroundColor: '#ffffff', borderRadius: 2, boxShadow: 1}}>*/}
                        {modelList}

                        <Typography variant="h5" component="h3" sx={{mt: 3}} gutterBottom><SlChemistry /> Main Publication for Chebifier</Typography>
                        <Typography variant="body1">
                            Glauer, Martin, et al.: Chebifier: Automating Semantic Classification in ChEBI to Accelerate
                            Data-driven Discovery; Digital Discovery 3.5 (2024), <Link
                            href="https://doi.org/10.1039/D3DD00238A">Link</Link>
                        </Typography>
                    </Paper>
                </Box>
            </header>
        </div>

    );
};

export default About;
