import React from 'react';
import axios from 'axios'

import {useEffect, useRef} from "react";
import {Network} from "vis-network";

import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import Box from '@mui/material/Box';
import Button from '@mui/material/Button';
import CancelIcon from '@mui/icons-material/Close';
import Divider from '@mui/material/Divider';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import FormControl from '@mui/material/FormControl';
import Grid from '@mui/material/Grid';
import Paper from '@mui/material/Paper';
import Skeleton from '@mui/material/Skeleton';
import Stack from '@mui/material/Stack';
import Tab from '@mui/material/Tab';
import TabContext from '@mui/lab/TabContext';
import TabList from '@mui/lab/TabList';
import TabPanel from '@mui/lab/TabPanel';
import TextField from '@mui/material/TextField';
import Typography from '@mui/material/Typography';
import {styled} from '@mui/material/styles';

import {plot_ontology} from "./ontology-utils";
import {DetailsElectra} from "./details-electra";
import {DetailsChemlog} from "./details-page-chemlog";




// is this used anywhere?
const Item = styled(Paper)(({theme}) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.textsecondary,
}));


const Legend = ({ colors }) => {
	return (
		<Box>
			{Object.entries(colors).map(([color, text]) => (
				<Box key={color} sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
					<Box sx={{ width: 30, height: 10, backgroundColor: color, borderRadius: 1 }} />
					<Typography>{text}</Typography>
				</Box>
			))}
		</Box>
	);
};

export function DetailsPerModel(data) {
	const models_info = data.models_info;
	return (
		<Box>
			{Object.entries(models_info).map(([model_name, model_data]) => (
				<Box>
					<h2>Insights for {model_name}:</h2>
					{model_data.model_type === "ElectraPredictor" ? (
						<DetailsElectra model_data={model_data} />
					) : model_data.model_type === "ChemLogPredictor" ? (
						<DetailsChemlog model_data={model_data} />
					) : (
						<Typography>Model type {model_data.model_type} is not supported for explanations.</Typography>
					)}
				</Box>
			))
			}
		</Box>
	);
};

export default function DetailsPage(data) {
    const handleClose = data.handleClose;
    data = data.detail;

    return (
        <Box sx={{ height: '100%'}}>
            <Box sx={{ height: '100%'}}>
                <Paper sx={{ height: '100%'}}>
                    <Button color="primary" onClick={handleClose} startIcon={<CancelIcon/>}/>
                    <Box sx={{
                        overflowX: "scroll",
                        overflowY: "scroll",
                        overflow: "auto",
                        display: "flex",
                        flexDirection: "column",
                        position: 'relative',
                        width: '94%',
                        height: '94%',
                        top: '1%',
                        left: '1%'
                    }}>
                        <Box sx={{ height: '100%'}}>
                            <Box>
                                <Typography><h2>Molecular structure</h2></Typography>
                                <div style={{
                                    display: 'flex',
                                    alignItems: 'left',
                                    justifyContent: 'center',
                                }}>
                                <img src={`data:image/jpeg;base64,${data.plain_molecule}`} width="300" height="300"/></div>
                            </Box>
                            <Box>
                                <Typography><h2>ChEBI Classification</h2></Typography>
                                <Typography><h3>What am I seeing?</h3>
                                    Chebifier uses an <i>Ensemble</i> of Machine Learning (ML) and rule-based methods for
                                    each class. The ensemble prediction is influenced by several factors:
                                    <ol>
                                        <li>Trust: We are testing our ML models on a validation set. If they achieve a high F1-score for a given class, we trust them more.
                                            The rule-based methods have a higher trust than the ML models.</li>
                                        <li>Confidence: The ML models can report themselves if they are sure of their prediction for a given class and sample.
                                            If a model is very sure of its prediction, we will give it more weight.</li>
                                        <li>Logical consistency: ChEBI has logical relations between classes that we can use.
                                            For instance, <a target="_blank" rel="noopener noreferrer" href={"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:25699"}>organic ion</a> is a subclass of <a target="_blank" rel="noopener noreferrer" href={"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:50860"}>organic molecular entity</a>.
                                            Therefore, if the ensemble predicts organic ion, but not organic molecular entity, this will be corrected. Similarly, disjoint relations
                                            (e.g., a molecule cannot be an organic cation and an organic anion) are resolved automatically.</li>
                                    </ol>
                                    In the graph below, the color indicates by which
                                    method each prediction has been made.</Typography>
                                <Typography><b>Legend</b></Typography>
                                <Legend colors={data.chebi_legend} />
                                <Box>
                                    {plot_ontology(data.chebi)}
                                </Box>
                            </Box>
                            {DetailsPerModel(data)}
                        </Box>
                    </Box>
                </Paper>
            </Box>
        </Box>
    )
}


