import * as React from 'react';
import axios from "axios";
import Alert from '@mui/material/Alert';
import Box from '@mui/material/Box';
import Paper from '@mui/material/Paper';
import Divider from '@mui/material/Divider';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import TextField from '@mui/material/TextField';
import Chip from '@mui/material/Chip';
import LightbulbIcon from '@mui/icons-material/Lightbulb';
import StartIcon from '@mui/icons-material/Start';
import Modal from '@mui/material/Modal';
import FormLabel from '@mui/material/FormLabel';
import FormControl from '@mui/material/FormControl';
import Tooltip from '@mui/material/Tooltip';
import { randomId } from '@mui/x-data-grid-generator';

import DetailsPage from "./details-page";
import {plot_ontology} from "./ontology-utils";
import {CircularProgress} from "@mui/material";
import Select from '@mui/material/Select';
import MenuItem from '@mui/material/MenuItem';



export default function ClassificationGrid() {
    const [rows, setRows] = React.useState([]);
    const [detail, setDetail] = React.useState(null);

    const [ontology, setOntology] = React.useState(null);

    const [availableModels, setAvailableModels] = React.useState([]);
	const [availableModelsInfoTexts, setAvailableModelsInfoTexts] = React.useState([]);
    const [selectedModel, setSelectedModel] = React.useState('Ensemble');
    const [modelsLoaded, setModelsLoaded] = React.useState(false);

    const [inputText, setInputText] = React.useState("");
    const [predictionsLoading, setPredictionsLoading] = React.useState(false);

    const buildSelectedModels = () => {
        // Backend expects an object map of modelName -> boolean
        if (selectedModel === 'Ensemble') {
            const allTrue = {};
            availableModels.forEach(m => { allTrue[m] = true; });
            return allTrue;
        } else {
            const map = {};
            availableModels.forEach(m => { map[m] = (m === selectedModel); });
            return map;
        }
    };

    const selectedModels = buildSelectedModels();


    if (availableModels.length === 0) {
        axios.get('/api/modelinfo').then(response => {
            setAvailableModels(response.data.available_models);
            setAvailableModelsInfoTexts(response.data.available_models_info_texts);
            setModelsLoaded(true);

        });
    }

    const renderClasses = (params) => {
        const data = params.value;
        const violations = params.row.violations;
        if (data == null){
          return  <Alert severity="error">Could not process input!</Alert>
        } else {
            // todo the violation handling is too much logic for the front end - this should be done in the backend
            // (or be removed entirely as our ensemble - by design - resolves violations before outputting them)
            return (
            	<Box sx={{ display: 'flex', flexWrap: 'wrap' }}>
                {data.map((x) => {
                	var isViolation = false;
                	var tooltipText = "";
                	for (var i = 0; i < violations.length; i++) {
                		const this_violation = violations[i][0];
                        // each violation has the structure {violated_cls_id: (violated_cls_label, [(child1_id, child1_label), (child2_id, child2_label), ...])}
                		var children_str = [];
						for (const [violated_cls, violated_cls_label_and_children] of Object.entries(this_violation)) {
                            const violated_cls_label = violated_cls_label_and_children[0];
                            const children = violated_cls_label_and_children[1];
							children_str.push(`${children.map(child => child[1]).join(', ')} ${(children.length !== 1) ? "are subclasses" : "is a subclass"} of ${violated_cls_label}`);
						}
                        const violated_cls_labels = Object.values(this_violation).map(v => v[0]);
                        tooltipText = `This prediction is inconsistent: ${violated_cls_labels.join(' and ')} are marked as disjoint in ChEBI. ${children_str.join(', ')}`;
                		for (const [violated_cls, violated_cls_label_and_children] of Object.entries(this_violation)) {
                            const children_ids = violated_cls_label_and_children[1].map(child => child[0]);
                			if (children_ids.includes(x[0])) {
                				isViolation = true;
                                break;
                			}
                		}

					}
                	if (isViolation) {
						return (
							<Tooltip key={x[0]} title={tooltipText} placement="top" arrow>
								<Chip
									component="a"
									href={"http://purl.obolibrary.org/obo/" + x[0].replace(":", "_")}
									label={x[1]}
									clickable
									target="_blank"
									sx={{ backgroundColor: isViolation ? 'red' : 'default' }}
								/>
							</Tooltip>
						);
                	} else {
                		return (
                			<Chip component="a" href={"http://purl.obolibrary.org/obo/" + x[0].replace(":", "_")} label={x[1]} clickable target="_blank"/>
                		);
                	}
                })}
            	</Box>
			);
        }
    };

    const [open, setOpen] = React.useState(false);
    const [detailsLoading, setDetailsLoading] = React.useState(false);
    const handleOpen = (id) => () => {
        const thisRow = rows.find((row) => row.id === id);
        setDetailsLoading(thisRow.id);
        axios.post('/api/details', {smiles: thisRow.smiles, selectedModels: buildSelectedModels()}).then(response => {
            setDetail({
                plain_molecule: response.data.figures.plain_molecule,
                models_info: response.data.models,
                chebi: response.data.classification,
                chebi_legend: response.data.color_legend
            });
            setOpen(true);
        }).finally(() => {
            setDetailsLoading(false);
        });

    }
    const handleClose = () => setOpen(false);

    return (
      <div className="App">
        <header className="App-header">
        <Box sx={{width: "100%"}}>

            <Box sx={{ padding: 2, backgroundColor: '#f0f0f0', marginBottom: 2, borderRadius: 1, marginLeft: 2, marginRight: 2 }}>
                <Typography variant="h6" align="left" color="textPrimary" gutterBottom>
                    If you like Chebifier, please cite: Glauer, Martin, et al. "Chebifier: Automating Semantic
                    Classification in ChEBI to Accelerate Data-driven Discovery."
                    <a href={"https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00238a"}>Digital Discovery, 2024, 3, 896</a>.
                </Typography>
            </Box>

            <Paper sx={{width: "100%"}}>
                <Box
                    sx={{
                        height: 600,
                        width: '100%',
                        '& .actions': {
                            color: 'text.secondary',
                        },
                        '& .textPrimary': {
                            color: 'text.primary',
                        },
                    }}
                >
                    <Box sx={{ p: 2 }}>
                        <TextField
                          label="Enter SMILES (one per line)"
                          placeholder="C1=CC=CC=C1\nCC(=O)O"
                          value={inputText}
                          onChange={(e) => setInputText(e.target.value)}
                          onKeyDown={(e) => {
                            if (e.key === 'Enter' && !e.shiftKey) {
                              e.preventDefault();
                              if (!modelsLoaded || predictionsLoading) return;
                              const smiles = inputText.trim().replace(/\r/g, '').split('\n').map(s => s.trim()).filter(Boolean);
                              if (smiles.length === 0) return;
                              // initialize rows and run classification
                              const ids = smiles.map(() => randomId());
                              setRows(smiles.map((s, i) => ({ id: ids[i], smiles: s, direct_parents: [], predicted_parents: [] })));
                              setPredictionsLoading(true);
                              axios({
                                url: '/api/classify',
                                method: 'post',
                                data: {
                                  smiles: smiles,
                                  ontology: true,
                                  selectedModels: selectedModels
                                }
                              }).then(response => {
                                setRows((old) => old.map((row, i) => ({
                                  ...row,
                                  direct_parents: response.data.direct_parents[i],
                                  predicted_parents: response.data.predicted_parents[i],
                                  violations: response.data.violations[i]
                                })));
                                setOntology(response.data.ontology);
                              }).finally(() => setPredictionsLoading(false));
                            }
                          }}
                          fullWidth
                          multiline
                          minRows={3}
                        />
                        <Box sx={{ mt: 2, display: 'flex', alignItems: 'flex-end', gap: 2, flexWrap: 'wrap' }}>
                          <FormControl sx={{ minWidth: 240 }}>
                            <Select
                              value={selectedModel}
                              onChange={(e) => setSelectedModel(e.target.value)}
                              disabled={!modelsLoaded}
                            >
                              <MenuItem value={'Ensemble'}>Ensemble</MenuItem>
                              {availableModels.map((m, idx) => (
                                <MenuItem key={m} value={m}>{m}</MenuItem>
                              ))}
                            </Select>
                          </FormControl>
                          <Button
                            variant="contained"
                            onClick={() => {
                              const smiles = inputText.trim().replace(/\r/g, '').split('\n').map(s => s.trim()).filter(Boolean);
                              if (smiles.length === 0) return;
                              const ids = smiles.map(() => randomId());
                              setRows(smiles.map((s, i) => ({ id: ids[i], smiles: s, direct_parents: [], predicted_parents: [] })));
                              setPredictionsLoading(true);
                              axios({
                                url: '/api/classify',
                                method: 'post',
                                data: {
                                  smiles: smiles,
                                  ontology: true,
                                  selectedModels: selectedModels
                                }
                              }).then(response => {
                                setRows((old) => old.map((row, i) => ({
                                  ...row,
                                  direct_parents: response.data.direct_parents[i],
                                  predicted_parents: response.data.predicted_parents[i],
                                  violations: response.data.violations[i]
                                })));
                                setOntology(response.data.ontology);
                              }).finally(() => setPredictionsLoading(false));
                            }}
                            disabled={predictionsLoading || !modelsLoaded}
                            startIcon={predictionsLoading ? <CircularProgress size={20}/> : <StartIcon/>}
                          >
                            Predict
                          </Button>
                        </Box>
                        <Divider sx={{ my: 2 }} />
                        <Box>
                          {rows.length > 0 && rows.map((row) => (
                            <Paper key={row.id} sx={{ p: 2, mb: 1 }}>
                              <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                                <Typography variant="subtitle2">{row.smiles}</Typography>
                                <Button
                                  size="small"
                                  onClick={handleOpen(row.id)}
                                  disabled={!(row.direct_parents && row.direct_parents.length > 0) || predictionsLoading}
                                  startIcon={detailsLoading === row.id ? <CircularProgress size={16}/> : <LightbulbIcon/>}
                                >
                                  Details
                                </Button>
                              </Box>
                              <Box sx={{ mt: 1 }}>
                                {renderClasses({ value: row.direct_parents, row })}
                              </Box>
                            </Paper>
                          ))}
                        </Box>
                    </Box>

                </Box>
            </Paper>

            <Paper>
                {plot_ontology(ontology,true,false)}
            </Paper>




            <Modal
              open={open}
              onClose={handleClose}
              aria-labelledby="modal-modal-title"
              aria-describedby="modal-modal-description"
            >
                <Box sx={{
                  mb: 2,
                  display: "flex",
                  flexDirection: "column",
                  width: '95%',
                  height: '95%',
                  position: 'fixed',
                  left: '2.5%',
                  top: '2.5%',
                }}>

                    <DetailsPage detail={detail} handleClose={handleClose}/>
                </Box>
            </Modal>

        </Box>
        </header>
  </div>
    );
}
