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
import { SlChemistry } from "react-icons/sl";
// import Modal from '@mui/material/Modal';
import FormLabel from '@mui/material/FormLabel';
import FormControl from '@mui/material/FormControl';
import Tooltip from '@mui/material/Tooltip';
import {randomId} from '@mui/x-data-grid-generator';

import DetailsPage from "./details-page";
import {plot_ontology, MoleculeStructure} from "./ontology-utils";
import {Molecules} from "./ontology-utils";
import {CircularProgress} from "@mui/material";
import Select from '@mui/material/Select';
import MenuItem from '@mui/material/MenuItem';
import Collapse from '@mui/material/Collapse';
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';

export default function ClassificationGrid() {
  const [rows, setRows] = React.useState([]);
  const [detailsByRow, setDetailsByRow] = React.useState({});

  const [availableModels, setAvailableModels] = React.useState([]);
  const [availableModelsInfoTexts, setAvailableModelsInfoTexts] = React.useState([]);
  const [selectedModel, setSelectedModel] = React.useState('Ensemble');
  const [modelsLoaded, setModelsLoaded] = React.useState(false);

  const [inputText, setInputText] = React.useState("");
  const [predictionsLoading, setPredictionsLoading] = React.useState(false);
  const [hasPredicted, setHasPredicted] = React.useState(false);
  const [expandedRowId, setExpandedRowId] = React.useState(null);
  // map of `${rowId}-${classIdx}` -> boolean for per-chip "Why this class?" panel
  const [openWhyMap, setOpenWhyMap] = React.useState({});
  const toggleWhy = (rowId, classIdx) => () => {
    const key = `${rowId}-${classIdx}`;
    setOpenWhyMap(prev => ({...prev, [key]: !prev[key]}));
  };
  const isWhyOpen = (rowId, classIdx) => !!openWhyMap[`${rowId}-${classIdx}`];

  const buildSelectedModels = () => {
    // Backend expects an object map of modelName -> boolean
    if (selectedModel === 'Ensemble') {
      const allTrue = {};
      availableModels.forEach(m => {
        allTrue[m] = true;
      });
      return allTrue;
    } else {
      const map = {};
      availableModels.forEach(m => {
        map[m] = (m === selectedModel);
      });
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
    const row = params.row || {};
    const isExpanded = row.id === expandedRowId;

    if (data === null) {
      return <Alert severity="error">Could not process input!</Alert>
    }

    // Collapsed view: chips flowing inline
    const collapsedView = (
      <Box sx={{display: 'flex', flexWrap: 'wrap', gap: 1}}>
        {data.map((x, idx) => (
          <Box key={`class-collapsed-${row.id}-${idx}`} sx={{display: 'flex', alignItems: 'center', gap: 1}}>
            <Chip
              component="a"
              href={`http://purl.obolibrary.org/obo/CHEBI_${x[0]}`}
              label={x[1]}
              clickable
              target="_blank"
            />
          </Box>
        ))}
      </Box>
    );

    // Expanded view: one class per row with a Why button and collapsible table
    const expandedView = (
      <Box sx={{display: 'flex', flexDirection: 'column', gap: 1}}>
        {data.map((x, idx) => {
          const whyOpen = isWhyOpen(row.id, idx);
          return (
            <Box key={`class-expanded-${row.id}-${idx}`} sx={{borderRadius: 1, border: '1px solid #eee'}}>
              <Box sx={{display: 'flex', alignItems: 'center', gap: 1, p: 1}}>
                <Chip
                  component="a"
                  href={`http://purl.obolibrary.org/obo/CHEBI_${x[0]}`}
                  label={x[1]}
                  clickable
                  target="_blank"
                />
                <Button size="small" onClick={toggleWhy(row.id, idx)} sx={{ml: 'auto'}}>
                  {whyOpen ? 'Hide' : 'Why this class?'}
                </Button>
              </Box>
              <Collapse in={whyOpen} timeout="auto" unmountOnExit>
                <Box sx={{px: 1, pb: 1, overflowX: 'auto'}}>
                  <Table size="small" aria-label="why-this-class" sx={{ minWidth: 700 }}>
                    <TableHead>
                      <TableRow>
                        <TableCell>Model</TableCell>
                        <TableCell>Prediction</TableCell>
                        <TableCell>Confidence</TableCell>
                        <TableCell>Trust</TableCell>
                        <TableCell>Model weight</TableCell>
                        <TableCell>Model score</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {Object.entries(x[2] || {}).map(([k, v]) => (
                        <TableRow key={`why-${row.id}-${idx}-${k}`} style={{backgroundColor: v?.model_score > 0 ? '#a8f099' : '#f0a699'}}>
                          <TableCell>{k}</TableCell>
                          <TableCell>{String(v?.prediction)}</TableCell>
                          <TableCell>{typeof v?.confidence === 'number' ? v.confidence.toFixed(3) : String(v?.confidence)}</TableCell>
                          <TableCell>{typeof v?.trust === 'number' ? v.trust.toFixed(3) : String(v?.trust)}</TableCell>
                          <TableCell>{typeof v?.model_weight === 'number' ? v.model_weight.toFixed(3) : String(v?.model_weight)}</TableCell>
                          <TableCell>{typeof v?.model_score === 'number' ? v.model_score.toFixed(3) : String(v?.model_score)}</TableCell>
                        </TableRow>
                      ))}
                      <TableRow>
                        <TableCell><b>Ensemble</b></TableCell>
                        <TableCell><b>true</b></TableCell>
                        <TableCell></TableCell>
                        <TableCell></TableCell>
                        <TableCell></TableCell>
                        <TableCell><b>{typeof x[3] === 'number' ? x[3].toFixed(3) : String(x[3])}</b></TableCell>
                      </TableRow>
                    </TableBody>
                  </Table>
                </Box>
              </Collapse>
            </Box>
          );
        })}
      </Box>
    );

    return (
      <Box>
        <Collapse in={!isExpanded} timeout="auto" unmountOnExit>
          {collapsedView}
        </Collapse>
        <Collapse in={isExpanded} timeout="auto" unmountOnExit>
          {expandedView}
        </Collapse>
      </Box>
    );
  };

  const [detailsLoading, setDetailsLoading] = React.useState(false);
  const handleToggleExpand = (id) => () => {
    if (expandedRowId === id) {
      setExpandedRowId(null);
      return;
    }
    setExpandedRowId(id);

    // fetch details if not cached
    if (!detailsByRow[id]) {
      const thisRow = rows.find((row) => row.id === id);
      if (!thisRow) return;
      setDetailsLoading(id);
      axios.post('/api/details', {smiles: thisRow.smiles, selectedModels: buildSelectedModels()}).then(response => {
        const detailObj = {
          models_info: response.data.models,
          chebi: response.data.classification,
          chebi_legend: response.data.color_legend
        };
        setDetailsByRow(prev => ({...prev, [id]: detailObj}));
      }).finally(() => {
        setDetailsLoading(false);
      });
    }
  }

  return (
    <div className="App">
      <header className="App-header">
        {(() => {
          const isCentered = !hasPredicted && rows.length === 0;
          return (
            <Box sx={{
              width: '100%',
              minHeight: '100vh',
              backgroundColor: '#f7f2e7',
              display: 'flex',
              flexDirection: 'column'
            }}>
              <Box sx={{
                padding: 2,
                backgroundColor: '#f0f0f0',
                marginBottom: 2,
                borderRadius: 1,
                marginLeft: 2,
                marginRight: 2
              }}>
                <Typography variant="h6" align="left" color="textPrimary" gutterBottom>
                  If you like Chebifier, please cite: Glauer, Martin, et al. "Chebifier: Automating Semantic
                  Classification in ChEBI to Accelerate Data-driven Discovery."
                  <a href={"https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00238a"}>Digital Discovery, 2024, 3,
                    896</a>.
                </Typography>
              </Box>

              <Paper sx={{
                width: 'fit-content',
                height: 'fit-content',
                backgroundColor: '#ffffff',
                boxShadow: 3,
                borderRadius: 2,
                flex: '0 0 auto',
                display: isCentered ? 'flex' : 'block',
                alignItems: isCentered ? 'center' : 'stretch',
                justifyContent: isCentered ? 'center' : 'flex-start',
                minHeight: 'auto',
                marginX: 'auto'
              }}>
                <Box sx={{width: 'auto', height: 'auto'}}>
                  <Box sx={{p: 2, width: 'auto', minWidth: '700px', display: 'inline-flex', flexDirection: 'column'}}>
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
                          setRows(smiles.map((s, i) => ({
                            id: ids[i],
                            smiles: s,
                            direct_parents: [],
                            predicted_parents: []
                          })));
                          setHasPredicted(true);
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
                              ontology: response.data.ontology[0][i],
                            })));
                          }).finally(() => setPredictionsLoading(false));
                        }
                      }}
                      fullWidth
                      multiline
                      minRows={3}
                    />
                    <Box sx={{mt: 2, display: 'flex', alignItems: 'flex-end', gap: 2, flexWrap: 'wrap', width: '100%'}}>
                      <FormControl>
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
                        sx={{ml: 'auto'}}
                        onClick={() => {
                          const smiles = inputText.trim().replace(/\r/g, '').split('\n').map(s => s.trim()).filter(Boolean);
                          if (smiles.length === 0) return;
                          const ids = smiles.map(() => randomId());
                          setRows(smiles.map((s, i) => ({
                            id: ids[i],
                            smiles: s,
                            direct_parents: [],
                            predicted_parents: []
                          })));
                          setHasPredicted(true);
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
                              ontology: response.data.ontology[0][i],
                            })));
                          }).finally(() => setPredictionsLoading(false));
                        }}
                        disabled={predictionsLoading || !modelsLoaded}
                        startIcon={predictionsLoading ? <CircularProgress size={20}/> : <SlChemistry/>}
                      >
                        Predict
                      </Button>
                    </Box>
                  </Box>

                </Box>
              </Paper>

              {hasPredicted && (
                <Paper sx={{
                  mt: 2,
                  width: '90%',
                  mx: 'auto',
                  p: 2,
                  backgroundColor: '#ffffff',
                  borderRadius: 2,
                  boxShadow: 1,
                  overflowX: 'auto'
                }}>
                  {rows.length > 0 && (
                    <Box>
                      {rows.map((row) => {
                        const canExpand = (row.direct_parents && row.direct_parents.length > 0) && !predictionsLoading;
                        return (
                        <Paper key={row.id} sx={{p: 2, mb: 1}}
                               onClick={(e) => {
                                 if (!canExpand) return;
                                 const t = e.target;
                                 const interactive = t.closest(
                                   'a, button, [role="button"], input, textarea, select, .MuiButtonBase-root, .MuiLink-root, .MuiChip-root'
                                 );
                                 if (interactive) return;
                                 const nestedPaper = t.closest('.MuiPaper-root');
                                 if (nestedPaper && nestedPaper !== e.currentTarget) return;
                                 handleToggleExpand(row.id)();
                               }}
                        >
                          <Box sx={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
                            <Typography variant="subtitle2">{row.smiles}</Typography>
                            <Button
                              size="small"
                              onClick={(e) => {e.stopPropagation(); handleToggleExpand(row.id)(); }}
                              disabled={!canExpand}
                              startIcon={detailsLoading === row.id && expandedRowId !== row.id ? <CircularProgress size={16}/> : <LightbulbIcon/>}
                            >
                              {expandedRowId === row.id ? 'Hide' : 'Details'}
                            </Button>
                          </Box>
                          {expandedRowId !== row.id && (
                            <Box sx={{mt: 1}}>
                              {renderClasses({value: row.direct_parents, row})}
                            </Box>
                          )}
                          {expandedRowId === row.id && (
                            <Box sx={{mt: 2, display: 'flex', flexWrap: 'wrap', gap: 2}}>
                              <Paper sx={{p: 2, flex: '1 1 320px', minWidth: 400}}>
                                <Typography variant="subtitle2" gutterBottom>Predicted classes</Typography>
                                {renderClasses({value: row.direct_parents, row})}
                              </Paper>
                              <Paper sx={{p: 2, flex: '1 1 480px', minWidth: 200, overflowX: 'auto'}}>
                                <Typography variant="subtitle2" gutterBottom>Molecular graph</Typography>
                                <MoleculeStructure smiles={row.smiles} height={200} width={200} />
                              </Paper>
                              <Paper sx={{p: 2, flex: '1 1 480px', minWidth: 820, overflowX: 'auto'}}>
                                <Typography variant="subtitle2" gutterBottom>Ontology graph</Typography>
                                {plot_ontology(row.ontology, true, false)}
                              </Paper>
                              <Paper sx={{p: 2, flex: '1 1 480px', minWidth: 420, overflow: 'hidden'}}>
                                <Typography variant="subtitle2" gutterBottom>Model-specific insights</Typography>
                                {detailsByRow[row.id] ? (
                                  <DetailsPage detail={detailsByRow[row.id]}
                                               handleClose={() => setExpandedRowId(null)}/>
                                ) : (
                                  <Box sx={{display: 'flex', alignItems: 'center', gap: 1}}>
                                    <CircularProgress size={20}/>
                                    <Typography variant="body2">Loading model-specific insights…</Typography>
                                  </Box>
                                )}
                              </Paper>

                            </Box>
                          )}
                        </Paper>
                        );
                      })}
                      <Divider sx={{my: 2}}/>
                    </Box>
                  )}
                </Paper>
              )}

            </Box>
          );
        })()}
      </header>
    </div>
  );
}
