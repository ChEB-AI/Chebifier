import * as React from 'react';
import { useEffect } from 'react';
import {Text} from 'react-native';
import axios from "axios";
import Alert from '@mui/material/Alert';
import PropTypes from 'prop-types';
import Box from '@mui/material/Box';
import Paper from '@mui/material/Paper';
import Divider from '@mui/material/Divider';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import Switch from '@mui/material/Switch';
import AddIcon from '@mui/icons-material/Add';
import EditIcon from '@mui/icons-material/Edit';
import CloseIcon from '@mui/icons-material/Close';
import Chip from '@mui/material/Chip';
import DeleteIcon from '@mui/icons-material/DeleteOutlined';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import LightbulbIcon from '@mui/icons-material/Lightbulb';
import SaveIcon from '@mui/icons-material/Save';
import CancelIcon from '@mui/icons-material/Close';
import StartIcon from '@mui/icons-material/Start';
import Modal from '@mui/material/Modal';
import FormLabel from '@mui/material/FormLabel';
import FormControl from '@mui/material/FormControl';
import FormGroup from '@mui/material/FormGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormHelperText from '@mui/material/FormHelperText';
import Tooltip from '@mui/material/Tooltip';
import Link from '@mui/material/Link';

import {
    GridRowModes,
    DataGrid,
    GridToolbarContainer,
    GridActionsCellItem,
} from '@mui/x-data-grid';
import {
    randomId,
} from '@mui/x-data-grid-generator';

import DetailsPage from "./details-page";
import {plot_ontology} from "./ontology-utils";
import {CircularProgress} from "@mui/material";

const RenderDate = (props) => {
  const { hasFocus, value } = props;
  const buttonElement = React.useRef(null);
  const rippleRef = React.useRef(null);

  React.useLayoutEffect(() => {
    if (hasFocus) {
      const input = buttonElement.current?.querySelector('input');
      input?.focus();
    } else if (rippleRef.current) {
      // Only available in @mui/material v5.4.1 or later
      rippleRef.current.stop({});
    }
  }, [hasFocus]);

  return (
    <strong>
      {value?.getFullYear() ?? ''}
      <Button
        component="button"
        ref={buttonElement}
        touchRippleRef={rippleRef}
        variant="contained"
        size="small"
        style={{ marginLeft: 16 }}
        // Remove button from tab sequence when cell does not have focus
        tabIndex={hasFocus ? 0 : -1}
        onKeyDown={(event) => {
          if (event.key === ' ') {
            // Prevent key navigation when focus is on button
            event.stopPropagation();
          }
        }}
      >
        Open
      </Button>
    </strong>
  );
};

const Checkbox = ({label, value, onChange, checked=true}) => {
	return (
		<label key={value}>
			<input type="checkbox" name={label} checked={checked} onChange={onChange} />
			{label}
		</label>
	);
};



function EditToolbar(props) {
    const {setRows, setRowModesModel, rows, getLabel, setOntology, selectedModels} = props;

    const [predictionsLoading, setPredictionsLoading] = React.useState(false);

    const addRows = ((smiles) => {
            const ids = smiles.map((s) => randomId());
            setRows((oldRows) => [...oldRows, ...smiles.map((s, i) => ({id: ids[i], smiles: s, direct_parents: [], predicted_parents: []}))]);
            return ids
        }
    )

    const handleAdd = () => {
        const ids = addRows([''])
        setRowModesModel((oldModel) => ({
            ...oldModel,
            ...Object.fromEntries(ids.map(id => [id, {mode: GridRowModes.Edit, fieldToFocus: 'smiles'}])),
        }));
    };

    const handleUpload = (event) => {
        event.preventDefault();
        const reader = new FileReader()
        reader.onload = async (e) => {
            addRows(e.target.result.trim().replace(/\r/g, "").split("\n"))
        };
        reader.readAsText(event.target.files[0])
    };

    const handleDownload = (event) => {
        event.preventDefault();
        const fileData = JSON.stringify(rows.map((r) => ({"smiles": r["smiles"], "direct_parents": r["direct_parents"],"predicted_parents": r["predicted_parents"],})).filter((d) => d.direct_parents?.length >= 0));
        const blob = new Blob([fileData], { type: "text/plain" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.download = "chebifier-predictions.json";
        link.href = url;
        link.click();
    };

    const handleRun = () => {
        setPredictionsLoading(true);
        axios({
            url: '/api/classify',
            method: 'post',
            data: {
            	smiles: rows.map((r) => (r["smiles"])),
            	ontology: true,
            	selectedModels: selectedModels
            }
        }).then(response => {
            setRows((oldRows) => oldRows.map((row, i) => ({
                ...row, "direct_parents": response.data.direct_parents[i], "predicted_parents": response.data.predicted_parents[i], "violations": response.data.violations[i]
            })));
            setOntology(response.data.ontology)
        }).finally(() => {
            setPredictionsLoading(false);
        });
    };

    return (
        <GridToolbarContainer>
            <Button color="primary" startIcon={<AddIcon/>} onClick={handleAdd}>
                Add SMILES
            </Button>
            <Button color="primary" startIcon={<FileUploadIcon/>} component="label">
                Upload file
                <input
                    accept="text/plain"
                    style={{display: 'none'}}
                    id="file-upload"
                    type="file"
                    onChange={handleUpload}
                />
            </Button>
            <Divider/>
            <Button color="primary" startIcon={predictionsLoading ? <CircularProgress size={20}/> : <StartIcon/>} onClick={handleRun} disabled={predictionsLoading}>
                Predict classes
            </Button>
            <Button color="primary" startIcon={<FileDownloadIcon/>} onClick={handleDownload} disabled={rows.filter((d) => d.direct_parents?.length > 0).length === 0}>
                Download JSON
            </Button>
        </GridToolbarContainer>
    );
}

export default function ClassificationGrid() {
    const [rows, setRows] = React.useState([]);
    const [rowModesModel, setRowModesModel] = React.useState({});
    const [detail, setDetail] = React.useState(null);
    const [hierarchy, setHierarchy] = React.useState({});

    const [ontology, setOntology] = React.useState(null);

    const [availableModels, setAvailableModels] = React.useState([]);
	const [availableModelsInfoTexts, setAvailableModelsInfoTexts] = React.useState([]);
    const [selectedModels, setSelectedModels] = React.useState({});



    if (Object.keys(hierarchy).length === 0) {
        axios.get('/api/hierarchy').then(response => {
            setHierarchy(response.data.hierarchy);
            setAvailableModels(response.data.available_models);
            setAvailableModelsInfoTexts(response.data.available_models_info_texts);
           	var newSelectedModels = {};
           	for (var i = 0; i < response.data.available_models.length; i++) {
           		newSelectedModels[response.data.available_models[i]] = true;
           	}
			setSelectedModels(newSelectedModels);


        });
    }

    const handleDeleteClick = (id) => () => {
        setRows(rows.filter((row) => row.id !== id));
    };

    const processRowUpdate = (newRow) => {
        const oldRow = rows.find(row => (row.id === newRow.id))
        if (oldRow.smiles != newRow.smiles) {
            newRow.predicted_parents = [];
            newRow.direct_parents = [];
            newRow.isNew = false;
        }
        setRows(rows.map((row) => (row.id === newRow.id ? newRow : row)));
        return newRow;
    };

    const renderClasses = (params) => {
        const data = params.value;
        const violations = params.row.violations;
        if (data == null){
          return  <Alert severity="error">Could not process input!</Alert>
        } else {
            return (
            	<Box sx={{ display: 'flex', flexWrap: 'wrap' }}>
                {data.map((x) => {
                	var isViolation = false;
                	var tooltipText = "";
                	for (var i = 0; i < violations.length; i++) {
                		const this_violation = violations[i][0];
                		var children_str = [];
						for (const [violated_cls, children] of Object.entries(this_violation)) {
							children_str.push(`${children.map(v => hierarchy[v].label).join(', ')} ${(children.length !== 1) ? "are subclasses" : "is a subclass"} of ${hierarchy[violated_cls].label}`);
						}
                		for (const [violated_cls, children] of Object.entries(this_violation)) {
                			if (children.includes(x)) {
                				isViolation = true;

								tooltipText = `This prediction is inconsistent: ${Object.keys(this_violation).map(v => hierarchy[v].label).join(' and ')} are marked as disjoint in ChEBI. ${children_str.join(', ')}`;
								break;
                			}
                		}

					}
                	if (isViolation) {
						return (
							<Tooltip key={x} title={tooltipText} placement="top" arrow>
								<Chip
									component="a"
									href={"http://purl.obolibrary.org/obo/" + x.replace(":", "_")}
									label={hierarchy[x].label}
									clickable
									target="_blank"
									sx={{ backgroundColor: isViolation ? 'red' : 'default' }}
								/>
							</Tooltip>
						);
                	} else {
                		return (
                			<Chip component="a" href={"http://purl.obolibrary.org/obo/" + x.replace(":", "_")} label={hierarchy[x].label} clickable target="_blank"/>
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
        axios.post('/api/details', {smiles: thisRow.smiles, selectedModels: selectedModels}).then(response => {
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

    const handleCheckboxChange= (event) => {
    	const checkedModel = event.target.name;
    	setSelectedModels({...selectedModels, [checkedModel]: event.target.checked});
	};

	const modelList = availableModels.map((model,index) => (
		<>
			<Tooltip title={availableModelsInfoTexts[index]} placement="bottom-start" arrow>
			<FormControlLabel
				control={<Switch checked={selectedModels[model]} onChange={handleCheckboxChange} name={model} />}
				label={model}
			/>
			</Tooltip>
		</>

    ));




    const columns = [
        {
            field: 'smiles',
            headerName: 'SMILES',
            flex: 0.45,
            editable: true,
            preProcessEditCellProps: (params) => {
                if (params.hasChanged) {
                    const newRow = {...params.row, "predicted_parents": [], "direct_parents": [], "isNew": false, "smiles": params.props.value}
                    setRows(rows.map((row) => (row.id === newRow.id ? newRow : row)));
                }
                return { ...params.props,};
            },
        },
        {field: 'direct_parents', headerName: 'Predicted Class', flex: 0.45, editable: false, renderCell:renderClasses},
        {
            field: 'actions',
            type: 'actions',
            headerName: 'Actions',
            flex: 0.1,
            cellClassName: 'actions',
            getActions: ({id}) => {
                const isInEditMode = rowModesModel[id]?.mode === GridRowModes.Edit;
                const thisRow = rows.find((row) => row.id === id);
                const wasPredicted = thisRow.direct_parents?.length > 0;

                return [
                	<GridActionsCellItem
                        icon={detailsLoading === id ? <CircularProgress size={20}/> : <LightbulbIcon/>}
                        label="Details"
                        onClick={handleOpen(id)}
                        color="inherit"
                        disabled={!wasPredicted || detailsLoading !== false}
                    />,
                    <GridActionsCellItem
                        icon={<DeleteIcon/>}
                        label="Delete"
                        onClick={handleDeleteClick(id)}
                        color="inherit"
                    />,
                ];
            },
        },
    ];

    const getLabel = (x) => {
        return hierarchy[x]["label"]
    }

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
                    <DataGrid
                        rows={rows}
                        columns={columns}
                        editMode="row"
                        rowModesModel={rowModesModel}
                        onRowModesModelChange={(newModel) => setRowModesModel(newModel)}
                        processRowUpdate={processRowUpdate}
                        getRowHeight={() => 'auto'}
                        components={{
                            Toolbar: EditToolbar,
                        }}
                        componentsProps={{
                            toolbar: {setRows, setRowModesModel, rows, getLabel, setOntology, selectedModels},
                        }}
                        experimentalFeatures={{newEditingApi: true}}
                    />

                </Box>
            </Paper>

            <Paper>
                {plot_ontology(ontology,true,false)}
            </Paper>
			<Paper>
				<FormControl component="fieldset" variant="standard">
					<FormLabel component="legend">Select models:</FormLabel>
					<FormGroup>
						{modelList}
					</FormGroup>
				</FormControl>
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
