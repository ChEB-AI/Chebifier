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


const GLOBAL_MOL_PARAMS = {
	width: 500,
	height: 500,
}
window
      .initRDKitModule()
      .then(function (RDKit) {
        console.log("RDKit version: " + RDKit.version());
        window.RDKit = RDKit;
        /**
         * The RDKit module is now loaded.
         * You can use it anywhere.
         */
      })
      .catch(() => {
        // handle loading errors here...
      });


const NetworkElement = (data) => {
    const visJsRef = useRef(null);
    useEffect(() => {
        if(data.graph != null){
            const g = {
                nodes: data.graph.nodes,
                edges: data.graph.edges
            }
            const network = visJsRef.current && new Network(
                visJsRef.current, g, {
                physics: false,
                width: "100%",
                height: "100%",
                clickToUse: true
            });
            network.fit()
        }
    }, [visJsRef, data]);
    return <Box><div ref={visJsRef}/></Box>
}

const LayerComponent = (data) => {
    const [value, setValue] = React.useState(0);

    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    return (
        <Box sx={{width: '100%', typography: 'body1'}}>
            <TabContext value={value}>
                <Box sx={{borderBottom: 1, borderColor: 'divider'}}>
                    <TabList onChange={handleChange} aria-label="lab API tabs example">
                        {data.layer.highlights.map((highlight, i) => <Tab label={data.layer.name + " " + (i + 1)} value={i}/>)}
                    </TabList>
                </Box>
                {data.layer.highlights.map((g, i) => <TabPanel value={i}><div className="svg-mol" dangerouslySetInnerHTML={{__html: data.mol.get_svg_with_highlights(JSON.stringify({...GLOBAL_MOL_PARAMS, 'atoms': g}))}}></div></TabPanel>)}
            </TabContext>
        </Box>
    );
}


export function LayerTabs(data) {
    const [value, setValue] = React.useState(0);
    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    return (
        <Box sx={{width: '100%', typography: 'body1'}}>
            <TabContext value={value}>
                <Box sx={{borderBottom: 1, borderColor: 'divider'}}>
                    <TabList onChange={handleChange} aria-label="lab API tabs example">
                        {data.layers.map((layer, i) => <Tab label={layer.name} value={i} centered/>)}
                    </TabList>
                </Box>
                {data.layers.map((layer, i) => <TabPanel value={i}><LayerComponent mol={data.mol} layer={layer}/></TabPanel>)}
            </TabContext>
        </Box>
    );
}


export function HighlightsBlocks(data) {
	const [value, setValue] = React.useState(0);
    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    var blocks = data.highlights;
    var blocks_content = [];
    for (let i = 0; i < blocks.length; i++) {
    	var block_type = blocks[i][0];
    	var block_content = blocks[i][1];
    	console.log(block_type)
    	console.log(block_content)
    	if (block_type == "text") {
			blocks_content.push(
				<Box>
					<Typography>{block_content}</Typography>
				</Box>
			);
		} else if (block_type == "single") {
			var mdetails = {...GLOBAL_MOL_PARAMS};
			mdetails["atoms"] = block_content;
			var svg_mol = data.mol.get_svg_with_highlights(JSON.stringify(mdetails));
			blocks_content.push(
				<Box>
					<div className="svg-mol" dangerouslySetInnerHTML={{__html: svg_mol}}></div>
				</Box>
			);

		} else if (block_type == "tabs") {
			var layers = [];
			for (const[key, value] of Object.entries(block_content)) {
				layers.push({name: key, highlights: value});
			}
			console.log(layers)
			blocks_content.push(
				<LayerTabs layers={layers} mol={data.mol}/>
			);
		}
	}
	return blocks_content;

}

const Item = styled(Paper)(({theme}) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.textsecondary,
}));

export default function DetailsPage(data) {
    const handleClose = data.handleClose;
    data = data.detail;
    console.log(data)

    var smiles = data.smiles
  	var mol = window.RDKit.get_mol(smiles);
  	var svg_mol = mol.get_svg_with_highlights(JSON.stringify(GLOBAL_MOL_PARAMS));
  	svg_mol = svg_mol.substring(svg_mol.indexOf("<svg"));

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
								<Box>
									<div className="svg-mol" dangerouslySetInnerHTML={{__html: svg_mol}}></div>
								</Box>
                            </Box>

                            <Box>
                                <Typography><h2>ChemLog Classification</h2></Typography>
                                <Typography><h3>What am I seeing?</h3>
                                    Highlights!</Typography>
                                <HighlightsBlocks highlights={data.highlights} mol={mol}/>
                            </Box>
                        </Box>
                    </Box>
                </Paper>
            </Box>
        </Box>
    )
}


