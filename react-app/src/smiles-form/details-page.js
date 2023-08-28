import React from 'react';
import axios from 'axios'

import {useEffect, useRef} from "react";
import {Network} from "vis-network";

import Box from '@mui/material/Box';
import Button from '@mui/material/Button';
import CancelIcon from '@mui/icons-material/Close';
import Paper from '@mui/material/Paper';
import Tab from '@mui/material/Tab';
import TabContext from '@mui/lab/TabContext';
import TabList from '@mui/lab/TabList';
import TabPanel from '@mui/lab/TabPanel';
import Typography from '@mui/material/Typography';
import {styled} from '@mui/material/styles';

import {plot_ontology} from "./ontology-utils";

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
                        {data.layer.map((layer, i) => <Tab label={"Head " + (i + 1)} value={i}/>)}
                    </TabList>
                </Box>
                {data.layer.map((g, i) => <TabPanel value={i}>< NetworkElement graph={g} /></TabPanel>)}
            </TabContext>
        </Box>
    );
}


export function LayerTabs(layers) {
    const [value, setValue] = React.useState(0);
    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    return (
        <Box sx={{width: '100%', typography: 'body1'}}>
            <TabContext value={value}>
                <Box sx={{borderBottom: 1, borderColor: 'divider'}}>
                    <TabList onChange={handleChange} aria-label="lab API tabs example">
                        {layers.layers.map((layer, i) => <Tab label={"Layer " + (i + 1)} value={i} centered/>)}
                    </TabList>
                </Box>
                {layers.layers.map((layer, i) => <TabPanel value={i}><LayerComponent layer={layer}/></TabPanel>)}
            </TabContext>
        </Box>
    );
}

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
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                }}>
                                <img src={`data:image/jpeg;base64,${data.plain_molecule}`} alt="Molecular graph of the submitted chemical" width="300" height="300"/></div>
                            </Box>
                            <Box>
                                <Typography><h2>ChEBI Classification</h2></Typography>
                                <Typography><h3>What am I seeing?</h3>
                                    The graph below shows you the position in the ChEBI ontology that our system
                                    proposed. Not that this prediction is an estimate based on available data and may be
                                    prone to errors.</Typography>
                                <Box>
                                    {plot_ontology(data.chebi)}
                                </Box>
                            </Box>
                            <Box>
                                <Typography><h2>Attention</h2></Typography>
                                <Typography><h3>What am I seeing?</h3>
                                    The model iterates over all parts of the molecule. For each of this parts, the
                                    system is distributing its attention over all parts of the molecule. E.g. if an
                                    opening parenthesis is encountered, the system may try to identify the closing
                                    counterpart. The parts of the molecule that the system is paying attention to, are
                                    indicated by lines. Darker shades indicate stronger attention.</Typography>
                                <Box sx={{height: "400px"}}>
                                    <LayerTabs layers={data.graphs}/>
                                </Box>
                            </Box>
                        </Box>
                    </Box>
                </Paper>
            </Box>
        </Box>
    )
}


