import React from 'react';
import axios from 'axios'

import {useEffect, useRef} from "react";
import {Network} from "vis-network";

import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import Box from '@mui/material/Box';
import Button from '@mui/material/Button';
import Divider from '@mui/material/Divider';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
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

const VisNetwork = (data) => {

    const visJsRef = useRef(null);
    useEffect(() => {
        const network =
            visJsRef.current &&
            new Network(visJsRef.current, data.graph, {physics: {enabled: data.physics}, layout: data.layout, width:data.width || "100%", height:data.height || "100%", clickToUse: true});
        network.fit();
    }, [visJsRef, data]);

    return <Box>
        <div ref={visJsRef}/>
    </Box>;
};

const LayerComponent = (layer) => {
    const Item = styled(Paper)(({theme}) => ({
        backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
        ...theme.typography.body2,
        padding: theme.spacing(1),
        textAlign: 'center',
        color: theme.palette.text.secondary,
    }));

    return <Stack> {layer.layer.map((g) => <Item><VisNetwork graph={g} physics={false} layout={{}}/></Item>)} </Stack>
}


export default function LayerTabs(layers) {
    const [value, setValue] = React.useState('1');

    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    return (
        <Box sx={{width: '100%', typography: 'body1'}}>
            <TabContext value={value}>
                <Box sx={{borderBottom: 1, borderColor: 'divider'}}>
                    <TabList onChange={handleChange} aria-label="lab API tabs example">
                        {layers.layers.map((layer, i) => <Tab label={"Layer " + i} value={i}/>)}
                    </TabList>
                </Box>
                {layers.layers.map((layer, i) => <TabPanel value={i}><LayerComponent layer={layer}/></TabPanel>)}
            </TabContext>
        </Box>
    );
}

const Item = styled(Paper)(({theme}) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.text.secondary,
}));

const Molecule = ({data}) => <img src={`data:image/jpeg;base64,${data}`} hidden={(data === null)} width="200" height="200"/>;

export class SmilesForm extends React.Component {
    constructor(props) {
        super(props);
        this.state = {value: '', attention_fig: '', graphs: [], chebi: null};
        this.handleChange = this.handleChange.bind(this);
        this.handleSubmit = this.handleSubmit.bind(this);
    }

    handleChange(event) {
        this.setState({value: event.target.value});
    }

    handleSubmit(event) {
        event.preventDefault();
        axios.post('/api/classify', {smiles: this.state.value}).then(response => {
            this.setState({
                attention_fig: response.data.figures.attention_mol,
                graphs: response.data.graphs,
                chebi: response.data.classification
            });
        });
    }


    render() {


        return (
            <Box width={'100%'}>
                <Grid container spacing={2}>
                    <Grid item xs={8}>
                        <Item>
                            <Box>
                                <Box>
                                    <h1>ChEBIfier</h1>
                                    Use artifical intelligence to classify a molecular structure in ChEBI.
                                </Box>
                                <Box component="form"
                                     sx={{
                                         '& > :not(style)': {m: 1, width: '25ch'},
                                     }}
                                     noValidate
                                     onSubmit={this.handleSubmit}>
                                    <TextField defaultValue={this.state.value} onChange={this.handleChange}
                                               label="Your SMILES string"
                                               variant="outlined"/>
                                    <Button variant="outlined" type="submit">Classify</Button>
                                </Box>
                            </Box>
                        </Item>
                    </Grid>
                    <Grid item xs={4}>
                        <Box sx={{ display: (this.state.chebi === null ?'none':'inline') }}>
                            <Item><Molecule data={this.state.attention_fig}/></Item>
                        </Box>
                    </Grid>
                </Grid>
                <Box sx={{ display: (this.state.chebi === null ?'none':'inline') }}>
                    <Accordion>
                        <AccordionSummary
                            expandIcon={<ExpandMoreIcon/>}
                            aria-controls="panel1a-content"
                            id="panel1a-header"
                        >
                            <Typography><h2>ChEBI Classification</h2></Typography>
                        </AccordionSummary>
                        <AccordionDetails>
                            <Paper>
                                    <Typography><h3>What am I seeing?</h3>
                                    The graph below shows you the position in the ChEBI ontology that our system
                                    proposed. Not that this prediction is an estimate based on available data and may be
                                    prone to errors.</Typography>
                                </Paper>
                            <Box>
                                <VisNetwork graph={this.state.chebi} physics={false}
                                            layout={{
                                                hierarchical: {
                                                    enabled: true,
                                                    direction: "LR",
                                                    sortMethod: "directed",
                                                    levelSeparation: 250,
                                                }
                                            }}
                                            height={"400px"}
                                />
                            </Box>
                        </AccordionDetails>
                    </Accordion>
                    <Accordion>
                        <AccordionSummary
                            expandIcon={<ExpandMoreIcon/>}
                            aria-controls="panel1a-content"
                            id="panel1b-header"
                        >
                            <Typography><h2>Attention</h2></Typography>
                        </AccordionSummary>
                        <AccordionDetails>
                                <Paper>
                                    <Typography><h3>What am I seeing?</h3>
                                    The model iterates over all parts of the molecule. For each of this parts, the
                                    system is distributing its attention over all parts of the molecule. E.g. if an
                                    opening parenthesis is encountered, the system may try to identify the closing
                                    counterpart. The parts of the molecule that the system is paying attention to, are
                                    indicated by lines. Darker shades indicate stronger attention.</Typography>
                                </Paper>
                                <LayerTabs layers={this.state.graphs}/>
                        </AccordionDetails>
                    </Accordion>
                </Box>
            </Box>
        );
    }
}