import React from 'react';

import {useEffect, useRef} from "react";

import axios from "axios";
import Grid from "@mui/material/Grid";
import Box from "@mui/material/Box";
import FormControl from "@mui/material/FormControl";
import TextField from "@mui/material/TextField";
import Button from "@mui/material/Button";
import Divider from '@mui/material/Divider';

import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemButton from '@mui/material/ListItemButton';
import ListItemIcon from '@mui/material/ListItemIcon';
import ListItemText from '@mui/material/ListItemText';

import Paper from '@mui/material/Paper';

import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';

import { DataGrid } from '@mui/x-data-grid';

const ParentsTable = (data) => {
    return (<TableContainer component={Paper} sx={{display: (data.parents.length === 0 ? 'none' : 'inline')}}>
        <Table sx={{minWidth: 650}} aria-label="simple table">
            <TableHead>
                <TableRow>
                    <TableCell>Smiles</TableCell>
                    <TableCell>Predicted parents</TableCell>
                </TableRow>
            </TableHead>
            <TableBody>
                {data.parents.map((parents, i) => <TableRow>
                    <TableCell>{data.smiles[i]}</TableCell><TableCell>{parents.join(", ")}</TableCell></TableRow>)}
            </TableBody>
        </Table>
    </TableContainer>)
}


export class ClassificationForm extends React.Component {
    constructor(props) {
        super(props);
        this.state = {text_in_single_field: '', pending_smiles: [], direct_parents: [], done_smiles: [], predicted_parents: []};
        this.handleChange = this.handleChange.bind(this);
        this.add_smiles = this.add_smiles.bind(this);
        this.process_upload = this.process_upload.bind(this);
        this.submit_smiles_list = this.submit_smiles_list.bind(this);
    }

    handleChange(event) {
        this.setState({text_in_single_field: event.target.text_in_single_field});
    }

    add_smiles(smiles_list) {
        var sm = this.state.pending_smiles.concat(smiles_list);
        this.setState({
            pending_smiles: sm,
        });
        console.log(this.state.pending_smiles)
    }

    process_upload(event) {
        event.preventDefault()
        const reader = new FileReader()
        reader.onload = async (e) => {
          this.add_smiles(e.target.result.split("\n"))
        };
        reader.readAsText(event.target.files[0])
    }

    add_single_smiles(event) {
        event.preventDefault();
        this.add_smiles([this.state.text_in_single_field]);
        this.setState({text_in_single_field: ""});
    }

    submit_smiles_list(event) {
        event.preventDefault();
        axios({
            url: '/api/classify',
            method: 'post',
            data: {smiles: this.state.pending_smiles}
        }).then(response => {
            this.setState({
                direct_parents: this.state.direct_parents.concat(response.data.direct_parents),
                pending_smiles: [],
                done_smiles: this.state.done_smiles.concat(this.state.pending_smiles),
                predicted_parents: this.state.predicted_parents.concat(response.data.predicted_parents)
            });

        });
    }

    render() {

        return (
            <Box>
                <Paper>
                    <Box>
                        <h1>ChEBIfier</h1>
                        Use artifical intelligence to classify a molecular structure in ChEBI.
                    </Box>
                    <Box>
                        <Grid container>
                            <Grid item xs>
                                <Grid container
                                      component="form"
                                      noValidate
                                      onSubmit={this.add_smiles}
                                      sx={{
                                          '& > :not(style)': {m: 1, width: '25ch'},
                                      }}
                                      spacing={2}>
                                    <Grid item xs={10} zeroMinWidth>
                                        <FormControl fullWidth variant="standard">
                                            <TextField defaultValue={this.state.text_in_single_field} onChange={this.handleChange}
                                                       label="Your SMILES string"
                                                       variant="outlined"/></FormControl></Grid>
                                    <Grid item xs={2} zeroMinWidth><FormControl><Button variant="outlined"
                                                                                        type="submit">Add
                                        smiles</Button></FormControl></Grid>
                                </Grid>
                            </Grid>
                            <Divider orientation="vertical" flexItem>
                                OR
                            </Divider>
                            <Grid item xs>
                                <input
                                    accept="text/plain"
                                    style={{display: 'none'}}
                                    id="raised-button-file"
                                    multiple
                                    type="file"
                                    onChange={this.process_upload}
                                />
                                <label htmlFor="raised-button-file">
                                    <Button variant="raised" component="span">
                                        Upload
                                    </Button>
                                </label>
                            </Grid>
                        </Grid>

                    </Box>
                    <Box sx={{ height: 400, width: '100%' }}>
                      <DataGrid
                        rows={this.state.pending_smiles.map((s) => ({"smiles": s}))}
                        columns={["Pending Smiles"]}
                        pageSize={5}
                        rowsPerPageOptions={[5]}
                        disableSelectionOnClick
                      />
                    </Box>
                    <List>{this.state.pending_smiles.map((s) => <ListItem disablePadding>
                        <ListItemButton><ListItemText secondary={s}/></ListItemButton></ListItem>)}

                    </List>
                    <Box component="form" noValidate onSubmit={this.submit_smiles_list}>
                        <Button variant="outlined" type="submit">Classify</Button>
                    </Box>
                    <ParentsTable parents={this.state.direct_parents} smiles={this.state.done_smiles}/>
                </Paper>
            </Box>)
    }

}