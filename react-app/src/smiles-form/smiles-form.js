import React from 'react';
import axios from 'axios'
import Box from '@mui/material/Box';
import Grid from '@mui/material/Grid';
import { useEffect, useRef } from "react";
import { Network } from "vis-network";
import Tab from '@mui/material/Tab';
import TabContext from '@mui/lab/TabContext';
import TabList from '@mui/lab/TabList';
import TabPanel from '@mui/lab/TabPanel';

const VisNetwork = (graph) => {

	const visJsRef = useRef(null);
	useEffect(() => {
		const network =
			visJsRef.current &&
			new Network(visJsRef.current, graph.graph, {physics: {enabled:false}} );
	}, [visJsRef, graph]);

	return <Box sx={{ width: '100%', typography: 'body1' }}><div ref={visJsRef} /></Box>;
};

const LayerComponent = (layer) => {
    console.log(layer);
    return <Grid> {layer.layer.map((g) => <VisNetwork graph={g} />)} </Grid>
}


export default function LayerTabs(layers) {
  const [value, setValue] = React.useState('1');

  const handleChange = (event, newValue) => {
    setValue(newValue);
  };

  return (
    <Box sx={{ width: '100%', typography: 'body1' }}>
      <TabContext value={value}>
        <Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
          <TabList onChange={handleChange} aria-label="lab API tabs example">
            {layers.layers.map((layer,i) => <Tab label={"Layer " + i} value={i} />)}
          </TabList>
        </Box>
          {layers.layers.map((layer,i) => <TabPanel value={i}><LayerComponent layer={layer} /></TabPanel>)}
      </TabContext>
    </Box>
  );
}


export class SmilesForm extends React.Component {
  constructor(props) {
    super(props);
    this.state = {value: '', attention_fig: '', graphs:[]};
    this.handleChange = this.handleChange.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
  }

  handleChange(event) {    this.setState({value: event.target.value});  }
  handleSubmit(event) {
    event.preventDefault();
    axios.post('/api/classify', { smiles: this.state.value }).then(response => {
        this.setState({attention_fig:response.data.figures.attention_mol, graphs:response.data.graphs});
    });
  }



  render() {
    const Example = ({ data }) => <img src={`data:image/jpeg;base64,${data}`} />

    return (
      <div>
          <form onSubmit={this.handleSubmit}>        <label>
              Smiles:
              <input type="text" defaultValue={this.state.value} onChange={this.handleChange} />        </label>
            <input type="submit" value="Submit" />
          </form>
          <Example data={this.state.attention_fig} />
          <LayerTabs layers={this.state.graphs}/>
      </div>
    );
  }
}