import {useEffect, useRef} from "react";
import * as React from 'react';
import {Network} from "vis-network";

import ArrowDownwardIcon from '@mui/icons-material/ArrowDownward';
import ArrowUpwardIcon from '@mui/icons-material/ArrowUpward';
import Box from '@mui/material/Box';
import Grid from '@mui/material/Grid';
import Link from '@mui/material/Link';
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemText from '@mui/material/ListItemText';
import ListSubheader from '@mui/material/ListSubheader';
import Typography from '@mui/material/Typography';


function buildNode(id, node, includeLabel=true){
    const d = {id:id}
    d["title"] = node["lbl"]
    if(!node.artificial){
        d["lbl"] = d["title"]
    } else {
        d["color"] = "#EEEEEC"
    }
    return d
}

function buildEdge(id, edge){
    return {from: edge[0], to:edge[1], arrows:{to:true}}
}

function renderClassListElement(s, node){
    const node_id = String(s)
    if(node_id.startsWith("CHEBI:")){
        return (<ListItemText><Link href={node_id.replace("CHEBI:", "http://purl.obolibrary.org/obo/CHEBI_")}>{node["title"] || node_id}</Link></ListItemText>)
    } else {
        return <ListItemText>{node["title"] || s}</ListItemText>
    }
}

function subheader(list){
    return list.map((e) => (
                (<ListItem>
                    {e}
                </ListItem>)
            ))
    }

function renderOverview(node, graph){

    if(graph == null || node == null){

        return <Typography> Select a class </Typography>

    }
    const nodeDict = Object.fromEntries(graph.nodes.map(x => [x["id"], x]));

    console.log(graph);
    const superclasses = graph.edges.filter((e) => (e["from"] == node)).map((e) => renderClassListElement(e["to"], nodeDict[e["to"]]))
    const subclasses = graph.edges.filter((e) => (e["to"] == node)).map((e) => renderClassListElement(e["from"], nodeDict[e["from"]]))

    return (<Box>
        <List dense={true}>
            <ListSubheader>
                This class
            </ListSubheader>
            <ListItem>
                <ListItemText>{nodeDict[node]["title"] || node}</ListItemText>
            </ListItem>
            <ListSubheader>
                <ArrowUpwardIcon/>Superclasses
            </ListSubheader>
            {subheader(superclasses)}
            <ListSubheader>
                <ArrowDownwardIcon/>Subclasses
            </ListSubheader>
            {subheader(subclasses)}
        </List>
    </Box>)
}

export function VisNetwork(data) {

    const visJsRef = useRef(null);
    const [selectedNode, setSelectedNode] = React.useState(null);
    const [graph, setGraph] = React.useState(null);
    useEffect(() => {
        if(data.graph != null){
            console.log(data.graph);
            const g = {
                nodes: Object.keys(data.graph.nodes).map(k => buildNode(k, data.graph.nodes[k])),
                edges: Object.keys(data.graph.edges).map(k => buildEdge(k, data.graph.edges[k]))
            }
            const network =
                visJsRef.current &&
                new Network(visJsRef.current, g, {
                    physics: {enabled: data.physics},
                    layout: data.layout,
                    interaction: data.interaction,
                    width: data.width || "100%",
                    height: data.height || "100%",
                    clickToUse: true
                });
            network.fit();
            network.on("selectNode", function (params) {
                console.log(params);
                setSelectedNode(params.nodes[0] || null);
            });
            setGraph(g)
        }

    }, [visJsRef, data]);
    return <Grid container spacing={2}>
          <Grid item xs={9}>
            <Box><div ref={visJsRef}/></Box>
          </Grid>
          <Grid item xs={3}>
            <Box>{renderOverview(selectedNode, graph)}</Box>
          </Grid>
        </Grid>
};

export function plot_ontology(graph, physics=false, hierarchical=true) {
    var layout = null;
    if(hierarchical){
        layout={
            hierarchical: {
                enabled: true,
                direction: "LR",
                sortMethod: "directed",
                levelSeparation: 250,
            }
        }
    } else {
        layout={
            hierarchical: {
                enabled: true,
                direction: "LR",
                sortMethod: "directed",
                levelSeparation: 250,
            }
        }
        physics = true;
    }

    var interaction = {
        dragNodes:true,
        dragView: true,
        hideEdgesOnDrag: false,
        hideEdgesOnZoom: false,
        hideNodesOnDrag: false,
        hover: true,
        hoverConnectedEdges: true,
        keyboard: {
          enabled: false,
          speed: {x: 10, y: 10, zoom: 0.02},
          bindToWindow: true,
          autoFocus: true,
        },
        multiselect: false,
        navigationButtons: false,
        selectable: true,
        selectConnectedEdges: true,
        tooltipDelay: 300,
        zoomSpeed: 1,
        zoomView: true
    }

    return <VisNetwork graph={graph} physics={physics}
                                        layout={layout}
                                        interaction={interaction}
                                        height={"400px"}
                            />

}