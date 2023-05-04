import {useEffect, useRef} from "react";

import {Network} from "vis-network";
import Box from '@mui/material/Box';

export function VisNetwork(data) {

    const visJsRef = useRef(null);
    useEffect(() => {
        const network =
            visJsRef.current &&
            new Network(visJsRef.current, data.graph, {
                physics: {enabled: data.physics},
                layout: data.layout,
                width: data.width || "100%",
                height: data.height || "100%",
                clickToUse: true
            });
        network.fit();
    }, [visJsRef, data]);

    return <Box>
        <div ref={visJsRef}/>
    </Box>;
};

export function plot_ontology(graph) {

   return <VisNetwork graph={graph} physics={false}
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

}