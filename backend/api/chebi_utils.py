from app import app
import json
import networkx as nx


PREDICTION_HEADERS = ["CHEBI:" + r.strip() for r in open(app.config["CLASS_HEADERS"])]
LABEL_HIERARCHY = json.load(open(app.config["CHEBI_JSON"], encoding="utf-8"))

def load_sub_ontology():
    g = nx.DiGraph()
    g.add_nodes_from([(c["ID"], dict(lbl=c["LABEL"][0])) for c in LABEL_HIERARCHY])
    g.add_edges_from([(t, c["ID"]) for c in LABEL_HIERARCHY for t in c.get("SubClasses",[])])
    return g

CHEBI_FRAGMENT = load_sub_ontology()
