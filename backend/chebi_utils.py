import queue

from app import app
import json
import networkx as nx


LABEL_HIERARCHY = json.load(open(app.config["CHEBI_JSON"], encoding="utf-8"))

def load_sub_ontology():
    g = nx.DiGraph()
    g.add_nodes_from([(c["ID"], dict(lbl=c["LABEL"][0])) for c in LABEL_HIERARCHY])
    g.add_edges_from([(t, c["ID"]) for c in LABEL_HIERARCHY for t in c.get("SubClasses",[])])
    return g

CHEBI_FRAGMENT = load_sub_ontology()

def get_transitive_predictions(predicted_classes):
    # get all parents of predicted classes
    all_predicted = [cls for clss in predicted_classes for cls in clss]
    q = queue.Queue()
    for cls in all_predicted:
        q.put(cls)
    while not q.empty():
        cls = q.get()
        for parent in CHEBI_FRAGMENT.successors(cls):
            if parent not in all_predicted:
                all_predicted.append(parent)
                q.put(parent)
    return all_predicted
