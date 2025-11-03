import queue

from app import app
import json
import networkx as nx
from chebifier.utils import load_chebi_graph, build_chebi_graph


CHEBI_FRAGMENT = build_chebi_graph(244)


def get_transitive_predictions(predicted_classes):
    # get all parents of predicted classes
    all_predicted = [cls for clss in predicted_classes for cls in clss]
    q = queue.Queue()
    for cls in all_predicted:
        q.put(cls)
    while not q.empty():
        cls = q.get()
        for parent in CHEBI_FRAGMENT.predecessors(cls):
            if parent not in all_predicted:
                all_predicted.append(parent)
                q.put(parent)
    return all_predicted
