from flask import Flask, send_from_directory
from flask_restful import Api, Resource, reqparse
from flask_cors import CORS  # comment this on deployment
import torch.multiprocessing as mp

import json

app = Flask(__name__, static_url_path='', static_folder='../react-app/build')
CORS(app) # comment this on deployment

app.config.from_file("config.template.json", load=json.load)
app.config.from_file("config.json", load=json.load)

api = Api(app)


@app.route("/", defaults={'path': ''})
def serve(path):
    return send_from_directory(app.static_folder, 'index.html')


def load_endpoints():
    from api.chemclass import PredictionDetailApiHandler, BatchPrediction, \
        HierarchyAPI
    api.add_resource(PredictionDetailApiHandler, '/api/details')
    api.add_resource(BatchPrediction, '/api/classify')
    api.add_resource(HierarchyAPI, '/api/hierarchy')


with app.app_context():
    mp.set_start_method("spawn")
    load_endpoints()
