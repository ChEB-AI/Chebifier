from flask import Flask, send_from_directory
from flask_restful import Api, Resource, reqparse
from flask_cors import CORS  # comment this on deployment

import json

app = Flask(__name__, static_url_path='', static_folder='../react-app/build')
CORS(app) # comment this on deployment

app.config.from_file("config.default.json", load=json.load)
app.config.from_file("config.json", load=json.load)

api = Api(app)

@app.route("/", defaults={'path': ''})
def serve(path):
    return send_from_directory(app.static_folder, 'index.html')

from api.chemclass import PredictionApiHandler
api.add_resource(PredictionApiHandler, '/api/classify')
#api.add_resource(ChEBI100ApiHandler, '/api/chebi100')