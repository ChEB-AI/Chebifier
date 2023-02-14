from flask_restful import Api, Resource, reqparse


class PredictionApiHandler(Resource):
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('type', type=str)
        parser.add_argument('smiles', type=str)

        args = parser.parse_args()

        print(args)
        # note, the post req from frontend needs to match the strings here (e.g. 'type and 'message')

        request_type = args['type']
        smiles = args['smiles']



        return smiles
