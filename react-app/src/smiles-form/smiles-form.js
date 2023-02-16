import React from 'react';
import ReactDOM from 'react-dom';
import axios from 'axios'

export class SmilesForm extends React.Component {
  constructor(props) {
    super(props);
    this.state = {value: '', attention_fig: ''};
    this.handleChange = this.handleChange.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
  }

  handleChange(event) {    this.setState({value: event.target.value});  }
  handleSubmit(event) {
    event.preventDefault();
    axios.post('/api/classify', { smiles: this.state.value }).then(response => {
        this.setState({attention_fig:response.data.attention_fig});
    });
  }



  render() {
    const Example = ({ data }) => <img src={`data:image/jpeg;base64,${data}`} />
    return (
      <form onSubmit={this.handleSubmit}>        <label>
          Smiles:
          <input type="text" defaultValue={this.state.value} onChange={this.handleChange} />        </label>
        <input type="submit" value="Submit" />
      </form>
      <Example data={this.state.attention_fig} />
    );
  }
}