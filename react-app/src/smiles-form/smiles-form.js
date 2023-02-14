import React from 'react';
import ReactDOM from 'react-dom';
import axios from 'axios'

export class SmilesForm extends React.Component {
  constructor(props) {
    super(props);
    this.state = {value: ''};
    this.handleChange = this.handleChange.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
  }

  handleChange(event) {    this.setState({value: event.target.value});  }
  handleSubmit(event) {
    axios.post('/api/classify', { smiles: this.state.value }).then(response => {
        console.log(response);
    });
    event.preventDefault();
  }

  render() {
    return (
      <form onSubmit={this.handleSubmit}>        <label>
          Smiles:
          <input type="text" defaultValue={this.state.value} onChange={this.handleChange} />        </label>
        <input type="submit" value="Submit" />
      </form>
    );
  }
}