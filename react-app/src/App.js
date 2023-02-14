import logo from './logo.svg';
import './App.css';
import {SmilesForm} from './smiles-form/smiles-form.js'

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <img src={logo} className="App-logo" alt="logo" />
        <SmilesForm />
      </header>
    </div>
  );
}

export default App;
