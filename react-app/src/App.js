import logo from './logo.svg';
import './App.css';
import {ClassificationForm} from './smiles-form/classification-form.js'

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <img src={logo} className="App-logo" alt="logo" />
        <ClassificationForm />
      </header>
    </div>
  );
}

export default App;
