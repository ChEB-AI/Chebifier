import logo from './logo.svg';
import './App.css';
import ClassificationGrid from './smiles-form/classification-form.js'

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <img src={logo} className="App-logo" alt="logo" />
        <ClassificationGrid />
      </header>
    </div>
  );
}

export default App;
