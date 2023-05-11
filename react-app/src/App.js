import logo from './logo.svg';
import './App.css';
import ClassificationGrid from './smiles-form/classification-form.js'
import '../node_modules/vis-network/styles/vis-network.css';

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <ClassificationGrid />
      </header>
    </div>
  );
}

export default App;
