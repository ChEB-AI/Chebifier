import './App.css';
import ClassificationGrid from './smiles-form/classification-form.js'
import '../node_modules/vis-network/styles/vis-network.css';
import Link from '@mui/material/Link';

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <ClassificationGrid />
        <footer><Link variant="body2" href="https://github.com/ChEB-AI/Chebifier/issues">Report an issue</Link></footer>
      </header>
    </div>
  );
}

export default App;
