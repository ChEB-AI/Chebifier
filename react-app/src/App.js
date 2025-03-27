import logo from './logo.svg';
import './App.css';
import Navbar from './navbar.js';
import ClassificationGrid from './smiles-form/classification-form.js'
import '../node_modules/vis-network/styles/vis-network.css';
import Box from '@mui/material/Box';
import Link from '@mui/material/Link';
import { BrowserRouter, Routes, Route } from "react-router-dom";
import About from "./About";

function App() {
  return (
  	<BrowserRouter>
      <Routes>
        <Route path="/" element={<Navbar />}>
          <Route index element={<ClassificationGrid />} />
          <Route path="about" element={<About />} />
        </Route>
      </Routes>
    </BrowserRouter>

  );
}

export default App;
