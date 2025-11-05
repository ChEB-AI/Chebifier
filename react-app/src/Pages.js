import React from 'react';
import { useLocation } from 'react-router-dom';
import Box from '@mui/material/Box';
import ClassificationGrid from './smiles-form/classification-form';
import About from './About';

// This component stays mounted across route changes (used under a wildcard route)
// and toggles visibility of pages based on the URL, so each page preserves state.
export default function Pages() {
  const location = useLocation();
  const isAbout = location.pathname.toLowerCase().includes('/about');

  return (
    <Box>
      <Box sx={{ display: isAbout ? 'none' : 'block' }}>
        <ClassificationGrid />
      </Box>
      <Box sx={{ display: isAbout ? 'block' : 'none' }}>
        <About />
      </Box>
    </Box>
  );
}
