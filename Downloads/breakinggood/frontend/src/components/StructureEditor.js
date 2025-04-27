import React, { useState, useEffect, useRef } from 'react';
import { makeStyles } from '@material-ui/core/styles';
import { 
  Paper, 
  Button, 
  CircularProgress, 
  Typography,
  TextField,
  FormControl,
  FormLabel,
  RadioGroup,
  FormControlLabel,
  Radio,
  Snackbar,
  IconButton,
  Grid,
  Tooltip,
  Box
} from '@material-ui/core';
import { 
  Save, 
  Refresh, 
  FileCopy, 
  Undo, 
  Redo, 
  ZoomIn, 
  ZoomOut,
  Close as CloseIcon
} from '@material-ui/icons';
import MoleculeViewer3D from './MoleculeViewer3D'; // Might need this again for 3D toggle
import axios from 'axios';

const useStyles = makeStyles((theme) => ({
  root: {
    padding: theme.spacing(3),
  },
  paper: {
    padding: theme.spacing(2),
    marginBottom: theme.spacing(3),
  },
  editorContainer: {
    border: `1px solid ${theme.palette.divider}`,
    borderRadius: theme.shape.borderRadius,
    minHeight: '400px',
    position: 'relative',
    backgroundColor: '#f9f9f9',
  },
  ketcher: { // Style for the iframe
    width: '100%',
    height: '100%',
    minHeight: '400px',
    border: 'none', // Ensure no iframe border
  },
  controlsBar: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2),
    display: 'flex',
    flexWrap: 'wrap',
    gap: theme.spacing(1),
  },
  button: {
    margin: theme.spacing(0.5),
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120,
  },
  inputField: {
    marginBottom: theme.spacing(2),
  },
  infoText: {
    margin: theme.spacing(1, 0),
  },
  propertyPaper: {
    padding: theme.spacing(2),
    height: '100%',
  },
  loading: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)',
  },
  errorMessage: {
    color: theme.palette.error.main,
    marginTop: theme.spacing(1),
  },
  viewerContainer: {
    marginTop: theme.spacing(3),
    minHeight: '400px', // Ensure container has height
  },
  propertyGrid: {
    marginTop: theme.spacing(2),
  },
  propertyLabel: {
    fontWeight: 'bold',
  },
  propertyValue: {
    marginLeft: theme.spacing(1),
  },
}));

const StructureEditor = ({ 
  initialMolecule = '',
  onMoleculeChange = () => {},
  readOnly = false,
  showControls = true,
  showProperties = true 
}) => {
  const classes = useStyles();
  const ketcherFrame = useRef(null); // Ref for the iframe
  const ketcherInitialized = useRef(false);
  const [isLoading, setIsLoading] = useState(true);
  const [ketcher, setKetcher] = useState(null); // Holds the ketcher instance from iframe
  const [smiles, setSmiles] = useState('');
  const [molfile, setMolfile] = useState('');
  const [inputType, setInputType] = useState('smiles'); // Default to SMILES
  const [inputValue, setInputValue] = useState('');
  const [errorMessage, setErrorMessage] = useState('');
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [snackbarMessage, setSnackbarMessage] = useState('');
  const [properties, setProperties] = useState(null);
  const [is3DViewActive, setIs3DViewActive] = useState(false); // State for 3D toggle

  // Initialize Ketcher via iframe
  useEffect(() => {
    console.log('[StructureEditor Iframe] Initializing...');
    console.log('[StructureEditor Iframe] process.env.PUBLIC_URL:', process.env.PUBLIC_URL);
    let intervalId = null;
    let timeoutId = null;
    let ketcherInstance = null;
    ketcherInitialized.current = false;
    setIsLoading(true);
    setErrorMessage(''); 

    const iframe = document.createElement('iframe');
    const ketcherSrc = '/ketcher/standalone/index.html';
    console.log(`[StructureEditor Iframe] Setting iframe src to: ${ketcherSrc}`);
    iframe.setAttribute('src', ketcherSrc);
    iframe.setAttribute('class', classes.ketcher);
    iframe.setAttribute('id', 'ketcher-frame');
    iframe.style.border = 'none';

    if (ketcherFrame.current) {
      ketcherFrame.current.innerHTML = '';
    }
    ketcherFrame.current?.appendChild(iframe);

    timeoutId = setTimeout(() => {
        if (!ketcherInitialized.current) {
            console.error('[StructureEditor Iframe] Ketcher loading timed out after 10 seconds.');
            setErrorMessage('Editor failed to load. Check Ketcher path and console.');
            setIsLoading(false);
            if (intervalId) clearInterval(intervalId);
        }
    }, 10000);

    iframe.onload = () => {
      console.log('[StructureEditor Iframe] iframe loaded.');
      const ketcherWindow = iframe.contentWindow;

      intervalId = setInterval(() => {
        if (ketcherWindow && ketcherWindow.ketcher) {
          console.log('[StructureEditor Iframe] Ketcher instance found.');
          clearInterval(intervalId);
          clearTimeout(timeoutId);
          
          ketcherInstance = ketcherWindow.ketcher;
          setKetcher(ketcherInstance); // Set the instance to state
          ketcherInitialized.current = true;
          setIsLoading(false);
          
          // --- Inspect Ketcher Instance --- 
          console.log('[StructureEditor Iframe] Ketcher Instance Structure:', ketcherInstance);
          // ---------------------------------
          
          // --- Temporarily Comment Out Listener --- 
          /*
          if (ketcherInstance && ketcherInstance.editor) {
            ketcherInstance.editor.on('change', () => {
              console.log('[StructureEditor Iframe] Ketcher content changed.');
              updateMoleculeData(ketcherInstance);
            });
          } else { 
              console.warn('[SE Iframe] Editor object not found for listener.'); 
          }
          */
          // --- End Comment Out ---

          // Load initial molecule
          if (initialMolecule) { 
              console.log(`[SE Iframe] Loading initial molecule: ${initialMolecule.substring(0, 30)}...`);
              setTimeout(() => loadMolecule(initialMolecule, ketcherInstance), 1500);
          } else {
              setTimeout(() => updateMoleculeData(ketcherInstance), 1500);
          }
        } else {
          console.log('[StructureEditor Iframe] Waiting for Ketcher instance...');
        }
      }, 500);
    };

    iframe.onerror = (err) => {
      console.error('[StructureEditor Iframe] iframe loading error:', err);
      setErrorMessage('Failed to load the structure editor component (iframe error).');
      setIsLoading(false);
      if (intervalId) clearInterval(intervalId);
      if (timeoutId) clearTimeout(timeoutId);
    };

    return () => {
      console.log('[StructureEditor Iframe] Cleanup.');
      if (intervalId) clearInterval(intervalId);
      if (timeoutId) clearTimeout(timeoutId);
      // Detach listener if needed, though instance might be gone
      ketcherInitialized.current = false;
    };
  }, [initialMolecule, classes.ketcher]); // Dependencies

  // Load molecule (remains largely the same, uses ketcher state)
  const loadMolecule = async (data, ketcherInst = ketcher) => {
     if (!ketcherInst) {
         console.error('[SE Iframe] loadMolecule called before ketcher ready.');
         setErrorMessage('Editor not ready.');
         return;
     }
     try {
       setIsLoading(true);
       setErrorMessage('');
       console.log(`[SE Iframe] Loading: ${data.substring(0, 30)}...`);
       if (data.startsWith('InChI=')) {
         await ketcherInst.setMolecule(data);
       } else if (data.includes('\n') && data.includes('M  END')) {
         await ketcherInst.setMolecule(data);
       } else {
         console.log('[SE Iframe] Converting SMILES to Mol...');
         try {
            const response = await axios.post('/api/simulation/convert', {
                input: data, inputFormat: 'smiles', outputFormat: 'mol'
            });
            if (response.data && response.data.output) {
                await ketcherInst.setMolecule(response.data.output);
            } else { throw new Error('Conversion empty.'); }
         } catch(convErr) {
            console.error('[SE Iframe] Conversion failed:', convErr);
            throw new Error(`SMILES conversion error: ${convErr.response?.data?.error || convErr.message}`);
         }
       }
       await updateMoleculeData(ketcherInst);
     } catch (error) {
       console.error('[SE Iframe] Error loading molecule:', error);
       setErrorMessage(`Load failed: ${error.message}`);
     } finally {
       setIsLoading(false);
     }
  };

  // Update molecule data (remains largely the same, uses ketcher state)
  const updateMoleculeData = async (ketcherInst = ketcher) => {
     if (!ketcherInst) return;
     try {
         // --- Log Molfile Before Sending --- 
         const mol = await ketcherInst.getMolfile();
         console.log('[SE Iframe] Molfile from Ketcher:', mol); 
         // ----------------------------------
         setMolfile(mol);

         // --- Add Check for Empty Molfile --- 
         if (!mol || mol.trim() === '' || !mol.includes('M  END')) {
             console.log('[SE Iframe] Empty or invalid Molfile retrieved from Ketcher. Skipping conversion.');
             setSmiles(''); setProperties(null); onMoleculeChange('');
             return; // Don't attempt conversion if Molfile is empty/invalid
         }
         // ---------------------------------

         const response = await axios.post('/api/simulation/convert', {
             input: mol, inputFormat: 'mol', outputFormat: 'smiles'
         });
         if (response.data && response.data.output) {
             const newSmiles = response.data.output.trim();
             setSmiles(newSmiles);
             onMoleculeChange(newSmiles);
             if (showProperties && newSmiles) calculateProperties(newSmiles);
             else setProperties(null);
         } else {
             setSmiles(''); setProperties(null); onMoleculeChange('');
         }
     } catch (error) {
         console.error('[SE Iframe] Update failed:', error);
         setErrorMessage(`Update failed: ${error.message}`);
         setSmiles(''); setProperties(null); onMoleculeChange('');
     }
  };

  // Calculate properties (remains the same)
  const calculateProperties = async (smilesString) => {
    setProperties(null);
    try {
      console.log(`[SE Iframe] Calculating properties for ${smilesString}`);
      const response = await axios.post('/api/simulation/properties', {
        smiles: smilesString
      });
      if (response.data && !response.data.error) {
        setProperties(response.data);
      } else {
        setProperties({ error: response.data?.error || 'Calculation failed' });
      }
    } catch (error) {
      console.error('[SE Iframe] Properties API call failed:', error);
      setProperties({ error: 'API call failed' });
    }
  };

  // Handlers (Input type, value, load, controls, snackbar) remain largely the same
  const handleInputTypeChange = (event) => setInputType(event.target.value);
  const handleInputValueChange = (event) => setInputValue(event.target.value);
  const handleLoadMolecule = () => {
      if (!inputValue.trim()) { /* ... */ return; }
      loadMolecule(inputValue.trim());
  };
  const handleUndo = () => ketcher?.undo();
  const handleRedo = () => ketcher?.redo();
  const handleClear = () => ketcher?.clear(); // Change listener handles update
  const copySmilesToClipboard = () => { /* ... */ };
  const handleSave = async () => { /* ... */ };
  const toggle3DView = () => setIs3DViewActive(!is3DViewActive);
  const handleCloseSnackbar = (/*...*/) => setSnackbarOpen(false);

  return (
    <div className={classes.root}>
      <Paper className={classes.paper}>
        <Typography variant="h6" gutterBottom>
          Molecular Structure Editor
        </Typography>
        
        {showControls && (
          <div className={classes.controlsBar}>
             {/* Add disabled={!ketcher || isLoading} to buttons */}
            <Tooltip title="Undo"><IconButton className={classes.button} onClick={handleUndo} disabled={!ketcher || isLoading}><Undo /></IconButton></Tooltip>
            <Tooltip title="Redo"><IconButton className={classes.button} onClick={handleRedo} disabled={!ketcher || isLoading}><Redo /></IconButton></Tooltip>
            <Tooltip title="Clear"><Button variant="outlined" className={classes.button} onClick={handleClear} disabled={!ketcher || isLoading}>Clear</Button></Tooltip>
            <Tooltip title="Copy SMILES"><IconButton className={classes.button} onClick={copySmilesToClipboard} disabled={!smiles || !ketcher || isLoading}><FileCopy /></IconButton></Tooltip>
            <Tooltip title="Save Molecule"><IconButton className={classes.button} onClick={handleSave} disabled={!smiles || !ketcher || isLoading} color="primary"><Save /></IconButton></Tooltip>
            {/* 3D View Toggle Button - Rework if MoleculeViewer3D is kept */}
             <Button
               variant="contained"
               color={is3DViewActive ? "primary" : "default"}
               className={classes.button}
               onClick={toggle3DView}
               disabled={!smiles || isLoading}
             >
               {is3DViewActive ? "2D Editor" : "3D View"}
             </Button>
          </div>
        )}
        
        <div className={classes.editorContainer}>
          {/* Conditional rendering for 2D Editor (iframe) or 3D Viewer */}
          {!is3DViewActive ? (
            <div ref={ketcherFrame} className={classes.ketcher} /> // Iframe container
          ) : (
            <div className={classes.viewerContainer}>
              {/* Ensure MoleculeViewer3D is imported and working */}
              {smiles ? (
                  <MoleculeViewer3D 
                      moleculeData={smiles} 
                      format="smiles" 
                      height={400} 
                  /> 
              ) : (
                  <Typography>No molecule loaded for 3D view.</Typography>
              )}
            </div>
          )}
          
          {isLoading && (
            <div className={classes.loading}><CircularProgress /></div>
          )}
        </div>
        
        {errorMessage && (
          <Typography className={classes.errorMessage}>{errorMessage}</Typography>
        )}
        
        {/* Import Structure Section */}
        <Box mt={3}>
          {/* ... (Import controls: RadioGroup, TextField, Button) ... */}
          <Typography variant="subtitle1" gutterBottom>Import Structure</Typography>
          <FormControl component="fieldset" className={classes.formControl}>
            <FormLabel component="legend">Input Type</FormLabel>
            <RadioGroup row value={inputType} onChange={handleInputTypeChange}>
              <FormControlLabel value="smiles" control={<Radio color="primary" />} label="SMILES" />
              <FormControlLabel value="molfile" control={<Radio color="primary" />} label="Molfile" />
              <FormControlLabel value="inchi" control={<Radio color="primary" />} label="InChI" />
            </RadioGroup>
          </FormControl>
          <TextField
            label={`Enter ${inputType.toUpperCase()}`}
            variant="outlined"
            fullWidth
            multiline={inputType === 'molfile'}
            rows={inputType === 'molfile' ? 4 : 1}
            value={inputValue}
            onChange={handleInputValueChange}
            className={classes.inputField}
          />
          <Button
            variant="contained"
            color="primary"
            onClick={handleLoadMolecule}
            disabled={!inputValue.trim() || isLoading}
          >
            Load Molecule
          </Button>
        </Box>
      </Paper>
      
      {/* Properties Section */}
      {showProperties && properties && (
        <Paper className={classes.paper}>
          <Typography variant="h6" gutterBottom>Molecular Properties</Typography>
          {properties.error ? (
             <Typography color="error">Error: {properties.error}</Typography>
          ) : (
            <Grid container spacing={2} className={classes.propertyGrid}>
              {Object.entries(properties)
                .filter(([key]) => key !== 'error') 
                .map(([key, value]) => (
                  <Grid item xs={12} sm={6} md={4} key={key}>
                    <Typography variant="body1">
                      <span className={classes.propertyLabel}>{key}:</span>
                      <span className={classes.propertyValue}>
                          {typeof value === 'number' ? value.toFixed(2) : String(value)}
                      </span>
                    </Typography>
                  </Grid>
              ))}
            </Grid>
           )}
        </Paper>
      )}
      
      {/* Snackbar */}
      <Snackbar
        anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
        open={snackbarOpen}
        autoHideDuration={3000}
        onClose={handleCloseSnackbar}
        message={snackbarMessage}
        action={
          <IconButton size="small" color="inherit" onClick={handleCloseSnackbar}>
            <CloseIcon fontSize="small" />
          </IconButton>
        }
      />
    </div>
  );
};

export default StructureEditor; 