import React, { useState, useEffect } from 'react';
import {
  Typography,
  Paper,
  Grid,
  Button,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  CircularProgress,
  Divider,
  Card,
  CardContent,
  makeStyles,
  Slider,
  Chip,
  IconButton,
  Radio,
  RadioGroup,
  FormControlLabel,
  FormLabel,
  Box,
  Snackbar
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import MemoryIcon from '@material-ui/icons/Memory';
import ShowChartIcon from '@material-ui/icons/ShowChart';
import SubjectIcon from '@material-ui/icons/Subject';
import AiIcon from '@material-ui/icons/Psychology';
import SearchIcon from '@material-ui/icons/Search';
import LibraryBooksIcon from '@material-ui/icons/LibraryBooks';
import axios from 'axios';
import { simulationAPI, api } from '../services/api';
import MoleculeViewer3D from './MoleculeViewer3DImproved';
import { Radar } from 'react-chartjs-2';
import {
  Chart as ChartJS, 
  RadialLinearScale, 
  PointElement, 
  LineElement, 
  Filler, 
  Tooltip, 
  Legend 
} from 'chart.js';

ChartJS.register(
  RadialLinearScale, 
  PointElement, 
  LineElement, 
  Filler, 
  Tooltip, 
  Legend
);

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    padding: theme.spacing(3),
  },
  paper: {
    padding: theme.spacing(3),
    marginBottom: theme.spacing(3),
  },
  formControl: {
    marginBottom: theme.spacing(2),
    minWidth: '100%',
  },
  button: {
    marginTop: theme.spacing(2),
  },
  simulationCard: {
    marginBottom: theme.spacing(2),
  },
  simulationHeader: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: theme.spacing(2),
  },
  tabs: {
    marginBottom: theme.spacing(2),
  },
  chip: {
    margin: theme.spacing(0.5),
  },
  simulationResults: {
    padding: theme.spacing(2),
    backgroundColor: '#f5f5f5',
    borderRadius: theme.shape.borderRadius,
    marginTop: theme.spacing(2),
  },
  chartPlaceholder: {
    height: 300,
    border: `1px solid ${theme.palette.divider}`,
    backgroundColor: '#eee',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    borderRadius: theme.shape.borderRadius,
    marginTop: theme.spacing(2),
  },
  thinkingBlock: {
    backgroundColor: '#e3f2fd',
    padding: theme.spacing(2),
    borderRadius: theme.shape.borderRadius,
    marginBottom: theme.spacing(2),
    whiteSpace: 'pre-wrap',
    maxHeight: '400px',
    overflowY: 'auto',
    fontFamily: 'monospace',
    fontSize: '0.9rem'
  },
  responseBlock: {
    backgroundColor: '#f5f5f5',
    padding: theme.spacing(2),
    borderRadius: theme.shape.borderRadius,
    marginBottom: theme.spacing(2),
    whiteSpace: 'pre-wrap'
  },
  promptField: {
    marginBottom: theme.spacing(2)
  },
  resultsContainer: {
    marginTop: theme.spacing(3),
  },
  propertiesTable: {
    width: '100%',
    marginTop: theme.spacing(1),
    borderCollapse: 'collapse',
    '& td': {
      padding: theme.spacing(0.75, 1),
      borderBottom: `1px solid ${theme.palette.divider}`,
      verticalAlign: 'top',
    },
    '& td:first-child': {
      fontWeight: 'bold',
      width: '40%',
    }
  },
  admetSection: {
    marginTop: theme.spacing(2),
  },
  viewerPlaceholder: {
    minHeight: '300px',
    border: `1px solid ${theme.palette.divider}`,
    borderRadius: theme.shape.borderRadius,
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    backgroundColor: '#f9f9f9',
    textAlign: 'center',
    padding: theme.spacing(2),
    marginBottom: theme.spacing(2),
  }
}));

function Simulations() {
  const classes = useStyles();
  const [simulationType, setSimulationType] = useState('properties');
  const [targetProtein, setTargetProtein] = useState('DAT');
  const [savedMolecules, setSavedMolecules] = useState([]);
  const [selectedMoleculeId, setSelectedMoleculeId] = useState('');
  const [manualSmiles, setManualSmiles] = useState('');
  const [inputMethod, setInputMethod] = useState('select');
  const [loading, setLoading] = useState(false);
  const [simulationResult, setSimulationResult] = useState(null);
  const [error, setError] = useState(null);
  const [snackbar, setSnackbar] = useState({ open: false, message: '', severity: 'info' });

  useEffect(() => {
    try {
      const saved = localStorage.getItem('savedMolecules');
      console.log('[Simulations] Raw data from localStorage:', saved);
      if (saved) {
        const parsedMolecules = JSON.parse(saved);
        console.log('[Simulations] Parsed saved molecules:', parsedMolecules);
        if (Array.isArray(parsedMolecules)) {
            setSavedMolecules(parsedMolecules);
        } else {
            console.error('[Simulations] Data loaded from localStorage is not an array:', parsedMolecules);
            setSavedMolecules([]);
        }
      }
    } catch (err) {
      console.error('Error loading saved molecules:', err);
      setSnackbar({ open: true, message: 'Failed to load saved molecules', severity: 'error' });
    }
  }, []);

  const handleSnackbarClose = (event, reason) => {
    if (reason === 'clickaway') return;
    setSnackbar({ ...snackbar, open: false });
  };

  const getSelectedMolecule = () => {
    if (inputMethod === 'select') {
      return savedMolecules.find(m => m.id === selectedMoleculeId) || null;
    } else {
      if (manualSmiles && manualSmiles.trim().length > 3) {
        try {
          if (!manualSmiles.includes('C') && !manualSmiles.includes('c')) throw new Error('Invalid SMILES (must contain C)');
          return { id: 'manual-'+Date.now(), name: 'Manually Entered', smiles: manualSmiles.trim() };
        } catch (e) {
          setError('Invalid SMILES string entered.');
          return null;
        }
      }
      return null;
    }
  };

  const runSimulation = async () => {
    const currentMolecule = getSelectedMolecule();
    
    if (!currentMolecule || !currentMolecule.smiles) {
      setError('Please select a valid molecule or enter a valid SMILES string.');
      setSnackbar({ open: true, message: 'Select or enter a valid molecule SMILES.', severity: 'warning' });
      return;
    }

    setLoading(true);
    setError(null);
    setSimulationResult(null);
    const smilesToSimulate = currentMolecule.smiles;
    const moleculeName = currentMolecule.name;

    console.log(`Running simulation: ${simulationType} for ${moleculeName} (${smilesToSimulate})`);

    try {
      let response = null;
      let resultPayload = null;

      switch (simulationType) {
        case 'properties':
          console.log(`POST /api/simulation/properties with SMILES: ${smilesToSimulate}`);
          response = await api.post('/simulation/properties', { smiles: smilesToSimulate });
          resultPayload = response.data?.properties;
          if (!resultPayload || typeof resultPayload !== 'object') throw new Error('Invalid properties data received from API');
          break;
        case 'admet':
          console.log(`POST /api/simulation/admet with SMILES: ${smilesToSimulate}`);
          response = await api.post('/api/simulation/admet', { smiles: smilesToSimulate });
          resultPayload = response.data?.admet;
          if (!resultPayload || typeof resultPayload !== 'object') throw new Error('Invalid ADMET data received from API');
          break;
        case 'docking':
          if (!targetProtein) throw new Error('Target protein must be selected for docking.');
          console.log(`POST /api/simulation/docking against ${targetProtein}`);
          response = await api.post('/simulation/docking', { 
            ligandSmiles: smilesToSimulate, 
            receptorPdb: targetProtein,
          });
          resultPayload = response.data?.dockingResults || response.data;
          if (!resultPayload || typeof resultPayload !== 'object') throw new Error('Invalid docking data received from API');
          break;
        default:
          throw new Error(`Unsupported simulation type: ${simulationType}`);
      }
      
      console.log('Simulation API Response Data:', resultPayload);
      setSimulationResult({ 
        molecule: currentMolecule,
        resultType: simulationType, 
        data: resultPayload
      });
      setSnackbar({ open: true, message: `${simulationType} simulation completed!`, severity: 'success' });

    } catch (err) {
      console.error(`Error running ${simulationType} simulation:`, err);
      const errorDetail = err.response?.data?.details || err.response?.data?.error || err.message || 'Unknown error';
      setError(`Simulation failed: ${errorDetail}`);
      setSnackbar({ open: true, message: `Simulation failed: ${errorDetail}`, severity: 'error' });
      setSimulationResult(null);
    } finally {
      setLoading(false);
    }
  };

  const simulationTypes = [
    { value: 'properties', label: 'Molecular Properties' },
    { value: 'admet', label: 'ADMET Prediction' },
    { value: 'docking', label: 'Molecular Docking' },
  ];
  
  const formatAdmetForRadar = (admetData) => {
    const radarLabels = ['Solubility', 'BBB Perm', 'CYP Inhib', 'Hepatotox', 'Cardiotox', 'Absorption'];
    const radarValues = [
      admetData?.solubility_logS ?? 0,
      admetData?.bbb_permeability_prob ?? 0,
      admetData?.cyp_inhibition_prob ?? 0,
      admetData?.hepatotoxicity_prob ?? 0,
      admetData?.cardiotoxicity_prob ?? 0,
      admetData?.oral_absorption_percent ?? 0,
    ];
    
    return {
      labels: radarLabels,
      datasets: [
        {
          label: 'Predicted ADMET Profile',
          data: radarValues,
          backgroundColor: 'rgba(223, 138, 101, 0.2)',
          borderColor: '#DF8A65',
          borderWidth: 1,
          pointBackgroundColor: '#DF8A65',
        },
      ],
    };
  };

  const renderResults = () => {
    if (!simulationResult) return null;
    
    const { molecule, resultType, data } = simulationResult;
    
    if (!data || typeof data !== 'object') {
      return <Typography color="error">Simulation data is missing or invalid.</Typography>;
    }

    return (
      <Box mt={2}>
        <Box mb={3}>
          <Typography variant="h6" gutterBottom>3D Structure</Typography>
          <Paper elevation={1} className={classes.viewerPlaceholder}>
            {molecule.smiles ? (
              <MoleculeViewer3D smiles={molecule.smiles} height={300}/>
            ) : (
              <Typography color="textSecondary">No SMILES available for 3D view.</Typography>
            )}
          </Paper>
        </Box>

        {resultType === 'properties' && data && (
          <Paper elevation={1} style={{ padding: '16px' }}>
            <Typography variant="h6" gutterBottom>Molecular Properties</Typography>
            <table className={classes.propertiesTable}>
              <tbody>
                {Object.entries(data).map(([key, value]) => (
                  <tr key={key}>
                    <td><Typography variant="body2">{key}:</Typography></td>
                    <td><Typography variant="body2">{typeof value === 'number' ? value.toFixed(2) : String(value)}</Typography></td>
                  </tr>
                ))}
              </tbody>
            </table>
          </Paper>
        )}

        {resultType === 'admet' && data && (
          <Paper elevation={1} style={{ padding: '16px', marginTop: '16px' }}>
            <Typography variant="h6" gutterBottom>ADMET Properties</Typography>
            <table className={classes.propertiesTable}>
              <tbody>
                {Object.entries(data).map(([key, value]) => (
                  <tr key={key}>
                    <td><Typography variant="body2">{key}:</Typography></td>
                    <td><Typography variant="body2">{typeof value === 'number' ? value.toFixed(2) : String(value)}</Typography></td>
                  </tr>
                ))}
              </tbody>
            </table>
            
            <Typography variant="subtitle1" style={{marginTop: '16px'}}>ADMET Visualization</Typography>
            <Box className={classes.chartPlaceholder} style={{ minHeight: '300px', padding: '16px' }}>
              <Radar data={formatAdmetForRadar(data)} options={{ maintainAspectRatio: false }} /> 
            </Box>
          </Paper>
        )}

        {resultType === 'docking' && data && (
          <Paper elevation={1} style={{ padding: '16px', marginTop: '16px' }}>
            <Typography variant="h6" gutterBottom>Docking Results</Typography>
            <pre style={{ maxHeight: '200px', overflowY: 'auto' }}>{JSON.stringify(data, null, 2)}</pre>
          </Paper>
        )}
      </Box>
    );
  };

  return (
    <div className={classes.root}>
      <Typography variant="h4" gutterBottom>
        Simulations & Analysis
      </Typography>
      
      <Paper className={classes.paper}>
        <Grid container spacing={3}>
          <Grid item xs={12} md={4}>
            <Typography variant="h6" gutterBottom>Simulation Parameters</Typography>

            <FormControl component="fieldset" className={classes.formControl}>
              <FormLabel component="legend">Molecule Input</FormLabel>
              <RadioGroup row value={inputMethod} onChange={(e) => { setInputMethod(e.target.value); setError(null); }}>
                <FormControlLabel value="select" control={<Radio color="primary" />} label="Select Saved" />
                <FormControlLabel value="enter" control={<Radio color="primary" />} label="Enter SMILES" />
              </RadioGroup>
            </FormControl>

            {inputMethod === 'select' ? (
              <FormControl variant="outlined" className={classes.formControl} disabled={savedMolecules.length === 0}>
                <InputLabel>Select Molecule</InputLabel>
                <Select
                  value={selectedMoleculeId}
                  onChange={(e) => { setSelectedMoleculeId(e.target.value); setError(null); }}
                  label="Select Molecule"
                >
                  <MenuItem value=""><em>Select...</em></MenuItem>
                  {console.log('[Simulations] Rendering dropdown, savedMolecules state:', savedMolecules)}
                  {Array.isArray(savedMolecules) && savedMolecules.map((molecule) => (
                    console.log('[Simulations] Mapping molecule to MenuItem:', molecule),
                    (molecule && molecule.id && molecule.name) ? (
                        <MenuItem key={molecule.id} value={molecule.id}>
                        {molecule.name} ({molecule.smiles ? molecule.smiles.substring(0, 20) : 'No SMILES'}...)
                        </MenuItem>
                    ) : null
                  ))}
                </Select>
                {savedMolecules.length === 0 && <Typography variant="caption" color="textSecondary">Load or generate molecules first.</Typography>}
              </FormControl>
            ) : (
              <TextField
                label="Enter SMILES String"
                variant="outlined"
                fullWidth
                value={manualSmiles}
                onChange={(e) => { setManualSmiles(e.target.value); setError(null); }}
                className={classes.formControl}
              />
            )}
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel>Simulation Type</InputLabel>
              <Select
                value={simulationType}
                onChange={(e) => setSimulationType(e.target.value)}
                label="Simulation Type"
              >
                {simulationTypes.map((type) => (
                  <MenuItem key={type.value} value={type.value}>
                    {type.label}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>

            {simulationType === 'docking' && (
              <TextField 
                label="Target (PDB ID / Name)"
                variant="outlined"
                fullWidth
                value={targetProtein} 
                onChange={(e) => setTargetProtein(e.target.value)}
                className={classes.formControl}
                helperText="Example: 6W7J. Backend needs PDB handling."
              />
            )}
            
            <Button
              variant="contained"
              color="primary"
              className={classes.button}
              onClick={runSimulation}
              disabled={loading || !getSelectedMolecule()}
              fullWidth
              startIcon={loading ? <CircularProgress size={20} /> : <MemoryIcon />}
            >
              {loading ? 'Running...' : 'Run Simulation'}
            </Button>
            <Button
              variant="outlined"
              className={classes.button}
              onClick={() => { setSimulationResult(null); setError(null); }}
              fullWidth
              disabled={!simulationResult && !error}
            >
              Clear Results
            </Button>

            {error && (
              <Alert severity="error" style={{ marginTop: 16 }}>{error}</Alert>
            )}
          </Grid>

          <Grid item xs={12} md={8}>
            <Typography variant="h6" gutterBottom>
              {simulationResult ? 
                `