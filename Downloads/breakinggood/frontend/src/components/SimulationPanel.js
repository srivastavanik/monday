import React, { useState, useEffect, useRef } from 'react';
import { 
  Typography, 
  Grid, 
  Paper, 
  TextField, 
  Button, 
  Tabs,
  Tab,
  Box,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  CircularProgress,
  Divider,
  makeStyles,
  useTheme
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import EqualizerIcon from '@material-ui/icons/Equalizer';
import AssessmentIcon from '@material-ui/icons/Assessment';
import MoleculeViewer3DImproved from './MoleculeViewer3DImproved';
import BindingAffinityBarChart from './BindingAffinityBarChart';
import AdmetRadarChart from './AdmetRadarChart';

// TabPanel component for tabbed interface
function TabPanel(props) {
  const { children, value, index, ...other } = props;

  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`simple-tabpanel-${index}`}
      aria-labelledby={`simple-tab-${index}`}
      {...other}
    >
      {value === index && (
        <Box p={3}>
          {children}
        </Box>
      )}
    </div>
  );
}

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    padding: theme.spacing(3),
  },
  title: {
    marginBottom: theme.spacing(4),
    fontWeight: 500,
  },
  paper: {
    padding: theme.spacing(3),
    height: '100%',
  },
  formContainer: {
    marginBottom: theme.spacing(3),
  },
  formControl: {
    minWidth: '100%',
    marginBottom: theme.spacing(2),
  },
  buttonWrapper: {
    position: 'relative',
    marginTop: theme.spacing(2),
  },
  buttonProgress: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12,
  },
  divider: {
    margin: theme.spacing(3, 0),
  },
  sectionTitle: {
    marginBottom: theme.spacing(2),
  },
  tabs: {
    marginBottom: theme.spacing(2),
  },
  resultsContainer: {
    marginTop: theme.spacing(3),
  },
  propertyItem: {
    display: 'flex',
    justifyContent: 'space-between',
    padding: theme.spacing(1, 0),
    borderBottom: `1px solid ${theme.palette.divider}`,
  },
  propertyLabel: {
    fontWeight: 500,
  },
  alertMargin: {
    marginBottom: theme.spacing(2),
  },
  emptyState: {
    textAlign: 'center',
    padding: theme.spacing(6),
  },
  emptyStateIcon: {
    fontSize: 64,
    color: theme.palette.text.disabled,
    marginBottom: theme.spacing(2),
  },
  resultsPaper: {
    padding: theme.spacing(3),
  },
  chartSection: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(2),
  },
  detailsSection: {
    marginTop: theme.spacing(2),
  }
}));

const SimulationPanel = () => {
  const classes = useStyles();
  const theme = useTheme();
  const [loading, setLoading] = useState(false);
  const [tabValue, setTabValue] = useState(0);
  const [simulationType, setSimulationType] = useState('binding');
  const [smiles, setSmiles] = useState('');
  const [simulationResults, setSimulationResults] = useState(null);
  const [error, setError] = useState(null);
  const [storedMolecules, setStoredMolecules] = useState([]);
  
  // Function to load molecules from localStorage
  const loadMoleculesFromStorage = () => {
      try {
          const moleculesString = localStorage.getItem('molecules') || '[]';
          console.log("Reloading molecules string from localStorage:", moleculesString);
          const loadedMolecules = JSON.parse(moleculesString);
          console.log("Parsed molecules on reload:", loadedMolecules);

          if (Array.isArray(loadedMolecules)) {
              const validMolecules = loadedMolecules.filter(m => m && m.id && m.name && m.smiles);
              console.log("Filtered valid molecules on reload:", validMolecules);
              setStoredMolecules(validMolecules);
              if (validMolecules.length === 0 && loadedMolecules.length > 0) {
                  console.warn("Some molecules were loaded but deemed invalid (missing id, name, or smiles).");
              }
          } else {
              console.error("Invalid data format found in localStorage for 'molecules': Expected an array.");
              setStoredMolecules([]);
          }
      } catch (e) {
          console.error("Error loading/parsing molecules from localStorage on reload:", e);
          setStoredMolecules([]);
          // Optionally set an error, but might be too noisy if storage clears often
          // setError("Could not load saved molecules. Check console for details.");
      }
  };

  // Initial load on mount
  useEffect(() => {
    loadMoleculesFromStorage();
  }, []); // Empty dependency array means run only once on mount

  // Listen for storage changes specifically for the 'molecules' key
  useEffect(() => {
      const handleStorageChange = (event) => {
          if (event.key === 'molecules') {
              console.log('Detected storage change for \'molecules\'. Reloading...');
              loadMoleculesFromStorage();
          }
      };

      window.addEventListener('storage', handleStorageChange);

      // Cleanup listener on component unmount
      return () => {
          window.removeEventListener('storage', handleStorageChange);
      };
  }, []); // Empty dependency array, listener setup once
  
  const handleTabChange = (event, newValue) => {
    setTabValue(newValue);
  };
  
  const handleSimulationTypeChange = (event) => {
    setSimulationType(event.target.value);
  };
  
  const handleSmilesChange = (event) => {
    setSmiles(event.target.value);
  };
  
  const handleRunSimulation = async () => {
    if (!smiles) {
      setError('Please select or enter a molecule SMILES notation');
      return;
    }

    setLoading(true);
    setError(null);
    setSimulationResults(null); // Clear previous results

    try {
      // --- Correct API Call Path ---
      const response = await fetch('/api/simulation/simulate', { // Corrected path
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles, simulationType }),
      });

      if (!response.ok) {
        // Handle HTTP errors
        let errorMsg = `Simulation failed with status: ${response.status}`;
        try {
          const errorData = await response.json();
          errorMsg = errorData.message || errorMsg; // Use backend error message if available
        } catch (e) {
          // Ignore if response body is not JSON
        }
        throw new Error(errorMsg);
      }

      const results = await response.json();

      // Validate response structure (basic check)
      if (!results || !results.type || !results.smiles || 
          (results.type === 'binding' && !results.bindingAffinities) || 
          (results.type === 'admet' && !results.admet)) {
          throw new Error('Received invalid simulation results from the server.');
      }
      
      // Set results from API
      console.log("Received simulation results from API:", JSON.stringify(results, null, 2)); // Log the raw results
      setSimulationResults(results);

      // Remove the mock logic below
      /*
      await new Promise(resolve => setTimeout(resolve, 2000));
      const moleculeName = storedMolecules.find(m => m.smiles === smiles)?.name || 'Custom Molecule';
      if (simulationType === 'binding') {
        setSimulationResults({
          type: 'binding',
          moleculeName: moleculeName,
          smiles: smiles,
          bindingAffinities: {
            'Dopamine Transporter': { score: 84, classification: 'Strong' },
            'Norepinephrine Transporter': { score: 72, classification: 'Moderate' },
            'Serotonin Transporter': { score: 35, classification: 'Weak' },
            'D1 Receptor': { score: 22, classification: 'Weak' },
            'D2 Receptor': { score: 45, classification: 'Moderate' }
          },
          properties: {
            molecularWeight: 165.23,
            logP: 1.86,
            tpsa: 29.54,
            hDonors: 1,
            hAcceptors: 2,
            rotatableBonds: 3,
            lipinskiViolations: 0
          }
        });
      } else if (simulationType === 'admet') {
        setSimulationResults({
          type: 'admet',
          moleculeName: moleculeName,
          smiles: smiles,
          admet: {
            absorption: {
              oral: { score: 85, classification: 'High' },
              bbb: { score: 62, classification: 'Medium' },
              pgpSubstrate: { score: 32, classification: 'Unlikely' },
              f_oral: { percent: 78, classification: 'High' }
            },
            distribution: {
              vd: { vd: 1.8, classification: 'Moderate' },
              plasmaProteinBinding: { percent: 75.5, classification: 'Moderate' },
              tissueDistribution: {
                brain: 'Medium',
                adipose: 'Low',
                liver: 'High',
                kidney: 'Moderate'
              }
            },
            metabolism: {
              cyp450Substrates: {
                CYP3A4: { score: 72, isSubstrate: true },
                CYP2D6: { score: 45, isSubstrate: false },
                CYP2C9: { score: 38, isSubstrate: false }
              },
              cyp450Inhibition: {
                CYP3A4: { score: 35, inhibition: 'Weak' },
                CYP2D6: { score: 22, inhibition: 'Weak' },
                CYP2C9: { score: 18, inhibition: 'Weak' }
              },
              halfLife: { hours: 5.6, classification: 'Moderate' }
            },
            excretion: {
              renalClearance: { score: 65, classification: 'Moderate' }
            },
            toxicity: {
              herg: { score: 28, risk: 'Low' },
              hepatotoxicity: { score: 35, risk: 'Low' },
              cardiotoxicity: { score: 25, risk: 'Low' },
              mutagenicity: { score: 18, risk: 'Low' }
            }
          },
          properties: {
            molecularWeight: 165.23,
            logP: 1.86,
            tpsa: 29.54,
            hDonors: 1,
            hAcceptors: 2,
            rotatableBonds: 3,
            lipinskiViolations: 0
          }
        });
      }
      */
    } catch (err) {
      console.error('Simulation API call failed:', err);
      // Use the error message from the thrown Error object
      setError(err.message || 'An error occurred while running the simulation. Please check the console and backend logs.');
      setSimulationResults(null); // Ensure results are cleared on error
    } finally {
      setLoading(false);
    }
  };
  
  const handleClearResults = () => {
    setSimulationResults(null);
  };

  return (
    <div className={classes.root}>
      <Typography variant="h4" className={classes.title}>
        Simulation Panel
      </Typography>
      
      <Grid container spacing={3}>
        <Grid item xs={12} md={4}>
          <Paper className={classes.paper}>
            <Typography variant="h6" gutterBottom>
              Simulation Parameters
            </Typography>
            
            <Tabs
              value={tabValue}
              onChange={handleTabChange}
              indicatorColor="primary"
              textColor="primary"
              variant="fullWidth"
              className={classes.tabs}
            >
              <Tab label="Select Molecule" />
              <Tab label="Enter SMILES" />
            </Tabs>
            
            <TabPanel value={tabValue} index={0}>
              <FormControl variant="outlined" className={classes.formControl}>
                <InputLabel id="molecule-select-label">Select Molecule</InputLabel>
                <Select
                  labelId="molecule-select-label"
                  id="molecule-select"
                  value={smiles}
                  onChange={handleSmilesChange}
                  label="Select Molecule"
                >
                  <MenuItem value="">
                    <em>Select a molecule</em>
                  </MenuItem>
                  {storedMolecules.map((molecule) => (
                    <MenuItem key={molecule.id} value={molecule.smiles}>
                      {molecule.name}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            </TabPanel>
            
            <TabPanel value={tabValue} index={1}>
              <TextField
                label="SMILES Notation"
                variant="outlined"
                fullWidth
                value={smiles}
                onChange={handleSmilesChange}
                placeholder="e.g., CC(N)CC1=CC=CC=C1"
                className={classes.formControl}
              />
            </TabPanel>
            
            <Divider className={classes.divider} />
            
            <Typography variant="subtitle1" gutterBottom>
              Simulation Type
            </Typography>
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel id="simulation-type-label">Simulation Type</InputLabel>
              <Select
                labelId="simulation-type-label"
                id="simulation-type"
                value={simulationType}
                onChange={handleSimulationTypeChange}
                label="Simulation Type"
              >
                <MenuItem value="binding">Receptor Binding Affinity</MenuItem>
                <MenuItem value="admet">ADMET Properties</MenuItem>
              </Select>
            </FormControl>
            
            {error && (
              <Alert severity="error" className={classes.alertMargin}>
                {error}
              </Alert>
            )}
            
            <div className={classes.buttonWrapper}>
              <Button
                variant="contained"
                color="primary"
                fullWidth
                onClick={handleRunSimulation}
                disabled={loading}
              >
                Run Simulation
              </Button>
              {loading && <CircularProgress size={24} className={classes.buttonProgress} />}
            </div>
            
            {simulationResults && (
              <Button
                variant="outlined"
                color="primary"
                fullWidth
                onClick={handleClearResults}
                style={{ marginTop: 16 }}
              >
                Clear Results
              </Button>
            )}
          </Paper>
        </Grid>
        
        <Grid item xs={12} md={8}>
          {loading ? (
            <Paper className={classes.paper} style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: 500 }}>
              <CircularProgress />
            </Paper>
          ) : simulationResults ? (
            <Paper className={classes.resultsPaper}>
              <Typography variant="h6" gutterBottom>
                {simulationResults.moleculeName} - {simulationType === 'binding' ? 'Receptor Binding Affinity' : 'ADMET Properties'} Results
              </Typography>
              
              <Grid container spacing={3} style={{ marginBottom: '16px' }}>
                <Grid item xs={12} md={6}>
                  <MoleculeViewer3DImproved smiles={simulationResults.smiles} height={300} />
                </Grid>
                
                <Grid item xs={12} md={6}>
                  <Typography variant="subtitle1" gutterBottom>
                    Molecular Properties
                  </Typography>
                  <Grid container spacing={0}>
                    {simulationResults.properties && Object.entries(simulationResults.properties)
                      .filter(([key, value]) => value !== null && value !== undefined && key !== 'png_base64')
                      .map(([key, value]) => (
                        <Grid item xs={12} key={key}>
                          <div className={classes.propertyItem} style={{ padding: '4px 0', alignItems: 'baseline'}}>
                            <Typography className={classes.propertyLabel} variant="body2" style={{ minWidth: '120px', flexShrink: 0, marginRight: '8px' }}>
                              { {logP: 'LogP', tpsa: 'TPSA', hDonors: 'H-Bond Donors', hAcceptors: 'H-Bond Acceptors', 
                                 lipinskiViolations: 'Lipinski Violations', rotatableBonds: 'Rotatable Bonds', 
                                 molecularWeight: 'Molecular Weight'}[key] || 
                                key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
                            </Typography>
                            <Typography variant="body2" style={{ wordBreak: 'break-all' }}>
                              {typeof value === 'string' && value.length > 100 ? `${value.substring(0, 100)}...` : 
                               (typeof value === 'number' && !['lipinskiViolations', 'rings', 'heavy_atoms', 'hba', 'hbd', 'rotatable_bonds'].includes(key) ? value.toFixed(2) : value)} 
                              {key === 'molecularWeight' ? ' g/mol' : 
                               key === 'tpsa' ? ' Å²' : ''}
                            </Typography>
                          </div>
                        </Grid>
                    ))}
                  </Grid>
                </Grid>
              </Grid>
              
              <Divider className={classes.divider} />
              
              {simulationResults.type === 'binding' && (
                <>
                  <div className={classes.chartSection}>
                    <Typography variant="subtitle1" gutterBottom>
                      Receptor Binding Affinities
                    </Typography>
                    <BindingAffinityBarChart data={simulationResults.bindingAffinities} />
                  </div>
                  
                  <div className={classes.detailsSection}>
                    <Grid container spacing={2}>
                      {Object.entries(simulationResults.bindingAffinities ?? {}).map(([receptor, value]) => (
                        <Grid item xs={12} sm={6} md={4} key={receptor}>
                          <Paper elevation={0} style={{ padding: '8px', border: `1px solid ${theme.palette.divider}`, height: '100%'}}>
                            <Typography className={classes.propertyLabel} display="block" variant="body2">
                              {receptor}
                            </Typography>
                            <Typography variant="body2">
                              {value.score ?? 'N/A'}/100 ({value.classification ?? 'N/A'})
                            </Typography>
                          </Paper>
                        </Grid>
                      ))}
                    </Grid>
                  </div>
                </>
              )}
              
              {simulationResults.type === 'admet' && (
                <>
                  <div className={classes.chartSection}>
                    <Typography variant="subtitle1" gutterBottom>
                      ADMET Properties
                    </Typography>
                    <AdmetRadarChart data={simulationResults.admet} />
                  </div>
                  
                  <div className={classes.detailsSection}>
                    <Typography variant="subtitle1" gutterBottom style={{ marginBottom: '12px' }}>
                      Detailed ADMET Parameters
                    </Typography>
                    <Grid container spacing={3}>
                      <Grid item xs={12} md={6}>
                        <Typography variant="subtitle2" gutterBottom>
                          Absorption
                        </Typography>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Oral Absorption</Typography>
                          <Typography>{simulationResults.admet?.absorption?.oral?.classification ?? 'N/A'} ({simulationResults.admet?.absorption?.oral?.score ?? 'N/A'}/100)</Typography>
                        </div>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Blood-Brain Barrier</Typography>
                          <Typography>{simulationResults.admet?.absorption?.bbb?.classification ?? 'N/A'} ({simulationResults.admet?.absorption?.bbb?.score ?? 'N/A'}/100)</Typography>
                        </div>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>P-gp Substrate</Typography>
                          <Typography>{simulationResults.admet?.absorption?.pgpSubstrate?.classification ?? 'N/A'}</Typography>
                        </div>
                        
                        <Typography variant="subtitle2" gutterBottom style={{ marginTop: 24 }}>
                          Distribution
                        </Typography>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Volume of Distribution</Typography>
                          <Typography>{simulationResults.admet?.distribution?.vd?.vd ?? 'N/A'} L/kg ({simulationResults.admet?.distribution?.vd?.classification ?? 'N/A'})</Typography>
                        </div>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Plasma Protein Binding</Typography>
                          <Typography>{simulationResults.admet?.distribution?.plasmaProteinBinding?.percent ?? 'N/A'}%</Typography>
                        </div>
                      </Grid>
                      
                      <Grid item xs={12} md={6}>
                        <Typography variant="subtitle2" gutterBottom>
                          Metabolism
                        </Typography>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>CYP3A4 Substrate</Typography>
                          <Typography>{simulationResults.admet?.metabolism?.cyp450Substrates?.CYP3A4?.isSubstrate ? 'Yes' : (simulationResults.admet?.metabolism?.cyp450Substrates?.CYP3A4 !== undefined ? 'No' : 'N/A')}</Typography>
                        </div>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Half-life</Typography>
                          <Typography>{simulationResults.admet?.metabolism?.halfLife?.hours ?? 'N/A'} hours ({simulationResults.admet?.metabolism?.halfLife?.classification ?? 'N/A'})</Typography>
                        </div>
                        
                        <Typography variant="subtitle2" gutterBottom style={{ marginTop: 24 }}>
                          Toxicity
                        </Typography>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>hERG Inhibition</Typography>
                          <Typography>{simulationResults.admet?.toxicity?.herg?.risk ?? 'N/A'} Risk</Typography>
                        </div>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Hepatotoxicity</Typography>
                          <Typography>{simulationResults.admet?.toxicity?.hepatotoxicity?.risk ?? 'N/A'} Risk</Typography>
                        </div>
                        
                        <div className={classes.propertyItem}>
                          <Typography className={classes.propertyLabel}>Cardiotoxicity</Typography>
                          <Typography>{simulationResults.admet?.toxicity?.cardiotoxicity?.risk ?? 'N/A'} Risk</Typography>
                        </div>
                      </Grid>
                    </Grid>
                  </div>
                </>
              )}
              
              <Divider className={classes.divider} />
              
              <Grid container spacing={2}>
                <Grid item xs={12} sm={6}>
                  <Button variant="outlined" color="primary" fullWidth>
                    Export Results
                  </Button>
                </Grid>
                <Grid item xs={12} sm={6}>
                  <Button variant="contained" color="primary" fullWidth>
                    Proceed to Regulatory Analysis
                  </Button>
                </Grid>
              </Grid>
            </Paper>
          ) : (
            <Paper className={classes.paper} style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', minHeight: 500, flexDirection: 'column' }}>
              <AssessmentIcon className={classes.emptyStateIcon} />
              <Typography variant="h6" color="textSecondary">
                No simulation results
              </Typography>
              <Typography variant="body2" color="textSecondary">
                Select a molecule and simulation type, then click "Run Simulation"
              </Typography>
            </Paper>
          )}
        </Grid>
      </Grid>
    </div>
  );
};

export default SimulationPanel; 