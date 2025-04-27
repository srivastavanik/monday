import React, { useState, useEffect, useRef } from 'react';
import { 
  Typography, 
  Grid, 
  Paper, 
  FormControl, 
  InputLabel, 
  Select, 
  MenuItem, 
  Button, 
  Divider, 
  CircularProgress, 
  Chip,
  makeStyles 
} from '@material-ui/core';
import CompareArrowsIcon from '@material-ui/icons/CompareArrows';
import DeleteIcon from '@material-ui/icons/Delete';
import Alert from '@material-ui/lab/Alert';
import { simulationAPI } from '../services/api';

// This would be imported from a third-party library in a real implementation
const MoleculeViewer = ({ smiles, height = 250 }) => {
  const viewerRef = useRef(null);
  
  useEffect(() => {
    if (smiles && viewerRef.current) {
      // In a real implementation, this would initialize 3Dmol.js
      const placeholder = document.createElement('div');
      placeholder.style.width = '100%';
      placeholder.style.height = `${height}px`;
      placeholder.style.backgroundColor = '#f5f5f5';
      placeholder.style.display = 'flex';
      placeholder.style.alignItems = 'center';
      placeholder.style.justifyContent = 'center';
      
      const text = document.createElement('div');
      text.innerHTML = `<strong>SMILES:</strong> ${smiles}<br><br>3D Molecule Viewer would render here`;
      text.style.textAlign = 'center';
      text.style.padding = '20px';
      
      placeholder.appendChild(text);
      
      viewerRef.current.innerHTML = '';
      viewerRef.current.appendChild(placeholder);
    }
  }, [smiles, height]);
  
  return <div ref={viewerRef} style={{ width: '100%', height: `${height}px`, border: '1px solid #e0e0e0' }}></div>;
};

// This would be imported from a third-party library in a real implementation
const SimilarityChart = ({ molecules, similarity }) => {
  const chartRef = useRef(null);
  
  useEffect(() => {
    if (molecules && similarity && chartRef.current) {
      // In a real implementation, this would initialize a Plotly or D3 chart
      const placeholder = document.createElement('div');
      placeholder.style.width = '100%';
      placeholder.style.height = '300px';
      placeholder.style.backgroundColor = '#f5f5f5';
      placeholder.style.display = 'flex';
      placeholder.style.alignItems = 'center';
      placeholder.style.justifyContent = 'center';
      placeholder.style.flexDirection = 'column';
      
      const text = document.createElement('div');
      text.innerHTML = '<strong>Tanimoto Similarity</strong><br><br>';
      text.style.textAlign = 'center';
      text.style.padding = '20px';
      
      const similarityText = document.createElement('div');
      similarityText.style.fontSize = '2rem';
      similarityText.style.fontWeight = 'bold';
      similarityText.innerHTML = `${(similarity * 100).toFixed(1)}%`;
      
      const description = document.createElement('div');
      description.style.marginTop = '16px';
      description.innerHTML = `Similarity between ${molecules[0].name} and ${molecules[1].name}<br>Based on Morgan fingerprint Tanimoto coefficient`;
      
      placeholder.appendChild(text);
      placeholder.appendChild(similarityText);
      placeholder.appendChild(description);
      
      chartRef.current.innerHTML = '';
      chartRef.current.appendChild(placeholder);
    }
  }, [molecules, similarity]);
  
  return <div ref={chartRef} style={{ width: '100%', height: '300px', border: '1px solid #e0e0e0' }}></div>;
};

// This would be imported from a third-party library in a real implementation
const PropertyComparisonChart = ({ data }) => {
  const chartRef = useRef(null);
  
  useEffect(() => {
    if (data && chartRef.current) {
      // In a real implementation, this would initialize a Plotly or D3 chart
      const placeholder = document.createElement('div');
      placeholder.style.width = '100%';
      placeholder.style.height = '400px';
      placeholder.style.backgroundColor = '#f5f5f5';
      placeholder.style.display = 'flex';
      placeholder.style.alignItems = 'center';
      placeholder.style.justifyContent = 'center';
      
      const text = document.createElement('div');
      text.innerHTML = '<strong>Property Comparison Chart</strong><br><br>Bar chart comparing properties would render here';
      text.style.textAlign = 'center';
      text.style.padding = '20px';
      
      placeholder.appendChild(text);
      
      chartRef.current.innerHTML = '';
      chartRef.current.appendChild(placeholder);
    }
  }, [data]);
  
  return <div ref={chartRef} style={{ width: '100%', height: '400px', border: '1px solid #e0e0e0' }}></div>;
};

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
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120,
    width: '100%',
  },
  divider: {
    margin: theme.spacing(3, 0),
  },
  moleculeSelector: {
    display: 'flex',
    alignItems: 'center',
    marginBottom: theme.spacing(2),
  },
  compareButton: {
    marginTop: theme.spacing(2),
  },
  progress: {
    display: 'flex',
    justifyContent: 'center',
    margin: theme.spacing(4, 0),
  },
  propertyContainer: {
    marginTop: theme.spacing(2),
  },
  propertyRow: {
    padding: theme.spacing(1, 0),
    borderBottom: `1px solid ${theme.palette.divider}`,
    display: 'flex',
  },
  propertyLabel: {
    fontWeight: 500,
    flex: 1,
  },
  propertyValue: {
    flex: 1,
    textAlign: 'center',
  },
  better: {
    color: theme.palette.success.main,
    fontWeight: 500,
  },
  worse: {
    color: theme.palette.error.main,
    fontWeight: 500,
  },
  neutral: {
    color: theme.palette.text.primary,
  },
  sectionTitle: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(2),
  },
  chip: {
    margin: theme.spacing(0.5),
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
  deleteButton: {
    marginLeft: theme.spacing(1),
  },
  alertMargin: {
    marginBottom: theme.spacing(2),
  },
}));

const ComparisonTool = () => {
  const classes = useStyles();
  const [loading, setLoading] = useState(false);
  const [selectedMolecules, setSelectedMolecules] = useState([null, null]);
  const [comparisonResults, setComparisonResults] = useState(null);
  const [error, setError] = useState(null);
  
  const [availableMolecules, setAvailableMolecules] = useState([
    { id: 1, name: 'Methylphenidate', smiles: 'CN(C)C(C1=CC=CC=C1)C(C)OC(=O)C' },
    { id: 2, name: 'Amphetamine', smiles: 'CC(N)CC1=CC=CC=C1' },
    { id: 3, name: 'Atomoxetine', smiles: 'CC(C)NCC1=CC=CC(OC2=CC=CC=C2)=C1' },
    { id: 4, name: 'Novel Amphetamine Derivative', smiles: 'CC(CC1=CC=C(C=C1)O)NC' },
    { id: 5, name: 'Dextroamphetamine', smiles: 'CC(N)CC1=CC=CC=C1' },
    { id: 6, name: 'Lisdexamfetamine', smiles: 'NC(=O)C(N)CCCCN.CC(N)CC1=CC=CC=C1' }
  ]);
  
  const handleMoleculeChange = (index, event) => {
    const value = event.target.value;
    const newSelectedMolecules = [...selectedMolecules];
    newSelectedMolecules[index] = value ? availableMolecules.find(m => m.id === value) : null;
    setSelectedMolecules(newSelectedMolecules);
    setComparisonResults(null);
    setError(null);
  };
  
  const handleClearMolecule = (index) => {
    const newSelectedMolecules = [...selectedMolecules];
    newSelectedMolecules[index] = null;
    setSelectedMolecules(newSelectedMolecules);
    setComparisonResults(null);
    setError(null);
  };
  
  const handleCompare = async () => {
    if (!selectedMolecules[0]?.smiles || !selectedMolecules[1]?.smiles) {
      setError('Please select two molecules to compare');
      return;
    }
    
    setLoading(true);
    setError(null);
    setComparisonResults(null);
    
    try {
      const response = await simulationAPI.compareMolecules(
        selectedMolecules[0].smiles,
        selectedMolecules[1].smiles
      );
      
      if (response.data && !response.data.error) {
        setComparisonResults({
          molecules: selectedMolecules,
          similarity: response.data.similarity,
          properties1: response.data.properties1,
          properties2: response.data.properties2
        });
      } else {
        throw new Error(response.data.error || 'Failed to compare molecules');
      }
    } catch (err) {
      console.error('Error comparing molecules:', err);
      setError(err.response?.data?.error || 'An error occurred during comparison. Please try again.');
    } finally {
      setLoading(false);
    }
  };
  
  const handleClearComparison = () => {
    setComparisonResults(null);
    setError(null);
  };
  
  const formatPropertyName = (key) => {
    return key
        .replace(/([A-Z])/g, ' $1')
        .replace(/^./, (str) => str.toUpperCase());
  };
  
  const getPropertyUnits = (key) => {
    if (key === 'molWeight') return ' g/mol';
    if (key === 'tpsa') return ' Å²';
    return '';
  };
  
  const getComparisonClass = (key, value1, value2) => {
    if (typeof value1 !== 'number' || typeof value2 !== 'number') return classes.neutral;
    const difference = value2 - value1;

    if (key === 'qed') {
      return difference > 0.01 ? classes.better : difference < -0.01 ? classes.worse : classes.neutral;
    }
    if (['logP', 'tpsa', 'molWeight', 'rotatableBondCount', 'lipinskiViolations'].includes(key)) {
      return difference < -0.01 ? classes.better : difference > 0.01 ? classes.worse : classes.neutral;
    }
    return classes.neutral;
  };

  return (
    <div className={classes.root}>
      <Typography variant="h4" className={classes.title}>
        Molecule Comparison Tool
      </Typography>
      
      <Paper className={classes.paper}>
        <Typography variant="h6" gutterBottom>
          Select Molecules to Compare
        </Typography>
        
        <Grid container spacing={3}>
          <Grid item xs={12} md={5}>
            <div className={classes.moleculeSelector}>
              <FormControl variant="outlined" className={classes.formControl}>
                <InputLabel id="molecule1-label">Molecule 1</InputLabel>
                <Select
                  labelId="molecule1-label"
                  value={selectedMolecules[0]?.id || ''}
                  onChange={(e) => handleMoleculeChange(0, e)}
                  label="Molecule 1"
                >
                  <MenuItem value="">
                    <em>Select a molecule</em>
                  </MenuItem>
                  {availableMolecules.map((molecule) => (
                    <MenuItem 
                      key={molecule.id} 
                      value={molecule.id}
                      disabled={selectedMolecules[1]?.id === molecule.id}
                    >
                      {molecule.name} ({molecule.smiles.substring(0, 15)}...)
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
              
              {selectedMolecules[0] && (
                <Button 
                  className={classes.deleteButton}
                  onClick={() => handleClearMolecule(0)}
                  aria-label="clear molecule 1"
                >
                  <DeleteIcon fontSize="small" />
                </Button>
              )}
            </div>
            
            {selectedMolecules[0] && (
              <Typography variant="caption">{selectedMolecules[0].smiles}</Typography>
            )}
          </Grid>
          
          <Grid item xs={12} md={2} style={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
            <CompareArrowsIcon style={{ fontSize: 48, color: '#757575' }} />
          </Grid>
          
          <Grid item xs={12} md={5}>
            <div className={classes.moleculeSelector}>
              <FormControl variant="outlined" className={classes.formControl}>
                <InputLabel id="molecule2-label">Molecule 2</InputLabel>
                <Select
                  labelId="molecule2-label"
                  value={selectedMolecules[1]?.id || ''}
                  onChange={(e) => handleMoleculeChange(1, e)}
                  label="Molecule 2"
                >
                  <MenuItem value="">
                    <em>Select a molecule</em>
                  </MenuItem>
                  {availableMolecules.map((molecule) => (
                    <MenuItem 
                      key={molecule.id} 
                      value={molecule.id}
                      disabled={selectedMolecules[0]?.id === molecule.id}
                    >
                      {molecule.name} ({molecule.smiles.substring(0, 15)}...)
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
              
              {selectedMolecules[1] && (
                <Button 
                  className={classes.deleteButton}
                  onClick={() => handleClearMolecule(1)}
                  aria-label="clear molecule 2"
                >
                  <DeleteIcon fontSize="small" />
                </Button>
              )}
            </div>
            
            {selectedMolecules[1] && (
              <Typography variant="caption">{selectedMolecules[1].smiles}</Typography>
            )}
          </Grid>
        </Grid>
        
        {error && (
          <Alert severity="error" className={classes.alertMargin} style={{ marginTop: 16 }}>
            {error}
          </Alert>
        )}
        
        <div style={{ textAlign: 'center', marginTop: 24 }}>
          {!comparisonResults ? (
            <Button
              variant="contained"
              color="primary"
              className={classes.compareButton}
              onClick={handleCompare}
              disabled={!selectedMolecules[0] || !selectedMolecules[1] || loading}
            >
              Compare Molecules
            </Button>
          ) : (
            <Button
              variant="outlined"
              color="primary"
              className={classes.compareButton}
              onClick={handleClearComparison}
            >
              Clear Comparison
            </Button>
          )}
        </div>
        
        {loading && (
          <div className={classes.progress}>
            <CircularProgress />
          </div>
        )}
        
        {comparisonResults && !loading && (
          <>
            <Divider className={classes.divider} />
            
            <Typography variant="h6" gutterBottom>
              Comparison Results: {comparisonResults.molecules[0]?.name || 'Molecule 1'} vs. {comparisonResults.molecules[1]?.name || 'Molecule 2'}
            </Typography>
            
            <Grid container spacing={3}>
              <Grid item xs={12} md={6}>
                <Typography variant="subtitle1" gutterBottom>Structural Similarity</Typography>
                <Paper variant="outlined" style={{ padding: 16, textAlign: 'center' }}>
                  <Typography variant="h4">{(comparisonResults.similarity.tanimoto * 100).toFixed(1)}%</Typography>
                  <Typography variant="caption">Tanimoto Similarity (Morgan Fingerprints)</Typography>
                  {comparisonResults.similarity.mcs_smarts && (
                    <Typography variant="body2" style={{ marginTop: 8 }}>
                      Max Common Substructure (MCS): 
                      <Chip size="small" label={`${comparisonResults.similarity.mcs_atom_count} atoms, ${comparisonResults.similarity.mcs_bond_count} bonds`} style={{marginLeft: 4}}/> <br/>
                      <code style={{fontSize: '0.8em', wordBreak: 'break-all'}}>{comparisonResults.similarity.mcs_smarts}</code>
                    </Typography>
                  )}
                </Paper>
              </Grid>
            </Grid>
            
            <Typography variant="subtitle1" className={classes.sectionTitle}>
              Detailed Property Comparison
            </Typography>
            
            <div className={classes.propertyContainer}>
              <div className={classes.propertyRow} style={{ fontWeight: 500, backgroundColor: '#f5f5f5' }}>
                <div className={classes.propertyLabel}>Property</div>
                <div className={classes.propertyValue}>{comparisonResults.molecules[0]?.name || 'Molecule 1'}</div>
                <div className={classes.propertyValue}>{comparisonResults.molecules[1]?.name || 'Molecule 2'}</div>
                <div className={classes.propertyValue}>Difference</div>
              </div>
              
              {comparisonResults.properties1 && Object.keys(comparisonResults.properties1).map((key) => {
                const val1 = comparisonResults.properties1[key];
                const val2 = comparisonResults.properties2 ? comparisonResults.properties2[key] : undefined;
                const diff = (typeof val1 === 'number' && typeof val2 === 'number') ? (val2 - val1) : 'N/A';
                const units = getPropertyUnits(key);
                
                return (
                  <div className={classes.propertyRow} key={key}>
                    <div className={classes.propertyLabel}>{formatPropertyName(key)}</div>
                    <div className={classes.propertyValue}>{(typeof val1 === 'number' ? val1.toFixed(2) : val1)}{units}</div>
                    <div className={classes.propertyValue}>{(typeof val2 === 'number' ? val2.toFixed(2) : val2)}{units}</div>
                    <div className={`${classes.propertyValue} ${getComparisonClass(key, val1, val2)}`}>
                      {typeof diff === 'number' ? `${diff > 0 ? '+' : ''}${diff.toFixed(2)}` : diff}
                    </div>
                  </div>
                );
              })}
            </div>
            
            <Divider className={classes.divider} />
            
            <Grid container spacing={2}>
              <Grid item xs={12} sm={6}>
                <Button variant="outlined" color="primary" fullWidth>
                  Export Comparison
                </Button>
              </Grid>
              <Grid item xs={12} sm={6}>
                <Button variant="contained" color="primary" fullWidth>
                  Add to Report
                </Button>
              </Grid>
            </Grid>
          </>
        )}
        
        {!comparisonResults && !loading && (
          <div className={classes.emptyState}>
            <CompareArrowsIcon className={classes.emptyStateIcon} />
            <Typography variant="h6" color="textSecondary">
              Select molecules to compare
            </Typography>
            <Typography variant="body2" color="textSecondary">
              Choose two molecules and click "Compare Molecules" to see a detailed comparison
            </Typography>
          </div>
        )}
      </Paper>
    </div>
  );
};

export default ComparisonTool; 