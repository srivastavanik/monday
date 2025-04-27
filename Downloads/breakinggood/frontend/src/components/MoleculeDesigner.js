import React, { useState, useEffect, useRef } from 'react';
import { 
  Typography, 
  Grid, 
  Paper, 
  TextField, 
  Button, 
  Divider, 
  CircularProgress,
  Tabs,
  Tab,
  Box,
  Chip,
  FormControlLabel,
  Checkbox,
  makeStyles,
  Card,
  CardContent,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Snackbar,
  IconButton,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  LinearProgress,
  CardHeader
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import BuildIcon from '@material-ui/icons/Build'; // Using Build icon instead of Science
import FormatQuoteIcon from '@material-ui/icons/FormatQuote';
import SaveIcon from '@material-ui/icons/Save';
import DeleteIcon from '@material-ui/icons/Delete';
import EditIcon from '@material-ui/icons/Edit';
import CompareArrowsIcon from '@material-ui/icons/CompareArrows';
import AddCircleOutlineIcon from '@material-ui/icons/AddCircleOutline';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import CloudDownloadIcon from '@material-ui/icons/CloudDownload';
import PsychologyAltIcon from '@material-ui/icons/EmojiObjects'; // Using EmojiObjects icon instead of Psychology
import DescriptionIcon from '@material-ui/icons/Description';
import { ToggleButton, ToggleButtonGroup } from '@material-ui/lab';
import MoleculeViewer3D from './MoleculeViewer3DImproved';
import AIChatInterfaceRedesigned from './AIChatInterfaceRedesigned'; // Import the REDESIGNED chat interface
import { drugDesignAPI, simulationAPI, claudeAPI } from '../services/api';
import StructureEditor from './StructureEditor'; // Import the actual editor
import ReactMarkdown from 'react-markdown'; // Make sure this is imported

// This would be imported from a third-party library in a real implementation
const MoleculeViewer = ({ smiles, height = 400 }) => {
  const viewerRef = useRef(null);
  
  useEffect(() => {
    if (smiles && viewerRef.current) {
      // In a real implementation, this would initialize 3Dmol.js
      // and render the molecule
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
  inputField: {
    marginBottom: theme.spacing(2),
  },
  divider: {
    margin: theme.spacing(2, 0),
  },
  viewerContainer: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(3),
  },
  resultContainer: {
    marginTop: theme.spacing(3),
  },
  chip: {
    margin: theme.spacing(0.5),
  },
  progress: {
    display: 'flex',
    justifyContent: 'center',
    marginTop: theme.spacing(4),
    marginBottom: theme.spacing(4),
  },
  section: {
    marginBottom: theme.spacing(3),
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
  referenceItem: {
    padding: theme.spacing(1),
    marginBottom: theme.spacing(1),
    backgroundColor: '#f9f9f9',
  },
  tabs: {
    marginBottom: theme.spacing(2),
  },
  toggleContainer: {
    margin: theme.spacing(2, 0),
  },
  buttonProgress: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12,
  },
  buttonWrapper: {
    position: 'relative',
    display: 'inline-block',
  },
  formControl: {
    marginBottom: theme.spacing(2),
    minWidth: '100%',
  },
  button: {
    marginTop: theme.spacing(2),
  },
  canvasContainer: {
    height: 400,
    border: '1px solid #ddd',
    borderRadius: theme.shape.borderRadius,
    backgroundColor: '#f9f9f9',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    marginTop: theme.spacing(2),
  },
  loadingContainer: {
    display: 'flex',
    flexDirection: 'column',
    justifyContent: 'center',
    alignItems: 'center',
    padding: theme.spacing(3),
  },
  resultCard: {
    marginTop: theme.spacing(2),
    position: 'relative',
  },
  propertyGrid: {
    marginTop: theme.spacing(2),
  },
  actionButtons: {
    position: 'absolute',
    top: theme.spacing(1),
    right: theme.spacing(1),
    display: 'flex',
  },
  tableContainer: {
    marginTop: theme.spacing(2),
    maxHeight: 440,
  },
  tableRow: {
    cursor: 'pointer',
    '&:hover': {
      backgroundColor: theme.palette.action.hover,
    },
  },
  editorContainer: {
    border: '1px solid #ddd',
    borderRadius: theme.shape.borderRadius,
    height: 400,
    marginTop: theme.spacing(2),
    backgroundColor: '#f9f9f9',
  },
  structureButtons: {
    display: 'flex',
    justifyContent: 'space-between',
    marginTop: theme.spacing(2),
  },
  uploadInput: {
    display: 'none',
  },
  buttonIconMargin: {
    marginRight: theme.spacing(1),
  },
  actionsContainer: {
    display: 'flex',
    justifyContent: 'flex-end',
    marginTop: theme.spacing(2),
  },
  aiGenerationCard: {
    marginTop: theme.spacing(3),
    padding: theme.spacing(2),
    backgroundColor: '#f0f7ff',
  },
  thinkingContainer: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(3),
  },
  aiControlsContainer: {
    display: 'flex',
    flexDirection: 'column',
    gap: theme.spacing(2),
    marginBottom: theme.spacing(2),
  },
  aiParameterCard: {
    padding: theme.spacing(2),
    marginBottom: theme.spacing(2),
  },
  generationRequirements: {
    marginBottom: theme.spacing(2),
  },
  apiResponseSection: {
    width: '100%',
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(3),
  },
  apiResponseTitle: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText,
    '& .MuiCardHeader-title': {
      fontSize: '1.25rem',
      fontWeight: 500,
    },
  },
  apiResponseContent: {
    fontSize: '1rem', 
    lineHeight: 1.6,
    '& h1, & h2, & h3': {
      color: theme.palette.primary.dark,
      marginTop: theme.spacing(3),
      marginBottom: theme.spacing(1),
    },
    '& p': {
      marginBottom: theme.spacing(2),
    },
    '& ul, & ol': {
      paddingLeft: theme.spacing(3),
      marginBottom: theme.spacing(2),
    },
    '& code': {
      fontFamily: '"Roboto Mono", monospace',
      backgroundColor: 'rgba(0, 0, 0, 0.04)',
      padding: theme.spacing(0.5),
      borderRadius: 4,
    },
  },
  fullWidthApiResponseContainer: {
    marginTop: theme.spacing(3),
  },
  fullWidthApiResponsePaper: {
    padding: theme.spacing(2),
    maxHeight: '600px',
    overflowY: 'auto',
  },
  apiResponseGridItem: {
    marginTop: theme.spacing(3),
  },
  apiResponseDivider: {
    marginBottom: theme.spacing(2),
  },
  stickyContainer: {
    position: 'sticky',
    top: 0,
    backgroundColor: theme.palette.background.default,
  },
}));

const MoleculeDesigner = () => {
  const classes = useStyles();
  const [activeTab, setActiveTab] = useState(0);
  const [loading, setLoading] = useState(false);
  const [structureInput, setStructureInput] = useState('');
  const [targetReceptor, setTargetReceptor] = useState('');
  const [generatedMolecules, setGeneratedMolecules] = useState([]);
  const [savedMolecules, setSavedMolecules] = useState([]);
  const [error, setError] = useState(null);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [snackbarMessage, setSnackbarMessage] = useState('');
  const [snackbarSeverity, setSnackbarSeverity] = useState('success');
  const [editDialogOpen, setEditDialogOpen] = useState(false);
  const [editMolecule, setEditMolecule] = useState(null);
  const [moleculeName, setMoleculeName] = useState('');
  const [smilesStructure, setSmilesStructure] = useState('');
  const [simulationDialogOpen, setSimulationDialogOpen] = useState(false);
  const [simulationType, setSimulationType] = useState('docking');
  
  // AI generation states
  const [aiRequirements, setAiRequirements] = useState('');
  const [aiTargetReceptors, setAiTargetReceptors] = useState('dopamine/norepinephrine transporters');
  const [aiMinimizeSideEffects, setAiMinimizeSideEffects] = useState(true);
  const [aiIncludeReferences, setAiIncludeReferences] = useState(true);
  const [aiGenerating, setAiGenerating] = useState(false);
  const [aiGeneratedMolecules, setAiGeneratedMolecules] = useState([]);
  const [aiRequestId, setAiRequestId] = useState(null);
  const [rawApiResponse, setRawApiResponse] = useState(null); // New state for raw response
  
  const [editorSmiles, setEditorSmiles] = useState('');
  const [chatContext, setChatContext] = useState('General molecule design'); // Context for the chat
  
  // Load saved molecules from localStorage on component mount
  useEffect(() => {
    const loadSavedMolecules = () => {
      try {
        const saved = localStorage.getItem('savedMolecules');
        if (saved) {
          setSavedMolecules(JSON.parse(saved));
        }
      } catch (err) {
        console.error('Error loading saved molecules:', err);
        showSnackbar('Failed to load saved molecules', 'error');
      }
    };
    
    loadSavedMolecules();
  }, []);
  
  // Save to localStorage when savedMolecules changes
  useEffect(() => {
    if (savedMolecules.length > 0) {
      try {
        localStorage.setItem('savedMolecules', JSON.stringify(savedMolecules));
      } catch (err) {
        console.error('Error saving molecules:', err);
      }
    }
  }, [savedMolecules]);

  const handleTabChange = (event, newValue) => {
    setActiveTab(newValue);
  };
  
  const showSnackbar = (message, severity = 'success') => {
    setSnackbarMessage(message);
    setSnackbarSeverity(severity);
    setSnackbarOpen(true);
  };
  
  const handleSnackbarClose = (event, reason) => {
    if (reason === 'clickaway') {
      return;
    }
    setSnackbarOpen(false);
  };

  const handleGenerate = () => {
    setLoading(true);
    setError(null);
    
    // Simulate API call
    setTimeout(() => {
      setLoading(false);
      
      const newMolecules = [
        {
          id: Date.now(),
          name: `Compound-${Math.floor(Math.random() * 1000)}`,
          smiles: 'CC1=C(C(=O)OC1=O)C2=CC=CC=C2',
          properties: {
            molecularWeight: '228.24 g/mol',
            logP: '1.8',
            hBondDonors: '0',
            hBondAcceptors: '3',
            rotableBonds: '2',
            psa: '52.6 Å²',
            targetAffinity: targetReceptor ? `${(Math.random() * 10 - 12).toFixed(1)} kcal/mol` : 'N/A',
          },
          timestamp: new Date().toISOString(),
        },
        {
          id: Date.now() + 1,
          name: `Compound-${Math.floor(Math.random() * 1000)}`,
          smiles: 'CC1=CC=C(C=C1)C2=CC(=NN2C)C(=O)N',
          properties: {
            molecularWeight: '215.25 g/mol',
            logP: '2.1',
            hBondDonors: '1',
            hBondAcceptors: '4',
            rotableBonds: '2',
            psa: '71.2 Å²',
            targetAffinity: targetReceptor ? `${(Math.random() * 10 - 12).toFixed(1)} kcal/mol` : 'N/A',
          },
          timestamp: new Date().toISOString(),
        }
      ];
      
      setGeneratedMolecules(newMolecules);
      setSelectedMolecule(newMolecules[0]);
    }, 2000);
  };
  
  // Generate molecules using Claude AI
  const handleAiGenerate = async () => {
    try {
      if (!aiRequirements.trim()) {
        showSnackbar('Please enter generation requirements', 'error');
        return;
      }
      
      setAiGenerating(true);
      setError(null);
      setAiGeneratedMolecules([]); // Ensure state is cleared BEFORE the API call
      setSelectedMolecule(null); // Clear selected molecule too
      setAiRequestId(null);
      setRawApiResponse(null); // Reset raw response
      
      const response = await claudeAPI.generateMolecule({
        requirements: aiRequirements,
        targetReceptors: aiTargetReceptors,
        minimizeSideEffects: aiMinimizeSideEffects,
        includeLiteratureReferences: aiIncludeReferences
      });
      
      // Process response - Backend should provide validated molecules
      const backendMolecules = response.data.molecules || [];
      
      // Map backend molecules to frontend state, ensuring unique IDs if possible
      // Let's assume backend `id` or `smiles` can be used for uniqueness
      const molecules = backendMolecules.map((molData, index) => ({
          // Use a more robust ID if backend provides one, otherwise fallback
          id: molData.id || `${Date.now()}-${index}`,
          name: molData.name || `AI-Compound-${index+1}`,
          smiles: molData.smiles,
          properties: molData.properties || {},
          admet: molData.admet || null,
          timestamp: molData.dateCreated || new Date().toISOString(), // Use backend date if available
          aiGenerated: true,
      }));

      if (molecules.length === 0) {
        // Handle case where backend returned no valid molecules
        throw new Error('No valid molecules were generated or extracted by the AI.');
      }
      
      // Set generated molecules
      setAiGeneratedMolecules(molecules);
      setSelectedMolecule(molecules[0]); // Select the first valid one
      
      // Save request ID and raw response text
      setAiRequestId(response.data.requestId);
      setRawApiResponse(response.data.rawClaudeResponse || 'No text response received.');
      showSnackbar(`AI generated ${molecules.length} molecules successfully`, 'success');
      
      console.log("[MoleculeDesigner] AI Generation Success:", { molecules, requestId: response.data.requestId, rawResponseLength: response.data.rawClaudeResponse?.length });
      
    } catch (err) {
      console.error('Error generating molecules with AI:', err);
      const errorMsg = err.response?.data?.error || err.message || 'Failed to generate molecules with AI';
      setError(errorMsg);
      showSnackbar(errorMsg, 'error');
      setRawApiResponse(null); // Clear response on error
      setAiGeneratedMolecules([]); // Clear molecules on error
      setSelectedMolecule(null);
    } finally {
      setAiGenerating(false);
    }
  };
  
  const handleSaveMolecule = (molecule) => {
    // Check if molecule already exists in savedMolecules
    const exists = savedMolecules.some(m => m.id === molecule.id);
    
    if (!exists) {
      const updatedSavedMolecules = [...savedMolecules, molecule];
      setSavedMolecules(updatedSavedMolecules);
      showSnackbar(`${molecule.name} saved successfully`);
    } else {
      showSnackbar(`${molecule.name} is already saved`, 'warning');
    }
  };
  
  const handleDeleteMolecule = (id) => {
    const updatedMolecules = savedMolecules.filter(molecule => molecule.id !== id);
    setSavedMolecules(updatedMolecules);
    showSnackbar('Molecule deleted successfully');
  };
  
  const handleMoleculeSelection = (molecule) => {
    setSelectedMolecule(molecule);
  };
  
  const handleEditDialogOpen = (molecule) => {
    setEditMolecule(molecule);
    setMoleculeName(molecule.name);
    setSmilesStructure(molecule.smiles);
    setEditDialogOpen(true);
  };
  
  const handleEditDialogClose = () => {
    setEditDialogOpen(false);
    setEditMolecule(null);
    setMoleculeName('');
    setSmilesStructure('');
  };
  
  const handleEditSave = () => {
    if (!moleculeName.trim() || !smilesStructure.trim()) {
      showSnackbar('Please fill in all fields', 'error');
      return;
    }
    
    const updatedMolecules = savedMolecules.map(molecule => 
      molecule.id === editMolecule.id 
        ? { ...molecule, name: moleculeName, smiles: smilesStructure }
        : molecule
    );
    
    setSavedMolecules(updatedMolecules);
    showSnackbar('Molecule updated successfully');
    handleEditDialogClose();
  };
  
  const handleSimulationDialogOpen = (molecule) => {
    setSelectedMolecule(molecule);
    setSimulationDialogOpen(true);
  };
  
  const handleSimulationDialogClose = () => {
    setSimulationDialogOpen(false);
  };
  
  const handleRunSimulation = () => {
    // In a real app, this would trigger a simulation in the Simulations component
    showSnackbar(`Simulation request for ${selectedMolecule.name} sent to Simulations`);
    handleSimulationDialogClose();
    
    // This is where integration with the Simulations component would happen
    // You would dispatch an action or use a context to send the molecule data
    console.log('Running simulation with:', {
      molecule: selectedMolecule,
      simulationType: simulationType
    });
  };
  
  const handleExportMolecules = () => {
    try {
      const dataStr = JSON.stringify(savedMolecules, null, 2);
      const dataUri = 'data:application/json;charset=utf-8,'+ encodeURIComponent(dataStr);
      
      const exportFileDefaultName = 'molecules.json';
      
      const linkElement = document.createElement('a');
      linkElement.setAttribute('href', dataUri);
      linkElement.setAttribute('download', exportFileDefaultName);
      linkElement.click();
      
      showSnackbar('Molecules exported successfully');
    } catch (err) {
      console.error('Error exporting molecules:', err);
      showSnackbar('Failed to export molecules', 'error');
    }
  };
  
  const handleImportMolecules = (event) => {
    const file = event.target.files[0];
    if (!file) return;
    
    const reader = new FileReader();
    reader.onload = (e) => {
      try {
        const importedMolecules = JSON.parse(e.target.result);
        if (Array.isArray(importedMolecules)) {
          // Merge with existing molecules, avoiding duplicates by ID
          const existingIds = savedMolecules.map(m => m.id);
          const newMolecules = importedMolecules.filter(m => !existingIds.includes(m.id));
          
          if (newMolecules.length > 0) {
            setSavedMolecules([...savedMolecules, ...newMolecules]);
            showSnackbar(`Imported ${newMolecules.length} new molecules successfully`);
          } else {
            showSnackbar('No new molecules to import', 'info');
          }
        } else {
          showSnackbar('Invalid import file format', 'error');
        }
      } catch (err) {
        console.error('Error importing molecules:', err);
        showSnackbar('Failed to import molecules: Invalid JSON', 'error');
      }
    };
    reader.readAsText(file);
    
    // Reset the file input
    event.target.value = '';
  };

  const receptors = [
    { value: 'D1', label: 'Dopamine D1 Receptor' },
    { value: 'D2', label: 'Dopamine D2 Receptor' },
    { value: 'DAT', label: 'Dopamine Transporter (DAT)' },
    { value: 'NET', label: 'Norepinephrine Transporter (NET)' },
    { value: '5HT2A', label: 'Serotonin 5-HT2A Receptor' },
    { value: 'AMPA', label: 'AMPA Receptor' },
  ];
  
  const simulationTypes = [
    { value: 'docking', label: 'Molecular Docking' },
    { value: 'md', label: 'Molecular Dynamics' },
    { value: 'admet', label: 'ADMET Prediction' },
  ];

  const formatDate = (dateString) => {
    const options = { year: 'numeric', month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' };
    return new Date(dateString).toLocaleDateString(undefined, options);
  };

  // Function to handle changes from StructureEditor
  const handleEditorMoleculeChange = (newSmiles) => {
    setEditorSmiles(newSmiles);
    // Optionally, select this molecule for viewing/saving
    // This requires converting SMILES back to a molecule object structure
    // For now, just update the SMILES state
  };

  // Handler for when the chat mentions a molecule SMILES
  const handleMoleculeMentionedInChat = (smiles) => {
    console.log('Chat mentioned SMILES:', smiles);
    // Potentially load this SMILES into the viewer or editor
    // Example: If a molecule object can be created/found, update selectedMolecule
    // For now, just log it
    showSnackbar(`AI mentioned molecule: ${smiles}`, 'info');
  };

  // Update chat context when a molecule is selected
  useEffect(() => {
    if (selectedMolecule) {
      setChatContext(`Discussing molecule: ${selectedMolecule.name} (SMILES: ${selectedMolecule.smiles})`);
    } else {
      setChatContext('General molecule design');
    }
  }, [selectedMolecule]);

  return (
    <div className={classes.root}>
      <Typography variant="h4" className={classes.title}>
        Molecule Designer
      </Typography>
      
      <Paper className={classes.paper}>
        <Tabs
          value={activeTab}
          onChange={handleTabChange}
          indicatorColor="primary"
          textColor="primary"
        >
          <Tab label="Generator" />
          <Tab label="AI Designer" />
          <Tab label="Structure Editor" />
          <Tab label="Saved Molecules" />
        </Tabs>
        
        <TabPanel value={activeTab} index={0}>
          <Grid container spacing={3}>
            <Grid item xs={12} md={6}>
              <Typography variant="h6" gutterBottom>
                Molecule Generation Parameters
              </Typography>
              
              <FormControl variant="outlined" className={classes.formControl}>
                <InputLabel>Target Receptor</InputLabel>
                <Select
                  value={targetReceptor}
                  onChange={(e) => setTargetReceptor(e.target.value)}
                  label="Target Receptor"
                >
                  <MenuItem value=""><em>None</em></MenuItem>
                  {receptors.map((receptor) => (
                    <MenuItem key={receptor.value} value={receptor.value}>
                      {receptor.label}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
              
              <TextField
                className={classes.formControl}
                label="Starting Structure (SMILES or Scaffold)"
                variant="outlined"
                multiline
                rows={4}
                value={structureInput}
                onChange={(e) => setStructureInput(e.target.value)}
                placeholder="Enter SMILES notation or describe desired scaffold..."
              />
              
              <Button
                variant="contained"
                color="primary"
                className={classes.button}
                onClick={handleGenerate}
                disabled={loading}
                fullWidth
                startIcon={<BuildIcon />}
              >
                {loading ? <CircularProgress size={24} /> : 'Generate Molecules'}
              </Button>
              
              {error && (
                <Alert severity="error" style={{ marginTop: 16 }}>
                  {error}
                </Alert>
              )}
            </Grid>
            
            <Grid item xs={12} md={6}>
              <Typography variant="h6" gutterBottom>
                Molecule Viewer
              </Typography>
              
              {loading ? (
                <div className={classes.loadingContainer}>
                  <CircularProgress />
                  <Typography variant="body1" style={{ marginTop: 16 }}>
                    Generating molecules...
                  </Typography>
                </div>
              ) : (
                <div>
                  {selectedMolecule ? (
                    <div>
                      <MoleculeViewer3D smiles={selectedMolecule.smiles} />
                      
                      <Card className={classes.resultCard}>
                        <div className={classes.actionButtons}>
                          <IconButton 
                            size="small" 
                            onClick={() => handleSaveMolecule(selectedMolecule)}
                            title="Save molecule"
                          >
                            <SaveIcon fontSize="small" />
                          </IconButton>
                          <IconButton 
                            size="small" 
                            onClick={() => handleSimulationDialogOpen(selectedMolecule)}
                            title="Run simulation"
                          >
                            <CompareArrowsIcon fontSize="small" />
                          </IconButton>
                        </div>
                        <CardContent>
                          <Typography variant="h6">{selectedMolecule.name}</Typography>
                          <Typography variant="body2" color="textSecondary" gutterBottom>
                            SMILES: {selectedMolecule.smiles}
                          </Typography>
                          
                          <Typography variant="subtitle1" gutterBottom style={{ marginTop: 8 }}>
                            Properties:
                          </Typography>
                          
                          <Grid container spacing={2} className={classes.propertyGrid}>
                            {selectedMolecule.properties && typeof selectedMolecule.properties === 'object' && Object.entries(selectedMolecule.properties).map(([key, value]) => (
                              <Grid item xs={6} sm={4} key={key}>
                                <Typography variant="body2" color="textSecondary">
                                  {key}:
                                </Typography>
                                <Typography variant="body1">
                                  {value === null ? 'null' : typeof value === 'object' ? '[Object]' : String(value)}
                                </Typography>
                              </Grid>
                            ))}
                          </Grid>
                        </CardContent>
                      </Card>
                      
                      {generatedMolecules.length > 1 && (
                        <div style={{ marginTop: 16 }}>
                          <Typography variant="subtitle1" gutterBottom>
                            Other Generated Molecules:
                          </Typography>
                          <div style={{ display: 'flex', flexWrap: 'wrap' }}>
                            {generatedMolecules
                              .filter(m => m.id !== selectedMolecule.id)
                              .map(molecule => (
                                <Chip
                                  key={molecule.id}
                                  label={molecule.name}
                                  clickable
                                  onClick={() => handleMoleculeSelection(molecule)}
                                  className={classes.chip}
                                />
                              ))}
                          </div>
                        </div>
                      )}
                    </div>
                  ) : generatedMolecules.length > 0 ? (
                    <div>
                      <div className={classes.canvasContainer}>
                        <Typography variant="body1" color="textSecondary">
                          Select a molecule to view details
                        </Typography>
                      </div>
                      <div style={{ marginTop: 16 }}>
                        <Typography variant="subtitle1" gutterBottom>
                          Generated Molecules:
                        </Typography>
                        <div style={{ display: 'flex', flexWrap: 'wrap' }}>
                          {generatedMolecules.map(molecule => (
                            <Chip
                              key={molecule.id}
                              label={molecule.name}
                              clickable
                              onClick={() => handleMoleculeSelection(molecule)}
                              className={classes.chip}
                            />
                          ))}
                        </div>
                      </div>
                    </div>
                  ) : (
                    <div className={classes.canvasContainer}>
                      <Typography variant="body1" color="textSecondary">
                        No molecules generated yet. Set parameters and click "Generate Molecules".
                      </Typography>
                    </div>
                  )}
                </div>
              )}
            </Grid>
          </Grid>
        </TabPanel>
        
        <TabPanel value={activeTab} index={1}>
          <Grid container spacing={3}>
            <Grid item xs={12} md={5}>
              <Typography variant="h6" gutterBottom>
                AI-Powered Molecule Design
              </Typography>
              
              <Card className={classes.aiParameterCard}>
                <CardContent>
                  <Typography variant="subtitle1" gutterBottom>
                    Design Requirements
                  </Typography>
                  
                  <TextField
                    className={classes.generationRequirements}
                    label="Design Goals and Requirements"
                    variant="outlined"
                    multiline
                    rows={4}
                    fullWidth
                    value={aiRequirements}
                    onChange={(e) => setAiRequirements(e.target.value)}
                    placeholder="Describe what you want in the molecule..."
                  />
                  
                  <Typography variant="subtitle1" gutterBottom>
                    Target Parameters
                  </Typography>
                  
                  <TextField
                    className={classes.formControl}
                    label="Target Receptors"
                    variant="outlined"
                    fullWidth
                    value={aiTargetReceptors}
                    onChange={(e) => setAiTargetReceptors(e.target.value)}
                    placeholder="e.g., dopamine/norepinephrine transporters"
                  />
                  
                  <FormControlLabel
                    control={
                      <Checkbox
                        checked={aiMinimizeSideEffects}
                        onChange={(e) => setAiMinimizeSideEffects(e.target.checked)}
                        color="primary"
                      />
                    }
                    label="Minimize side effects"
                  />
                  
                  <FormControlLabel
                    control={
                      <Checkbox
                        checked={aiIncludeReferences}
                        onChange={(e) => setAiIncludeReferences(e.target.checked)}
                        color="primary"
                      />
                    }
                    label="Include literature references"
                  />
                  
                  <Button
                    variant="contained"
                    color="primary"
                    fullWidth
                    onClick={handleAiGenerate}
                    disabled={aiGenerating || !aiRequirements.trim()}
                    startIcon={<BuildIcon />}
                    className={classes.button}
                  >
                    {aiGenerating ? <CircularProgress size={24} /> : 'Generate with AI'}
                  </Button>
                </CardContent>
              </Card>
              
              {aiGenerating && (
                <Paper className={classes.paper} style={{ marginTop: 16 }}>
                  <Typography variant="body1" align="center" gutterBottom>
                    AI is designing molecules based on your requirements...
                  </Typography>
                  <LinearProgress />
                  <Typography variant="body2" align="center" style={{ marginTop: 8 }}>
                    This may take a minute or two...
                  </Typography>
                </Paper>
              )}
            </Grid>
            
            <Grid item xs={12} md={7}>
              <Typography variant="h6" gutterBottom>
                Generated Molecules
              </Typography>
              {aiGeneratedMolecules.length > 0 ? (
                <div>
                  {selectedMolecule && (
                    <div>
                      <MoleculeViewer3D smiles={selectedMolecule.smiles} />
                      <Card className={classes.resultCard}>
                        <div className={classes.actionButtons}>
                          <IconButton 
                            size="small" 
                            onClick={() => handleSaveMolecule(selectedMolecule)}
                            title="Save molecule"
                          >
                            <SaveIcon fontSize="small" />
                          </IconButton>
                          <IconButton 
                            size="small" 
                            onClick={() => handleSimulationDialogOpen(selectedMolecule)}
                            title="Run simulation"
                          >
                            <CompareArrowsIcon fontSize="small" />
                          </IconButton>
                        </div>
                        <CardContent>
                          <Typography variant="h6">{selectedMolecule.name}</Typography>
                          <Typography variant="body2" color="textSecondary" gutterBottom>
                            SMILES: {selectedMolecule.smiles}
                          </Typography>
                          <Typography variant="subtitle1" gutterBottom style={{ marginTop: 8 }}>
                            Properties:
                          </Typography>
                          <Grid container spacing={2} className={classes.propertyGrid}>
                            {selectedMolecule.properties && typeof selectedMolecule.properties === 'object' && Object.entries(selectedMolecule.properties).map(([key, value]) => (
                              <Grid item xs={6} sm={4} key={key}>
                                <Typography variant="body2" color="textSecondary">{key}:</Typography>
                                <Typography variant="body1">{value === null ? 'null' : typeof value === 'object' ? '[Object]' : String(value)}</Typography>
                              </Grid>
                            ))}
                          </Grid>
                        </CardContent>
                      </Card>
                      {aiGeneratedMolecules.filter(m => m.id !== selectedMolecule.id).length > 0 && (
                        <div style={{ marginTop: 16 }}>
                          <Typography variant="subtitle1" gutterBottom>
                            Other AI-Generated Molecules:
                          </Typography>
                          <div style={{ display: 'flex', flexWrap: 'wrap' }}>
                            {[...new Map(aiGeneratedMolecules.filter(m => m.id !== selectedMolecule.id).map(item => [item.smiles, item])).values()]
                              .map(molecule => (
                                <Chip
                                  key={molecule.smiles}
                                  label={molecule.name}
                                  clickable
                                  onClick={() => handleMoleculeSelection(molecule)}
                                  className={classes.chip}
                                  color="primary"
                                  variant="outlined"
                                />
                              ))}
                          </div>
                        </div>
                      )}
                    </div>
                  )}
                  {!selectedMolecule && aiGeneratedMolecules.length > 0 && (
                    <div className={classes.canvasContainer}>
                      <Typography variant="body1" color="textSecondary" align="center">
                        Select a generated molecule to view details.
                      </Typography>
                    </div>
                  )}
                </div>
              ) : (
                <div className={classes.canvasContainer}>
                  <Typography variant="body1" color="textSecondary" align="center">
                    {aiGenerating ? (
                      'Working on your molecular designs...'
                    ) : (
                      'No AI-generated molecules yet. Set your requirements and click "Generate with AI".'
                    )}
                  </Typography>
                </div>
              )}

              {/* API Response Section - Full Width Below Columns */}
              {rawApiResponse && (
                <Grid item xs={12} className={classes.apiResponseGridItem}>
                  <Divider className={classes.apiResponseDivider}/> 
                  <Typography variant="h6" gutterBottom>
                    Claude API Response Text
                  </Typography>
                  {/* Render markdown directly without Paper/Card wrapper */}
                  <ReactMarkdown className={classes.apiResponseContent}>
                    {rawApiResponse}
                  </ReactMarkdown>
                </Grid>
              )}

              {/* --- CHATBOT MOVED HERE --- */}
              <Grid item xs={12} className={classes.stickyContainer}>
                 <Divider style={{marginBottom: '16px'}}/>
                 <Typography variant="h6" gutterBottom>
                   Chat with AI Assistant
                 </Typography>
                 <AIChatInterfaceRedesigned 
                   initialContext={chatContext} 
                   onMoleculeMentioned={handleMoleculeMentionedInChat}
                   selectedMolecule={selectedMolecule} 
                 /> 
              </Grid>
            </Grid>
          </Grid>
        </TabPanel>
        
        <TabPanel value={activeTab} index={2}>
          <StructureEditor 
            initialMolecule={editorSmiles} // Pass current editor SMILES or initial value
            onMoleculeChange={handleEditorMoleculeChange} // Handle updates from the editor
            showControls={true}
            showProperties={true}
          />
        </TabPanel>
        
        <TabPanel value={activeTab} index={3}>
          <Grid container spacing={3}>
            <Grid item xs={12}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 16 }}>
                <Typography variant="h6">
                  Saved Molecules
                </Typography>
                
                <div>
                  <input
                    accept=".json"
                    className={classes.uploadInput}
                    id="import-molecules-button"
                    type="file"
                    onChange={handleImportMolecules}
                  />
                  <label htmlFor="import-molecules-button">
                    <Button
                      variant="outlined"
                      component="span"
                      startIcon={<CloudUploadIcon />}
                      style={{ marginRight: 8 }}
                    >
                      Import
                    </Button>
                  </label>
                  
                  <Button 
                    variant="outlined" 
                    startIcon={<CloudDownloadIcon />}
                    onClick={handleExportMolecules}
                    disabled={savedMolecules.length === 0}
                  >
                    Export
                  </Button>
                </div>
              </div>
              
              {savedMolecules.length === 0 ? (
                <Paper style={{ padding: 24, textAlign: 'center' }}>
                  <Typography variant="body1" color="textSecondary">
                    No molecules saved yet. Generate and save molecules to see them here.
                  </Typography>
                </Paper>
              ) : (
                <TableContainer component={Paper} className={classes.tableContainer}>
                  <Table stickyHeader aria-label="saved molecules table">
                    <TableHead>
                      <TableRow>
                        <TableCell>Name</TableCell>
                        <TableCell>SMILES</TableCell>
                        <TableCell>Molecular Weight</TableCell>
                        <TableCell>Date Added</TableCell>
                        <TableCell align="center">Actions</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {savedMolecules.map((molecule) => (
                        <TableRow 
                          key={molecule.id} 
                          className={classes.tableRow}
                          hover
                          onClick={() => handleMoleculeSelection(molecule)}
                        >
                          <TableCell>{molecule.name}</TableCell>
                          <TableCell>
                            {molecule.smiles.length > 20 
                              ? `${molecule.smiles.substring(0, 20)}...` 
                              : molecule.smiles}
                          </TableCell>
                          <TableCell>{molecule.properties?.molecularWeight || 'N/A'}</TableCell>
                          <TableCell>{molecule.timestamp ? formatDate(molecule.timestamp) : 'N/A'}</TableCell>
                          <TableCell align="center">
                            <IconButton 
                              size="small" 
                              onClick={(e) => {
                                e.stopPropagation();
                                handleEditDialogOpen(molecule);
                              }}
                              title="Edit molecule"
                            >
                              <EditIcon fontSize="small" />
                            </IconButton>
                            <IconButton 
                              size="small" 
                              onClick={(e) => {
                                e.stopPropagation();
                                handleSimulationDialogOpen(molecule);
                              }}
                              title="Run simulation"
                            >
                              <CompareArrowsIcon fontSize="small" />
                            </IconButton>
                            <IconButton 
                              size="small" 
                              onClick={(e) => {
                                e.stopPropagation();
                                handleDeleteMolecule(molecule.id);
                              }}
                              title="Delete molecule"
                            >
                              <DeleteIcon fontSize="small" />
                            </IconButton>
                          </TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </TableContainer>
              )}
            </Grid>
            
            {selectedMolecule && (
              <Grid item xs={12} md={6}>
                <Typography variant="h6" gutterBottom>
                  Selected Molecule Preview
                </Typography>
                <MoleculeViewer3D smiles={selectedMolecule.smiles} />
              </Grid>
            )}
          </Grid>
        </TabPanel>
      </Paper>
      
      {/* Edit Molecule Dialog */}
      <Dialog open={editDialogOpen} onClose={handleEditDialogClose} maxWidth="sm" fullWidth>
        <DialogTitle>Edit Molecule</DialogTitle>
        <DialogContent>
          <TextField
            autoFocus
            margin="dense"
            label="Molecule Name"
            type="text"
            fullWidth
            value={moleculeName}
            onChange={(e) => setMoleculeName(e.target.value)}
            variant="outlined"
            className={classes.inputField}
          />
          <TextField
            margin="dense"
            label="SMILES Structure"
            type="text"
            fullWidth
            value={smilesStructure}
            onChange={(e) => setSmilesStructure(e.target.value)}
            variant="outlined"
            className={classes.inputField}
            multiline
            rows={4}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={handleEditDialogClose} color="primary">
            Cancel
          </Button>
          <Button onClick={handleEditSave} color="primary" variant="contained">
            Save Changes
          </Button>
        </DialogActions>
      </Dialog>
      
      {/* Simulation Dialog */}
      <Dialog open={simulationDialogOpen} onClose={handleSimulationDialogClose} maxWidth="sm" fullWidth>
        <DialogTitle>Run Simulation</DialogTitle>
        <DialogContent>
          {selectedMolecule && (
            <>
              <Typography variant="subtitle1" gutterBottom>
                {selectedMolecule.name}
              </Typography>
              <Typography variant="body2" color="textSecondary" gutterBottom>
                SMILES: {selectedMolecule.smiles}
              </Typography>
              
              <FormControl variant="outlined" className={classes.formControl} style={{ marginTop: 16 }}>
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
              
              <FormControl variant="outlined" className={classes.formControl}>
                <InputLabel>Target Receptor</InputLabel>
                <Select
                  value={targetReceptor}
                  onChange={(e) => setTargetReceptor(e.target.value)}
                  label="Target Receptor"
                >
                  <MenuItem value=""><em>None</em></MenuItem>
                  {receptors.map((receptor) => (
                    <MenuItem key={receptor.value} value={receptor.value}>
                      {receptor.label}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            </>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={handleSimulationDialogClose} color="primary">
            Cancel
          </Button>
          <Button 
            onClick={handleRunSimulation} 
            color="primary" 
            variant="contained"
            startIcon={<CompareArrowsIcon />}
          >
            Run Simulation
          </Button>
        </DialogActions>
      </Dialog>
      
      {/* Snackbar for notifications */}
      <Snackbar 
        open={snackbarOpen} 
        autoHideDuration={6000} 
        onClose={handleSnackbarClose}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
      >
        <Alert onClose={handleSnackbarClose} severity={snackbarSeverity}>
          {snackbarMessage}
        </Alert>
      </Snackbar>
    </div>
  );
};

export default MoleculeDesigner; 