import React, { useState, useEffect } from 'react';
import { 
  Typography, 
  Paper, 
  Grid, 
  Button, 
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  List,
  ListItem,
  ListItemText,
  ListItemSecondaryAction,
  IconButton,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Divider,
  makeStyles,
  CircularProgress,
  Chip,
  Card,
  CardContent,
  Snackbar
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import GetAppIcon from '@material-ui/icons/GetApp';
import PublishIcon from '@material-ui/icons/Publish';
import DeleteIcon from '@material-ui/icons/Delete';
import ShareIcon from '@material-ui/icons/Share';
import EditIcon from '@material-ui/icons/Edit';
import SaveIcon from '@material-ui/icons/Save';
import FileCopyIcon from '@material-ui/icons/FileCopy';

// This would be imported from a third-party library in a real implementation
const MoleculeViewer = ({ smiles, height = 250 }) => {
  if (!smiles) return <div style={{ height: height, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>No molecule to display</div>;
  
  return (
    <div style={{ textAlign: 'center', height: height, display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center', backgroundColor: '#f5f5f5', padding: '10px' }}>
      <div style={{ marginBottom: '10px' }}>
        <strong>SMILES:</strong> {smiles}
      </div>
      <div>
        [Molecule Visualization Placeholder]
      </div>
    </div>
  );
};

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    padding: theme.spacing(3),
  },
  paper: {
    padding: theme.spacing(3),
    height: '100%',
  },
  title: {
    marginBottom: theme.spacing(3),
  },
  formControl: {
    marginBottom: theme.spacing(2),
    width: '100%',
  },
  button: {
    marginTop: theme.spacing(1),
    marginRight: theme.spacing(1),
  },
  chip: {
    margin: theme.spacing(0.5),
  },
  divider: {
    margin: theme.spacing(2, 0),
  },
  list: {
    maxHeight: 500,
    overflow: 'auto',
  },
  listItem: {
    borderLeft: `3px solid ${theme.palette.primary.main}`,
    marginBottom: theme.spacing(1),
    '&:hover': {
      backgroundColor: '#f5f5f5',
    },
  },
  selectedListItem: {
    borderLeft: `3px solid ${theme.palette.secondary.main}`,
    backgroundColor: '#e3f2fd',
    marginBottom: theme.spacing(1),
  },
  exportImportButtons: {
    display: 'flex',
    justifyContent: 'space-between',
    marginTop: theme.spacing(2),
  },
  fileInput: {
    display: 'none',
  },
  card: {
    marginBottom: theme.spacing(2),
  },
  moleculeDisplay: {
    marginTop: theme.spacing(2),
  },
  noResults: {
    padding: theme.spacing(2),
    textAlign: 'center',
    color: theme.palette.text.secondary,
  },
}));

const MoleculeStorage = ({ onSelectMolecule, currentMolecule }) => {
  const classes = useStyles();
  const [molecules, setMolecules] = useState([]);
  const [filteredMolecules, setFilteredMolecules] = useState([]);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [searchTerm, setSearchTerm] = useState('');
  const [loading, setLoading] = useState(true);
  const [categoryFilter, setCategoryFilter] = useState('all');
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [moleculeToDelete, setMoleculeToDelete] = useState(null);
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [snackbarMessage, setSnackbarMessage] = useState('');
  const [snackbarSeverity, setSnackbarSeverity] = useState('success');
  const [exportFormat, setExportFormat] = useState('json');
  const fileInputRef = React.useRef(null);
  
  // Simulated categories for molecules
  const categories = [
    'all', 
    'stimulants', 
    'adhd treatments', 
    'novel compounds', 
    'derivatives', 
    'under review'
  ];
  
  // Load molecules from localStorage on component mount
  useEffect(() => {
    const loadMolecules = () => {
      try {
        // Simulate a small delay to show loading state
        setTimeout(() => {
          const storedMolecules = localStorage.getItem('molecules');
          
          if (storedMolecules) {
            const parsedMolecules = JSON.parse(storedMolecules);
            setMolecules(parsedMolecules);
            setFilteredMolecules(parsedMolecules);
          } else {
            // If no molecules in storage, create some example molecules
            const exampleMolecules = [
              {
                id: '1',
                name: 'Methylphenidate',
                smiles: 'CN(C)C(C1=CC=CC=C1)C(C)OC(=O)C',
                dateCreated: '2023-01-15',
                categories: ['stimulants', 'adhd treatments'],
                properties: {
                  molecularWeight: '233.31 g/mol',
                  logP: '2.15',
                }
              },
              {
                id: '2',
                name: 'Amphetamine',
                smiles: 'CC(N)CC1=CC=CC=C1',
                dateCreated: '2023-02-20',
                categories: ['stimulants', 'adhd treatments'],
                properties: {
                  molecularWeight: '135.21 g/mol',
                  logP: '1.76',
                }
              },
              {
                id: '3',
                name: 'Novel Dopamine Modulator X-42',
                smiles: 'CC1=CC(=C(C=C1)NC(=O)C)NC2=CC=CC=C2',
                dateCreated: '2023-05-03',
                categories: ['novel compounds', 'under review'],
                properties: {
                  molecularWeight: '254.32 g/mol',
                  logP: '3.21',
                }
              }
            ];
            
            setMolecules(exampleMolecules);
            setFilteredMolecules(exampleMolecules);
            // Save to localStorage
            localStorage.setItem('molecules', JSON.stringify(exampleMolecules));
          }
          
          setLoading(false);
        }, 1000);
      } catch (error) {
        console.error('Error loading molecules:', error);
        setLoading(false);
        showSnackbar('Error loading molecules from storage', 'error');
      }
    };
    
    loadMolecules();
  }, []);
  
  // Add current molecule if it exists and is provided
  useEffect(() => {
    if (currentMolecule && currentMolecule.name && currentMolecule.smiles) {
      saveMolecule(currentMolecule);
    }
  }, [currentMolecule]);
  
  // Update filtered molecules when search term or category changes
  useEffect(() => {
    filterMolecules();
  }, [searchTerm, categoryFilter, molecules]);
  
  const filterMolecules = () => {
    let filtered = [...molecules];
    
    // Filter by search term
    if (searchTerm) {
      filtered = filtered.filter(molecule => 
        molecule.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
        molecule.smiles.toLowerCase().includes(searchTerm.toLowerCase())
      );
    }
    
    // Filter by category
    if (categoryFilter && categoryFilter !== 'all') {
      filtered = filtered.filter(molecule => 
        molecule.categories && molecule.categories.includes(categoryFilter)
      );
    }
    
    setFilteredMolecules(filtered);
  };
  
  const handleSearchChange = (event) => {
    setSearchTerm(event.target.value);
  };
  
  const handleCategoryChange = (event) => {
    setCategoryFilter(event.target.value);
  };
  
  const handleSelectMolecule = (molecule) => {
    setSelectedMolecule(molecule);
    
    // Notify parent component if callback provided
    if (onSelectMolecule) {
      onSelectMolecule(molecule);
    }
  };
  
  const handleDeleteClick = (molecule) => {
    setMoleculeToDelete(molecule);
    setDeleteDialogOpen(true);
  };
  
  const confirmDelete = () => {
    if (moleculeToDelete) {
      const updatedMolecules = molecules.filter(m => m.id !== moleculeToDelete.id);
      setMolecules(updatedMolecules);
      
      // Update localStorage
      localStorage.setItem('molecules', JSON.stringify(updatedMolecules));
      
      // Clear selection if deleted molecule was selected
      if (selectedMolecule && selectedMolecule.id === moleculeToDelete.id) {
        setSelectedMolecule(null);
      }
      
      showSnackbar('Molecule deleted successfully', 'success');
    }
    
    setDeleteDialogOpen(false);
    setMoleculeToDelete(null);
  };
  
  const cancelDelete = () => {
    setDeleteDialogOpen(false);
    setMoleculeToDelete(null);
  };
  
  const saveMolecule = (molecule) => {
    try {
      // Check if this is an update or new molecule
      const existingIndex = molecules.findIndex(m => m.id === molecule.id);
      let updatedMolecules;
      
      if (existingIndex >= 0) {
        // Update existing molecule
        updatedMolecules = [...molecules];
        updatedMolecules[existingIndex] = molecule;
      } else {
        // Add new molecule with generated ID if it doesn't have one
        const newMolecule = {
          ...molecule,
          id: molecule.id || `molecule-${Date.now()}`,
          dateCreated: molecule.dateCreated || new Date().toISOString().split('T')[0]
        };
        updatedMolecules = [...molecules, newMolecule];
      }
      
      setMolecules(updatedMolecules);
      
      // Update localStorage
      localStorage.setItem('molecules', JSON.stringify(updatedMolecules));
      
      showSnackbar('Molecule saved successfully', 'success');
      
      return true;
    } catch (error) {
      console.error('Error saving molecule:', error);
      showSnackbar('Error saving molecule', 'error');
      return false;
    }
  };
  
  const exportMolecules = () => {
    try {
      let dataToExport;
      let filename;
      let mimeType;
      
      if (exportFormat === 'json') {
        dataToExport = JSON.stringify(selectedMolecule ? [selectedMolecule] : molecules, null, 2);
        filename = selectedMolecule ? `${selectedMolecule.name}.json` : 'molecules.json';
        mimeType = 'application/json';
      } else if (exportFormat === 'csv') {
        // Convert to CSV
        const dataArray = selectedMolecule ? [selectedMolecule] : molecules;
        const headers = ['id', 'name', 'smiles', 'dateCreated', 'categories'];
        const csvRows = [headers.join(',')];
        
        for (const molecule of dataArray) {
          const categories = molecule.categories ? molecule.categories.join('|') : '';
          const row = [
            molecule.id,
            `"${molecule.name}"`,
            `"${molecule.smiles}"`,
            molecule.dateCreated,
            `"${categories}"`
          ];
          csvRows.push(row.join(','));
        }
        
        dataToExport = csvRows.join('\n');
        filename = selectedMolecule ? `${selectedMolecule.name}.csv` : 'molecules.csv';
        mimeType = 'text/csv';
      } else if (exportFormat === 'smiles') {
        // Export just SMILES strings
        const dataArray = selectedMolecule ? [selectedMolecule] : molecules;
        dataToExport = dataArray.map(m => `${m.smiles} ${m.name}`).join('\n');
        filename = selectedMolecule ? `${selectedMolecule.name}.smi` : 'molecules.smi';
        mimeType = 'text/plain';
      }
      
      // Create and trigger download
      const blob = new Blob([dataToExport], { type: mimeType });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      
      showSnackbar(`Exported ${selectedMolecule ? 'molecule' : 'molecules'} successfully`, 'success');
    } catch (error) {
      console.error('Error exporting molecules:', error);
      showSnackbar('Error exporting molecules', 'error');
    }
  };
  
  const handleFileInputChange = (event) => {
    const file = event.target.files[0];
    if (!file) return;
    
    const reader = new FileReader();
    
    reader.onload = (e) => {
      try {
        const fileContents = e.target.result;
        
        // Determine file type and parse accordingly
        if (file.name.endsWith('.json')) {
          const importedMolecules = JSON.parse(fileContents);
          importMolecules(importedMolecules);
        } else if (file.name.endsWith('.csv')) {
          // Parse CSV
          const lines = fileContents.split('\n');
          const headers = lines[0].split(',');
          const importedMolecules = [];
          
          for (let i = 1; i < lines.length; i++) {
            if (!lines[i].trim()) continue;
            
            const values = lines[i].split(',');
            const molecule = {};
            
            for (let j = 0; j < headers.length; j++) {
              if (headers[j] === 'categories') {
                molecule[headers[j]] = values[j].replace(/"/g, '').split('|').filter(c => c);
              } else {
                molecule[headers[j]] = values[j].replace(/"/g, '');
              }
            }
            
            importedMolecules.push(molecule);
          }
          
          importMolecules(importedMolecules);
        } else if (file.name.endsWith('.smi')) {
          // Parse SMILES file
          const lines = fileContents.split('\n');
          const importedMolecules = lines.map((line, index) => {
            const parts = line.trim().split(/\s+/);
            const smiles = parts[0];
            const name = parts.slice(1).join(' ') || `Imported Molecule ${index + 1}`;
            
            return {
              id: `imported-${Date.now()}-${index}`,
              name,
              smiles,
              dateCreated: new Date().toISOString().split('T')[0],
              categories: ['imported']
            };
          }).filter(m => m.smiles);
          
          importMolecules(importedMolecules);
        } else {
          showSnackbar('Unsupported file format. Please use .json, .csv, or .smi files.', 'error');
        }
      } catch (error) {
        console.error('Error importing file:', error);
        showSnackbar('Error parsing import file', 'error');
      }
    };
    
    reader.onerror = () => {
      showSnackbar('Error reading file', 'error');
    };
    
    if (file.name.endsWith('.json') || file.name.endsWith('.csv') || file.name.endsWith('.smi')) {
      reader.readAsText(file);
    } else {
      showSnackbar('Unsupported file format. Please use .json, .csv, or .smi files.', 'error');
    }
    
    // Reset file input
    event.target.value = null;
  };
  
  const importMolecules = (importedMolecules) => {
    if (!Array.isArray(importedMolecules) || importedMolecules.length === 0) {
      showSnackbar('No valid molecules found in import file', 'error');
      return;
    }
    
    try {
      // Validate each molecule has at least name and smiles
      const validMolecules = importedMolecules.filter(m => m.name && m.smiles);
      
      if (validMolecules.length === 0) {
        showSnackbar('No valid molecules found in import file', 'error');
        return;
      }
      
      // Ensure each molecule has an ID and dateCreated
      const processedMolecules = validMolecules.map(m => ({
        ...m,
        id: m.id || `imported-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`,
        dateCreated: m.dateCreated || new Date().toISOString().split('T')[0],
        categories: m.categories || ['imported']
      }));
      
      // Merge with existing molecules
      const updatedMolecules = [...molecules, ...processedMolecules];
      setMolecules(updatedMolecules);
      
      // Update localStorage
      localStorage.setItem('molecules', JSON.stringify(updatedMolecules));
      
      showSnackbar(`Imported ${processedMolecules.length} molecules successfully`, 'success');
    } catch (error) {
      console.error('Error importing molecules:', error);
      showSnackbar('Error importing molecules', 'error');
    }
  };
  
  const showSnackbar = (message, severity) => {
    setSnackbarMessage(message);
    setSnackbarSeverity(severity);
    setSnackbarOpen(true);
  };
  
  const handleSnackbarClose = () => {
    setSnackbarOpen(false);
  };
  
  const triggerFileInput = () => {
    if (fileInputRef.current) {
      fileInputRef.current.click();
    }
  };
  
  return (
    <div className={classes.root}>
      <Typography variant="h5" className={classes.title}>
        Molecule Database
      </Typography>
      
      <Grid container spacing={3}>
        <Grid item xs={12} md={4}>
          <Paper className={classes.paper}>
            <Typography variant="h6" gutterBottom>
              Stored Molecules
            </Typography>
            
            <TextField
              label="Search Molecules"
              variant="outlined"
              fullWidth
              className={classes.formControl}
              value={searchTerm}
              onChange={handleSearchChange}
              placeholder="Search by name or SMILES"
            />
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel>Filter by Category</InputLabel>
              <Select
                value={categoryFilter}
                onChange={handleCategoryChange}
                label="Filter by Category"
              >
                {categories.map((category) => (
                  <MenuItem key={category} value={category}>
                    {category.charAt(0).toUpperCase() + category.slice(1)}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
            
            {loading ? (
              <div style={{ textAlign: 'center', padding: '20px' }}>
                <CircularProgress />
                <Typography style={{ marginTop: '10px' }}>
                  Loading molecules...
                </Typography>
              </div>
            ) : (
              <div>
                {filteredMolecules.length > 0 ? (
                  <List className={classes.list}>
                    {filteredMolecules.map((molecule) => (
                      <ListItem
                        key={molecule.id}
                        button
                        className={selectedMolecule && selectedMolecule.id === molecule.id ? classes.selectedListItem : classes.listItem}
                        onClick={() => handleSelectMolecule(molecule)}
                      >
                        <ListItemText
                          primary={molecule.name}
                          secondary={
                            <React.Fragment>
                              <Typography variant="body2" component="span" color="textSecondary">
                                {molecule.dateCreated}
                              </Typography>
                              <div style={{ marginTop: '4px' }}>
                                {molecule.categories && molecule.categories.map((category) => (
                                  <Chip
                                    key={category}
                                    label={category}
                                    size="small"
                                    className={classes.chip}
                                  />
                                ))}
                              </div>
                            </React.Fragment>
                          }
                        />
                        <ListItemSecondaryAction>
                          <IconButton
                            edge="end"
                            aria-label="delete"
                            onClick={() => handleDeleteClick(molecule)}
                            size="small"
                          >
                            <DeleteIcon />
                          </IconButton>
                        </ListItemSecondaryAction>
                      </ListItem>
                    ))}
                  </List>
                ) : (
                  <div className={classes.noResults}>
                    <Typography>
                      No molecules found matching your criteria
                    </Typography>
                  </div>
                )}
                
                <Divider className={classes.divider} />
                
                <div className={classes.exportImportButtons}>
                  <div>
                    <FormControl variant="outlined" size="small" style={{ width: '100px', marginRight: '8px' }}>
                      <InputLabel>Format</InputLabel>
                      <Select
                        value={exportFormat}
                        onChange={(e) => setExportFormat(e.target.value)}
                        label="Format"
                        size="small"
                      >
                        <MenuItem value="json">JSON</MenuItem>
                        <MenuItem value="csv">CSV</MenuItem>
                        <MenuItem value="smiles">SMILES</MenuItem>
                      </Select>
                    </FormControl>
                    
                    <Button
                      variant="outlined"
                      color="primary"
                      size="small"
                      startIcon={<GetAppIcon />}
                      onClick={exportMolecules}
                      disabled={molecules.length === 0}
                      className={classes.button}
                    >
                      Export
                    </Button>
                  </div>
                  
                  <input
                    accept=".json,.csv,.smi"
                    className={classes.fileInput}
                    id="import-file"
                    type="file"
                    onChange={handleFileInputChange}
                    ref={fileInputRef}
                  />
                  
                  <Button
                    variant="outlined"
                    color="primary"
                    size="small"
                    startIcon={<PublishIcon />}
                    onClick={triggerFileInput}
                    className={classes.button}
                  >
                    Import
                  </Button>
                </div>
              </div>
            )}
          </Paper>
        </Grid>
        
        <Grid item xs={12} md={8}>
          <Paper className={classes.paper}>
            <Typography variant="h6" gutterBottom>
              Molecule Details
            </Typography>
            
            {selectedMolecule ? (
              <div>
                <Card className={classes.card}>
                  <CardContent>
                    <Grid container spacing={2}>
                      <Grid item xs={12} md={6}>
                        <Typography variant="h6">{selectedMolecule.name}</Typography>
                        <Typography variant="body2" color="textSecondary">
                          ID: {selectedMolecule.id}
                        </Typography>
                        <Typography variant="body2" color="textSecondary" gutterBottom>
                          Created: {selectedMolecule.dateCreated}
                        </Typography>
                        
                        <Typography variant="subtitle2" style={{ marginTop: '16px' }}>
                          Categories:
                        </Typography>
                        <div>
                          {selectedMolecule.categories && selectedMolecule.categories.map((category) => (
                            <Chip
                              key={category}
                              label={category}
                              className={classes.chip}
                            />
                          ))}
                        </div>
                        
                        <Typography variant="subtitle2" style={{ marginTop: '16px' }}>
                          Properties:
                        </Typography>
                        {selectedMolecule.properties && Object.entries(selectedMolecule.properties).map(([key, value]) => (
                          <Typography key={key} variant="body2">
                            <strong>{key}:</strong> {value}
                          </Typography>
                        ))}
                      </Grid>
                      
                      <Grid item xs={12} md={6}>
                        <Typography variant="subtitle2" gutterBottom>
                          SMILES:
                        </Typography>
                        <Typography variant="body2" style={{ wordBreak: 'break-all' }}>
                          {selectedMolecule.smiles}
                        </Typography>
                        
                        <div className={classes.moleculeDisplay}>
                          <MoleculeViewer smiles={selectedMolecule.smiles} height={200} />
                        </div>
                      </Grid>
                    </Grid>
                  </CardContent>
                </Card>
                
                <div>
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<FileCopyIcon />}
                    className={classes.button}
                    onClick={() => onSelectMolecule(selectedMolecule)}
                  >
                    Use This Molecule
                  </Button>
                  
                  <Button
                    variant="outlined"
                    color="primary"
                    startIcon={<EditIcon />}
                    className={classes.button}
                  >
                    Edit Details
                  </Button>
                  
                  <Button
                    variant="outlined"
                    startIcon={<ShareIcon />}
                    className={classes.button}
                  >
                    Share
                  </Button>
                </div>
              </div>
            ) : (
              <div style={{ textAlign: 'center', padding: '40px 20px' }}>
                <Typography color="textSecondary">
                  Select a molecule from the list to view details
                </Typography>
                
                {currentMolecule && currentMolecule.smiles && (
                  <div style={{ marginTop: '30px' }}>
                    <Typography variant="subtitle1" gutterBottom>
                      Currently Working With:
                    </Typography>
                    <Typography variant="h6">
                      {currentMolecule.name || 'Untitled Molecule'}
                    </Typography>
                    <div style={{ maxWidth: '400px', margin: '10px auto' }}>
                      <MoleculeViewer smiles={currentMolecule.smiles} height={150} />
                    </div>
                    <Button
                      variant="contained"
                      color="primary"
                      startIcon={<SaveIcon />}
                      onClick={() => saveMolecule(currentMolecule)}
                      style={{ marginTop: '10px' }}
                    >
                      Save Current Molecule
                    </Button>
                  </div>
                )}
              </div>
            )}
          </Paper>
        </Grid>
      </Grid>
      
      {/* Delete Confirmation Dialog */}
      <Dialog
        open={deleteDialogOpen}
        onClose={cancelDelete}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Confirm Deletion
        </DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            Are you sure you want to delete the molecule "{moleculeToDelete?.name}"? This action cannot be undone.
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={cancelDelete} color="primary">
            Cancel
          </Button>
          <Button onClick={confirmDelete} color="secondary" autoFocus>
            Delete
          </Button>
        </DialogActions>
      </Dialog>
      
      {/* Snackbar for notifications */}
      <Snackbar
        open={snackbarOpen}
        autoHideDuration={4000}
        onClose={handleSnackbarClose}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
      >
        <Alert onClose={handleSnackbarClose} severity={snackbarSeverity}>
          {snackbarMessage}
        </Alert>
      </Snackbar>
    </div>
  );
};

export default MoleculeStorage; 