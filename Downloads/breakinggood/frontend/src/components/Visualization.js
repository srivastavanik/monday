import React, { useState, useEffect, useRef } from 'react';
import { 
  Typography, 
  Paper, 
  Grid, 
  Button, 
  FormControl, 
  InputLabel, 
  Select, 
  MenuItem, 
  Slider,
  Tabs,
  Tab,
  Box,
  Checkbox,
  FormControlLabel,
  Divider,
  makeStyles,
  CircularProgress,
  IconButton
} from '@material-ui/core';
import GetAppIcon from '@material-ui/icons/GetApp';
import FullscreenIcon from '@material-ui/icons/Fullscreen';
import FullscreenExitIcon from '@material-ui/icons/FullscreenExit';
import RotateLeftIcon from '@material-ui/icons/RotateLeft';
import ZoomInIcon from '@material-ui/icons/ZoomIn';
import ZoomOutIcon from '@material-ui/icons/ZoomOut';
import ColorLensIcon from '@material-ui/icons/ColorLens';
import Alert from '@material-ui/lab/Alert';

// This would be implemented with 3Dmol.js, RDKit, or similar in a real application
// Here we'll simulate the viewer with a placeholder
const MoleculeViewer3D = ({ 
  smiles, 
  representation, 
  colors, 
  height,
  showHydrogens,
  showLabels,
  rotationSpeed
}) => {
  const containerRef = useRef(null);
  const [loading, setLoading] = useState(true);
  
  useEffect(() => {
    if (!smiles || !containerRef.current) return;
    
    // Simulate loading of molecule
    setLoading(true);
    const timer = setTimeout(() => {
      setLoading(false);
      
      // In a real implementation, this would initialize 3Dmol.js or another library
      // and render the molecule with the specified options
      const container = containerRef.current;
      
      // Clear previous content
      while (container.firstChild) {
        container.removeChild(container.firstChild);
      }
      
      // Create placeholder visualization
      const placeholder = document.createElement('div');
      placeholder.style.width = '100%';
      placeholder.style.height = '100%';
      placeholder.style.display = 'flex';
      placeholder.style.flexDirection = 'column';
      placeholder.style.alignItems = 'center';
      placeholder.style.justifyContent = 'center';
      placeholder.style.backgroundColor = '#f5f5f5';
      placeholder.style.position = 'relative';
      
      const info = document.createElement('div');
      info.innerHTML = `
        <div style="text-align: center; padding: 0 20px;">
          <strong>SMILES:</strong> ${smiles}<br><br>
          <strong>Visualization Options:</strong><br>
          Representation: ${representation}<br>
          Color Scheme: ${colors}<br>
          Show Hydrogens: ${showHydrogens ? 'Yes' : 'No'}<br>
          Show Labels: ${showLabels ? 'Yes' : 'No'}<br>
          Rotation Speed: ${rotationSpeed}
        </div>
      `;
      
      // Simulated molecule representation
      const canvas = document.createElement('canvas');
      canvas.width = container.clientWidth;
      canvas.height = 200;
      canvas.style.marginBottom = '20px';
      
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#e3f2fd';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      
      // Draw a simulated molecule based on representation
      ctx.strokeStyle = colors === 'element' ? '#4285F4' : (colors === 'chain' ? '#DB4437' : '#0F9D58');
      ctx.lineWidth = 2;
      ctx.beginPath();
      
      if (representation === 'ball-and-stick') {
        // Draw ball and stick representation
        const centerX = canvas.width / 2;
        const centerY = canvas.height / 2;
        const radius = 20;
        
        // Draw a simple molecule structure
        ctx.moveTo(centerX - 50, centerY);
        ctx.lineTo(centerX, centerY);
        ctx.lineTo(centerX + 50, centerY);
        
        ctx.moveTo(centerX, centerY);
        ctx.lineTo(centerX, centerY - 50);
        
        ctx.stroke();
        
        // Draw "atoms"
        ctx.beginPath();
        ctx.arc(centerX - 50, centerY, radius, 0, Math.PI * 2);
        ctx.fillStyle = '#DB4437';
        ctx.fill();
        
        ctx.beginPath();
        ctx.arc(centerX, centerY, radius, 0, Math.PI * 2);
        ctx.fillStyle = '#4285F4';
        ctx.fill();
        
        ctx.beginPath();
        ctx.arc(centerX + 50, centerY, radius, 0, Math.PI * 2);
        ctx.fillStyle = '#0F9D58';
        ctx.fill();
        
        ctx.beginPath();
        ctx.arc(centerX, centerY - 50, radius, 0, Math.PI * 2);
        ctx.fillStyle = '#F4B400';
        ctx.fill();
        
        // Add labels if enabled
        if (showLabels) {
          ctx.fillStyle = 'black';
          ctx.font = '12px Arial';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';
          ctx.fillText('C', centerX - 50, centerY);
          ctx.fillText('N', centerX, centerY);
          ctx.fillText('O', centerX + 50, centerY);
          ctx.fillText('H', centerX, centerY - 50);
        }
      } else if (representation === 'stick') {
        // Draw stick representation
        const centerX = canvas.width / 2;
        const centerY = canvas.height / 2;
        
        ctx.lineWidth = 8;
        
        // Draw a simple molecule structure
        ctx.moveTo(centerX - 60, centerY);
        ctx.lineTo(centerX, centerY);
        ctx.lineTo(centerX + 60, centerY);
        
        ctx.moveTo(centerX, centerY);
        ctx.lineTo(centerX, centerY - 60);
        
        ctx.stroke();
        
        // Add labels if enabled
        if (showLabels) {
          ctx.fillStyle = 'black';
          ctx.font = '12px Arial';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';
          ctx.fillText('C', centerX - 60, centerY - 15);
          ctx.fillText('N', centerX, centerY - 15);
          ctx.fillText('O', centerX + 60, centerY - 15);
          ctx.fillText('H', centerX + 15, centerY - 60);
        }
      } else if (representation === 'surface') {
        // Draw molecular surface representation
        const centerX = canvas.width / 2;
        const centerY = canvas.height / 2;
        
        // Draw a blob for the surface
        ctx.beginPath();
        ctx.ellipse(centerX, centerY, 100, 70, 0, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(66, 133, 244, 0.5)';
        ctx.fill();
        ctx.strokeStyle = 'rgba(66, 133, 244, 0.8)';
        ctx.lineWidth = 2;
        ctx.stroke();
      }
      
      placeholder.appendChild(canvas);
      placeholder.appendChild(info);
      
      container.appendChild(placeholder);
      
      // Simulate rotation based on speed
      if (rotationSpeed > 0) {
        const rotationIndicator = document.createElement('div');
        rotationIndicator.style.position = 'absolute';
        rotationIndicator.style.top = '10px';
        rotationIndicator.style.right = '10px';
        rotationIndicator.style.backgroundColor = 'rgba(0, 0, 0, 0.2)';
        rotationIndicator.style.color = 'white';
        rotationIndicator.style.padding = '5px 10px';
        rotationIndicator.style.borderRadius = '3px';
        rotationIndicator.textContent = `Rotating (${rotationSpeed}x)`;
        
        placeholder.appendChild(rotationIndicator);
      }
    }, 1000);
    
    return () => clearTimeout(timer);
  }, [smiles, representation, colors, showHydrogens, showLabels, rotationSpeed]);
  
  return (
    <div 
      ref={containerRef} 
      style={{
        width: '100%',
        height: height || 400,
        border: '1px solid #e0e0e0',
        position: 'relative'
      }}
    >
      {loading && (
        <div style={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          backgroundColor: 'rgba(255, 255, 255, 0.7)'
        }}>
          <CircularProgress />
        </div>
      )}
    </div>
  );
};

// This would be implemented with a 2D renderer like RDKit in a real application
const MoleculeViewer2D = ({ smiles, colors, showHydrogens, showLabels, height }) => {
  const containerRef = useRef(null);
  const [loading, setLoading] = useState(true);
  
  useEffect(() => {
    if (!smiles || !containerRef.current) return;
    
    // Simulate loading of molecule
    setLoading(true);
    const timer = setTimeout(() => {
      setLoading(false);
      
      // In a real implementation, this would initialize a 2D renderer
      const container = containerRef.current;
      
      // Clear previous content
      while (container.firstChild) {
        container.removeChild(container.firstChild);
      }
      
      // Create placeholder visualization
      const placeholder = document.createElement('div');
      placeholder.style.width = '100%';
      placeholder.style.height = '100%';
      placeholder.style.display = 'flex';
      placeholder.style.flexDirection = 'column';
      placeholder.style.alignItems = 'center';
      placeholder.style.justifyContent = 'center';
      placeholder.style.backgroundColor = '#f5f5f5';
      
      const info = document.createElement('div');
      info.innerHTML = `
        <div style="text-align: center; padding: 0 20px;">
          <strong>SMILES:</strong> ${smiles}<br><br>
          <strong>Visualization Options:</strong><br>
          Color Scheme: ${colors}<br>
          Show Hydrogens: ${showHydrogens ? 'Yes' : 'No'}<br>
          Show Labels: ${showLabels ? 'Yes' : 'No'}
        </div>
      `;
      
      // Simulated molecule representation
      const canvas = document.createElement('canvas');
      canvas.width = container.clientWidth;
      canvas.height = 200;
      canvas.style.marginBottom = '20px';
      
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#fff';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      
      // Draw a simulated 2D molecule
      const centerX = canvas.width / 2;
      const centerY = canvas.height / 2;
      
      // Draw a hexagon for benzene
      ctx.beginPath();
      ctx.strokeStyle = colors === 'element' ? 'black' : (colors === 'chain' ? '#DB4437' : '#4285F4');
      ctx.lineWidth = 2;
      
      const radius = 60;
      for (let i = 0; i < 6; i++) {
        const angle = (i * Math.PI) / 3;
        const x = centerX + radius * Math.cos(angle);
        const y = centerY + radius * Math.sin(angle);
        
        if (i === 0) {
          ctx.moveTo(x, y);
        } else {
          ctx.lineTo(x, y);
        }
      }
      ctx.closePath();
      ctx.stroke();
      
      // Draw the double bonds alternating
      ctx.beginPath();
      for (let i = 0; i < 3; i++) {
        const angle1 = ((i * 2) * Math.PI) / 3;
        const angle2 = ((i * 2 + 1) * Math.PI) / 3;
        
        const x1 = centerX + (radius - 8) * Math.cos(angle1);
        const y1 = centerY + (radius - 8) * Math.sin(angle1);
        const x2 = centerX + (radius - 8) * Math.cos(angle2);
        const y2 = centerY + (radius - 8) * Math.sin(angle2);
        
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
      }
      ctx.stroke();
      
      // Add labels if enabled
      if (showLabels) {
        ctx.fillStyle = 'black';
        ctx.font = '14px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        
        for (let i = 0; i < 6; i++) {
          const angle = (i * Math.PI) / 3;
          const x = centerX + (radius + 20) * Math.cos(angle);
          const y = centerY + (radius + 20) * Math.sin(angle);
          
          ctx.fillText('C', x, y);
          
          // Draw hydrogen labels if enabled
          if (showHydrogens) {
            const hAngle = angle + Math.PI / 6;
            const hx = x + 15 * Math.cos(hAngle);
            const hy = y + 15 * Math.sin(hAngle);
            
            ctx.fillText('H', hx, hy);
          }
        }
      }
      
      placeholder.appendChild(canvas);
      placeholder.appendChild(info);
      
      container.appendChild(placeholder);
    }, 1000);
    
    return () => clearTimeout(timer);
  }, [smiles, colors, showHydrogens, showLabels]);
  
  return (
    <div 
      ref={containerRef} 
      style={{
        width: '100%',
        height: height || 400,
        border: '1px solid #e0e0e0',
        position: 'relative'
      }}
    >
      {loading && (
        <div style={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          backgroundColor: 'rgba(255, 255, 255, 0.7)'
        }}>
          <CircularProgress />
        </div>
      )}
    </div>
  );
};

// TabPanel component for tabbed interface
function TabPanel(props) {
  const { children, value, index, ...other } = props;

  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`visualization-tabpanel-${index}`}
      aria-labelledby={`visualization-tab-${index}`}
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
  paper: {
    padding: theme.spacing(3),
  },
  title: {
    marginBottom: theme.spacing(3),
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120,
  },
  slider: {
    width: '90%',
    margin: theme.spacing(0, 'auto'),
  },
  controlsSection: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2),
  },
  viewerContainer: {
    position: 'relative',
    marginTop: theme.spacing(3),
  },
  viewerControls: {
    position: 'absolute',
    top: 10,
    right: 10,
    zIndex: 1,
    backgroundColor: 'rgba(255, 255, 255, 0.7)',
    borderRadius: theme.shape.borderRadius,
    padding: theme.spacing(0.5),
  },
  fullscreenViewer: {
    position: 'fixed',
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    zIndex: 1300,
    backgroundColor: 'white',
    padding: theme.spacing(2),
  },
  fullscreenClose: {
    position: 'absolute',
    top: 10,
    right: 10,
    zIndex: 1301,
  },
  divider: {
    margin: theme.spacing(2, 0),
  },
  section: {
    marginBottom: theme.spacing(3),
  },
  exportButton: {
    marginTop: theme.spacing(2),
  },
}));

const Visualization = ({ molecule }) => {
  const classes = useStyles();
  const [tabValue, setTabValue] = useState(0);
  const [representation, setRepresentation] = useState('ball-and-stick');
  const [colorScheme, setColorScheme] = useState('element');
  const [showHydrogens, setShowHydrogens] = useState(true);
  const [showLabels, setShowLabels] = useState(true);
  const [rotationSpeed, setRotationSpeed] = useState(1);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [viewerHeight, setViewerHeight] = useState(400);
  const [error, setError] = useState(null);
  
  // Example molecule if none provided
  const defaultMolecule = {
    name: 'Methylphenidate',
    smiles: 'CN(C)C(C1=CC=CC=C1)C(C)OC(=O)C'
  };
  
  const currentMolecule = molecule || defaultMolecule;
  
  const handleTabChange = (event, newValue) => {
    setTabValue(newValue);
  };
  
  const handleRepresentationChange = (event) => {
    setRepresentation(event.target.value);
  };
  
  const handleColorSchemeChange = (event) => {
    setColorScheme(event.target.value);
  };
  
  const handleShowHydrogensChange = (event) => {
    setShowHydrogens(event.target.checked);
  };
  
  const handleShowLabelsChange = (event) => {
    setShowLabels(event.target.checked);
  };
  
  const handleRotationSpeedChange = (event, newValue) => {
    setRotationSpeed(newValue);
  };
  
  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
    if (!isFullscreen) {
      setViewerHeight(window.innerHeight - 100);
    } else {
      setViewerHeight(400);
    }
  };
  
  const resetView = () => {
    // In a real implementation, this would reset the view of the 3D molecule
    // Here we'll just indicate that it was called
    console.log('Reset view called');
    setError(null);
  };
  
  const zoomIn = () => {
    // In a real implementation, this would zoom in on the molecule
    // Here we'll just indicate that it was called
    console.log('Zoom in called');
    setError(null);
  };
  
  const zoomOut = () => {
    // In a real implementation, this would zoom out on the molecule
    // Here we'll just indicate that it was called
    console.log('Zoom out called');
    setError(null);
  };
  
  const exportImage = () => {
    try {
      // In a real implementation, this would export the current view as an image
      // Here we'll just simulate a download
      
      // Create a placeholder download
      const link = document.createElement('a');
      link.download = `${currentMolecule.name || 'molecule'}_visualization.png`;
      link.href = '#';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      
      console.log('Export called for', currentMolecule.name);
      setError(null);
    } catch (err) {
      console.error('Error exporting image:', err);
      setError('Error exporting image. Please try again.');
    }
  };
  
  return (
    <div className={isFullscreen ? classes.fullscreenViewer : classes.root}>
      {isFullscreen && (
        <div className={classes.fullscreenClose}>
          <Button
            variant="contained"
            color="secondary"
            onClick={toggleFullscreen}
            startIcon={<FullscreenExitIcon />}
          >
            Exit Fullscreen
          </Button>
        </div>
      )}
      
      {!isFullscreen && (
        <Typography variant="h5" className={classes.title}>
          Molecule Visualization
        </Typography>
      )}
      
      <Paper className={classes.paper}>
        <Typography variant="h6" gutterBottom>
          {currentMolecule.name || 'Unnamed Molecule'}
        </Typography>
        
        {!isFullscreen && (
          <Tabs
            value={tabValue}
            onChange={handleTabChange}
            indicatorColor="primary"
            textColor="primary"
            variant="fullWidth"
          >
            <Tab label="3D View" />
            <Tab label="2D View" />
          </Tabs>
        )}
        
        {error && (
          <Alert severity="error" style={{ marginTop: 16, marginBottom: 16 }}>
            {error}
          </Alert>
        )}
        
        <div className={classes.viewerContainer}>
          <div className={classes.viewerControls}>
            <IconButton onClick={resetView} title="Reset View">
              <RotateLeftIcon />
            </IconButton>
            <IconButton onClick={zoomIn} title="Zoom In">
              <ZoomInIcon />
            </IconButton>
            <IconButton onClick={zoomOut} title="Zoom Out">
              <ZoomOutIcon />
            </IconButton>
            <IconButton onClick={toggleFullscreen} title="Toggle Fullscreen">
              {isFullscreen ? <FullscreenExitIcon /> : <FullscreenIcon />}
            </IconButton>
          </div>
          
          {isFullscreen ? (
            tabValue === 0 ? (
              <MoleculeViewer3D
                smiles={currentMolecule.smiles}
                representation={representation}
                colors={colorScheme}
                showHydrogens={showHydrogens}
                showLabels={showLabels}
                rotationSpeed={rotationSpeed}
                height={viewerHeight}
              />
            ) : (
              <MoleculeViewer2D
                smiles={currentMolecule.smiles}
                colors={colorScheme}
                showHydrogens={showHydrogens}
                showLabels={showLabels}
                height={viewerHeight}
              />
            )
          ) : (
            <TabPanel value={tabValue} index={0}>
              <MoleculeViewer3D
                smiles={currentMolecule.smiles}
                representation={representation}
                colors={colorScheme}
                showHydrogens={showHydrogens}
                showLabels={showLabels}
                rotationSpeed={rotationSpeed}
                height={viewerHeight}
              />
            </TabPanel>
          )}
          
          {!isFullscreen && (
            <TabPanel value={tabValue} index={1}>
              <MoleculeViewer2D
                smiles={currentMolecule.smiles}
                colors={colorScheme}
                showHydrogens={showHydrogens}
                showLabels={showLabels}
                height={viewerHeight}
              />
            </TabPanel>
          )}
        </div>
        
        <Divider className={classes.divider} />
        
        <div className={classes.controlsSection}>
          <Typography variant="subtitle1" gutterBottom>
            Visualization Options
          </Typography>
          
          <Grid container spacing={2}>
            {tabValue === 0 && (
              <Grid item xs={12} sm={4}>
                <FormControl className={classes.formControl}>
                  <InputLabel>Representation</InputLabel>
                  <Select
                    value={representation}
                    onChange={handleRepresentationChange}
                  >
                    <MenuItem value="ball-and-stick">Ball and Stick</MenuItem>
                    <MenuItem value="stick">Stick</MenuItem>
                    <MenuItem value="surface">Surface</MenuItem>
                  </Select>
                </FormControl>
              </Grid>
            )}
            
            <Grid item xs={12} sm={4}>
              <FormControl className={classes.formControl}>
                <InputLabel>Color Scheme</InputLabel>
                <Select
                  value={colorScheme}
                  onChange={handleColorSchemeChange}
                >
                  <MenuItem value="element">Element</MenuItem>
                  <MenuItem value="chain">Chain</MenuItem>
                  <MenuItem value="residue">Residue</MenuItem>
                </Select>
              </FormControl>
            </Grid>
            
            <Grid item xs={12} sm={4}>
              <FormControlLabel
                control={
                  <Checkbox
                    checked={showHydrogens}
                    onChange={handleShowHydrogensChange}
                    color="primary"
                  />
                }
                label="Show Hydrogens"
              />
              
              <FormControlLabel
                control={
                  <Checkbox
                    checked={showLabels}
                    onChange={handleShowLabelsChange}
                    color="primary"
                  />
                }
                label="Show Labels"
              />
            </Grid>
            
            {tabValue === 0 && (
              <Grid item xs={12}>
                <Typography id="rotation-speed-slider" gutterBottom>
                  Rotation Speed: {rotationSpeed}x
                </Typography>
                <Slider
                  value={rotationSpeed}
                  onChange={handleRotationSpeedChange}
                  aria-labelledby="rotation-speed-slider"
                  min={0}
                  max={5}
                  step={1}
                  marks
                  className={classes.slider}
                />
              </Grid>
            )}
          </Grid>
          
          <Button
            variant="contained"
            color="primary"
            startIcon={<GetAppIcon />}
            onClick={exportImage}
            className={classes.exportButton}
          >
            Export as Image
          </Button>
        </div>
      </Paper>
    </div>
  );
};

export default Visualization; 