import React, { useEffect, useRef, useState } from 'react';
import PropTypes from 'prop-types';
import { makeStyles } from '@material-ui/core/styles';
import { Button, Typography, Slider, Grid, Paper, FormControl, InputLabel, Select, MenuItem, IconButton, Tooltip, Box } from '@material-ui/core';
import { FullscreenOutlined, FullscreenExitOutlined, RefreshOutlined, SaveAlt, CameraAlt } from '@material-ui/icons';
import * as $3Dmol from '3dmol/build/3Dmol-min.js';
import 'xterm/css/xterm.css';
import axios from 'axios';

const useStyles = makeStyles((theme) => ({
  viewerContainer: {
    position: 'relative',
    height: '400px',
    width: '100%',
    border: `1px solid ${theme.palette.divider}`,
    borderRadius: theme.shape.borderRadius,
    overflow: 'hidden',
  },
  fullscreenContainer: {
    position: 'fixed',
    top: 0,
    left: 0,
    width: '100vw',
    height: '100vh',
    zIndex: 9999,
    backgroundColor: theme.palette.background.paper,
  },
  controlsContainer: {
    padding: theme.spacing(2),
    marginTop: theme.spacing(2),
  },
  controlRow: {
    display: 'flex',
    alignItems: 'center',
    marginBottom: theme.spacing(2),
  },
  controlLabel: {
    minWidth: '120px',
  },
  viewerControls: {
    position: 'absolute',
    top: theme.spacing(1),
    right: theme.spacing(1),
    zIndex: 100,
    backgroundColor: 'rgba(255,255,255,0.7)',
    borderRadius: theme.shape.borderRadius,
    padding: theme.spacing(0.5),
  },
  slider: {
    width: '100%',
    maxWidth: '300px',
  },
  sliderContainer: {
    flexGrow: 1,
    marginLeft: theme.spacing(1),
    marginRight: theme.spacing(1),
  },
  formControl: {
    minWidth: 120,
  },
}));

const MoleculeViewer3D = ({
  moleculeData,
  format = 'smiles',
  viewerType = '3dmol',
  initialStyle = 'stick',
  showControls = true,
  height = 400,
  width = '100%'
}) => {
  const classes = useStyles();
  const viewerRef = useRef(null);
  const containerRef = useRef(null);
  const [viewer, setViewer] = useState(null);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [modelStyle, setModelStyle] = useState(initialStyle);
  const [surfaceOpacity, setSurfaceOpacity] = useState(0.5);
  const [surfaceType, setSurfaceType] = useState('VDW');
  const [backgroundColor, setBackgroundColor] = useState('#ffffff');
  const [rotationSpeed, setRotationSpeed] = useState(0);
  const [isSurfaceVisible, setIsSurfaceVisible] = useState(false);
  const [isLoading, setIsLoading] = useState(true);
  const [molData, setMolData] = useState(null);
  const [error, setError] = useState(null);
  const animationRef = useRef(null);

  // Fetch or convert molecule data whenever moleculeData or format changes
  useEffect(() => {
    setIsLoading(true);
    setError(null);
    setMolData(null); // Clear previous data

    if (!moleculeData) {
      setError('No molecule data provided');
      setIsLoading(false);
      return;
    }

    const fetchMoleculeData = async () => {
      try {
        let data = moleculeData;
        let currentFormat = format;

        // If SMILES is provided, convert it to 3D structure via our backend API
        if (format === 'smiles') {
          const response = await axios.post('/api/simulation/3d-structure', { 
            smiles: moleculeData 
          });
          
          if (response.data && response.data.molblock) {
            data = response.data.molblock;
            currentFormat = 'mol'; // The backend returns molblock format
          } else {
            throw new Error('Failed to convert SMILES to 3D structure');
          }
        }

        setMolData({ data: data, format: currentFormat });
        setIsLoading(false);
      } catch (err) {
        console.error('Error fetching molecule data:', err);
        setError(err.message || 'Failed to load molecule data');
        setIsLoading(false);
      }
    };

    fetchMoleculeData();
  // Add viewer to dependency array to re-run if viewer instance changes
  }, [moleculeData, format, viewer]); 

  // Initialize or update the 3Dmol.js viewer when molData or background changes
  useEffect(() => {
    if (!viewerRef.current || isLoading || !molData) return;

    let currentViewer = viewer;

    // Initialize viewer if it doesn't exist
    if (!currentViewer) {
      const viewerConfig = { backgroundColor };
      currentViewer = $3Dmol.createViewer(viewerRef.current, viewerConfig);
      setViewer(currentViewer);
    } else {
      // Clear existing models if updating
      currentViewer.removeAllModels();
      currentViewer.removeAllSurfaces();
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
        animationRef.current = null;
      }
      // Update background if it changed
      currentViewer.setBackgroundColor(backgroundColor); 
    }

    // Load the new molecule data
    try {
      let model;
      if (molData.format === 'pdb') {
        model = currentViewer.addModel(molData.data, 'pdb');
      } else if (molData.format === 'mol' || molData.format === 'sdf') {
        model = currentViewer.addModel(molData.data, 'mol');
      } else {
        // Attempt to add with the provided format, fallback or error handling might be needed
        model = currentViewer.addModel(molData.data, molData.format);
      }
      
      // Apply initial style
      applyStyle(currentViewer, model, modelStyle);
      
      // If surface is enabled, apply it
      if (isSurfaceVisible) {
        applySurface(currentViewer, model, surfaceType, surfaceOpacity);
      }
      
      currentViewer.zoomTo();
      currentViewer.render();
      
      // Restart animation if rotationSpeed > 0
      if (rotationSpeed > 0) {
        startAnimation(currentViewer, rotationSpeed);
      }
    } catch (err) {
      console.error('Error rendering molecule:', err);
      setError(`Failed to render molecule: ${err.message}`);
    }

    // Cleanup function for animation
    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  // Add dependencies that should trigger re-initialization or update
  }, [molData, isLoading, backgroundColor]);

  // Apply style changes
  useEffect(() => {
    if (!viewer) return;
    
    try {
      viewer.removeAllSurfaces();
      applyStyle(viewer, null, modelStyle);
      
      if (isSurfaceVisible) {
        applySurface(viewer, null, surfaceType, surfaceOpacity);
      }
      
      viewer.render();
    } catch (err) {
      console.error('Error applying style:', err);
    }
  }, [viewer, modelStyle, isSurfaceVisible, surfaceType, surfaceOpacity]);

  // Handle background color changes
  useEffect(() => {
    if (!viewer) return;
    viewer.setBackgroundColor(backgroundColor);
    viewer.render();
  }, [backgroundColor, viewer]);

  // Handle rotation speed changes
  useEffect(() => {
    if (!viewer) return;
    
    if (animationRef.current) {
      cancelAnimationFrame(animationRef.current);
      animationRef.current = null;
    }
    
    if (rotationSpeed > 0) {
      startAnimation(viewer, rotationSpeed);
    }
  }, [rotationSpeed, viewer]);

  // Apply style to the model
  const applyStyle = (viewer, model, style) => {
    if (!viewer) return;
    
    // Apply to all models if model is not specified
    const selector = model ? { model } : {};
    
    viewer.removeAllStyles();
    
    switch (style) {
      case 'stick':
        viewer.setStyle(selector, { stick: {} });
        break;
      case 'line':
        viewer.setStyle(selector, { line: {} });
        break;
      case 'cross':
        viewer.setStyle(selector, { cross: { linewidth: 5 } });
        break;
      case 'sphere':
        viewer.setStyle(selector, { sphere: { scale: 0.3 } });
        break;
      case 'cartoon':
        viewer.setStyle(selector, { cartoon: {} });
        break;
      case 'ballAndStick':
        viewer.setStyle(selector, { 
          stick: { radius: 0.15 }, 
          sphere: { scale: 0.3 } 
        });
        break;
      default:
        viewer.setStyle(selector, { stick: {} });
    }
  };

  // Apply surface to the model
  const applySurface = (viewer, model, surfaceType, opacity) => {
    if (!viewer) return;
    
    const selector = model ? { model } : {};
    const surfaceOptions = { opacity };
    
    switch (surfaceType) {
      case 'VDW':
        viewer.addSurface($3Dmol.SurfaceType.VDW, surfaceOptions, selector);
        break;
      case 'SAS':
        viewer.addSurface($3Dmol.SurfaceType.SAS, surfaceOptions, selector);
        break;
      case 'SES':
        viewer.addSurface($3Dmol.SurfaceType.SES, surfaceOptions, selector);
        break;
      case 'MS':
        viewer.addSurface($3Dmol.SurfaceType.MS, surfaceOptions, selector);
        break;
      default:
        viewer.addSurface($3Dmol.SurfaceType.VDW, surfaceOptions, selector);
    }
  };

  // Start rotation animation
  const startAnimation = (viewer, speed) => {
    if (!viewer) return;
    
    const rotate = () => {
      viewer.rotate(1 * (speed / 10), { x: 0, y: 1, z: 0 });
      viewer.render();
      animationRef.current = requestAnimationFrame(rotate);
    };
    
    animationRef.current = requestAnimationFrame(rotate);
  };

  // Toggle fullscreen mode
  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
  };

  // Reset view
  const resetView = () => {
    if (viewer) {
      viewer.zoomTo();
      viewer.render();
    }
  };

  // Take screenshot
  const takeScreenshot = () => {
    if (!viewer) return;
    
    try {
      const imgData = viewer.pngURI();
      const a = document.createElement('a');
      a.href = imgData;
      a.download = 'molecule.png';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } catch (err) {
      console.error('Error taking screenshot:', err);
    }
  };

  // Handle style change
  const handleStyleChange = (event) => {
    setModelStyle(event.target.value);
  };

  // Handle surface type change
  const handleSurfaceTypeChange = (event) => {
    setSurfaceType(event.target.value);
  };

  // Toggle surface visibility
  const toggleSurface = () => {
    setIsSurfaceVisible(!isSurfaceVisible);
  };

  return (
    <div>
      <div 
        ref={containerRef}
        className={isFullscreen ? classes.fullscreenContainer : classes.viewerContainer}
        style={!isFullscreen ? { height: `${height}px`, width } : {}}
      >
        {isLoading ? (
          <Box display="flex" justifyContent="center" alignItems="center" height="100%">
            <Typography variant="h6">Loading molecule...</Typography>
          </Box>
        ) : error ? (
          <Box display="flex" justifyContent="center" alignItems="center" height="100%">
            <Typography variant="h6" color="error">{error}</Typography>
          </Box>
        ) : (
          <>
            <div className={classes.viewerControls}>
              <Tooltip title={isFullscreen ? "Exit Fullscreen" : "Fullscreen"}>
                <IconButton size="small" onClick={toggleFullscreen}>
                  {isFullscreen ? <FullscreenExitOutlined /> : <FullscreenOutlined />}
                </IconButton>
              </Tooltip>
              <Tooltip title="Reset View">
                <IconButton size="small" onClick={resetView}>
                  <RefreshOutlined />
                </IconButton>
              </Tooltip>
              <Tooltip title="Save Image">
                <IconButton size="small" onClick={takeScreenshot}>
                  <CameraAlt />
                </IconButton>
              </Tooltip>
            </div>
            <div 
              ref={viewerRef} 
              style={{ width: '100%', height: '100%' }}
              data-testid="molecule-viewer-3d"
            ></div>
          </>
        )}
      </div>

      {showControls && !isLoading && !error && (
        <Paper className={classes.controlsContainer}>
          <Grid container spacing={3}>
            <Grid item xs={12} sm={6}>
              <div className={classes.controlRow}>
                <Typography className={classes.controlLabel}>
                  Representation:
                </Typography>
                <FormControl className={classes.formControl}>
                  <Select
                    value={modelStyle}
                    onChange={handleStyleChange}
                    displayEmpty
                  >
                    <MenuItem value="stick">Stick</MenuItem>
                    <MenuItem value="line">Line</MenuItem>
                    <MenuItem value="cross">Cross</MenuItem>
                    <MenuItem value="sphere">Sphere</MenuItem>
                    <MenuItem value="ballAndStick">Ball and Stick</MenuItem>
                    <MenuItem value="cartoon">Cartoon</MenuItem>
                  </Select>
                </FormControl>
              </div>
              
              <div className={classes.controlRow}>
                <Typography className={classes.controlLabel}>
                  Surface:
                </Typography>
                <Button 
                  variant={isSurfaceVisible ? "contained" : "outlined"} 
                  color={isSurfaceVisible ? "primary" : "default"}
                  onClick={toggleSurface}
                >
                  {isSurfaceVisible ? "Hide Surface" : "Show Surface"}
                </Button>
              </div>
              
              {isSurfaceVisible && (
                <>
                  <div className={classes.controlRow}>
                    <Typography className={classes.controlLabel}>
                      Surface Type:
                    </Typography>
                    <FormControl className={classes.formControl}>
                      <Select
                        value={surfaceType}
                        onChange={handleSurfaceTypeChange}
                      >
                        <MenuItem value="VDW">Van der Waals</MenuItem>
                        <MenuItem value="SAS">Solvent Accessible</MenuItem>
                        <MenuItem value="SES">Solvent Excluded</MenuItem>
                        <MenuItem value="MS">Molecular Surface</MenuItem>
                      </Select>
                    </FormControl>
                  </div>
                  
                  <div className={classes.controlRow}>
                    <Typography className={classes.controlLabel}>
                      Opacity:
                    </Typography>
                    <div className={classes.sliderContainer}>
                      <Slider
                        value={surfaceOpacity}
                        min={0}
                        max={1}
                        step={0.1}
                        onChange={(_, newValue) => setSurfaceOpacity(newValue)}
                        aria-labelledby="surface-opacity-slider"
                        className={classes.slider}
                      />
                    </div>
                    <Typography>{surfaceOpacity.toFixed(1)}</Typography>
                  </div>
                </>
              )}
            </Grid>

            <Grid item xs={12} sm={6}>
              <div className={classes.controlRow}>
                <Typography className={classes.controlLabel}>
                  Background:
                </Typography>
                <FormControl className={classes.formControl}>
                  <Select
                    value={backgroundColor}
                    onChange={(e) => setBackgroundColor(e.target.value)}
                  >
                    <MenuItem value="#ffffff">White</MenuItem>
                    <MenuItem value="#000000">Black</MenuItem>
                    <MenuItem value="#f5f5f5">Light Gray</MenuItem>
                    <MenuItem value="#121212">Dark Gray</MenuItem>
                    <MenuItem value="#e3f2fd">Light Blue</MenuItem>
                  </Select>
                </FormControl>
              </div>
              
              <div className={classes.controlRow}>
                <Typography className={classes.controlLabel}>
                  Auto-Rotate:
                </Typography>
                <div className={classes.sliderContainer}>
                  <Slider
                    value={rotationSpeed}
                    min={0}
                    max={10}
                    step={1}
                    onChange={(_, newValue) => setRotationSpeed(newValue)}
                    aria-labelledby="rotation-speed-slider"
                    className={classes.slider}
                  />
                </div>
                <Typography>{rotationSpeed}</Typography>
              </div>
            </Grid>
          </Grid>
        </Paper>
      )}
    </div>
  );
};

MoleculeViewer3D.propTypes = {
  moleculeData: PropTypes.string,
  format: PropTypes.oneOf(['smiles', 'pdb', 'mol2', 'sdf', 'cif', 'xyz']),
  viewerType: PropTypes.string,
  initialStyle: PropTypes.oneOf(['stick', 'line', 'cross', 'sphere', 'cartoon', 'ballAndStick']),
  showControls: PropTypes.bool,
  height: PropTypes.oneOfType([PropTypes.number, PropTypes.string]),
  width: PropTypes.oneOfType([PropTypes.number, PropTypes.string]),
};

export default MoleculeViewer3D; 