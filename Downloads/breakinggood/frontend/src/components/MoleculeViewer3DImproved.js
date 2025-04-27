import React, { useEffect, useRef, useState } from 'react';
import PropTypes from 'prop-types';
import { makeStyles } from '@material-ui/core/styles';
import { Button, Typography, Slider, Grid, Paper, FormControl, InputLabel, Select, MenuItem, IconButton, Tooltip, Box, CircularProgress } from '@material-ui/core';
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
  loadingOverlay: {
    position: 'absolute',
    top: 0,
    left: 0,
    width: '100%',
    height: '100%',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    backgroundColor: 'rgba(255,255,255,0.7)',
    zIndex: 10
  },
  errorOverlay: {
    position: 'absolute',
    top: 0,
    left: 0,
    width: '100%',
    height: '100%',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    backgroundColor: 'rgba(0,0,0,0.1)',
    zIndex: 10,
    padding: theme.spacing(2),
    textAlign: 'center'
  },
  viewerControls: {
    position: 'absolute',
    top: theme.spacing(1),
    right: theme.spacing(1),
    zIndex: 5,
    backgroundColor: 'rgba(255,255,255,0.7)',
    borderRadius: theme.shape.borderRadius,
    padding: theme.spacing(0.5),
  },
}));

const MoleculeViewer3DImproved = (props) => {
  const {
    smiles,               // For backward compatibility
    moleculeData,         // Main molecule data input
    format = 'smiles',    // Default to SMILES format
    initialStyle = 'stick',
    height = 400,
    width = '100%'
  } = props;
  
  const classes = useStyles();
  const viewerRef = useRef(null);
  const [viewer, setViewer] = useState(null);
  const [modelStyle, setModelStyle] = useState(initialStyle);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);
  const [molData, setMolData] = useState(null);
  const animationRef = useRef(null);
  
  // Support both smiles and moleculeData props
  const moleculeInput = smiles || moleculeData;
  
  // Fetch 3D structure from SMILES when input changes
  useEffect(() => {
    if (!moleculeInput) {
      setError('No molecule data provided');
      setIsLoading(false);
      return;
    }
    
    const fetchMoleculeData = async () => {
      setIsLoading(true);
      setError(null);
      
      try {
        // Only call API if input is SMILES format
        if (format === 'smiles') {
          // Make sure we're using the correct API endpoint
          const response = await axios.post('/api/simulation/3d-structure', { 
            smiles: moleculeInput 
          });
          
          if (response.data && response.data.molblock) {
            setMolData(response.data.molblock);
          } else {
            throw new Error('Failed to convert SMILES to 3D structure');
          }
        } else {
          // For other formats, use the input directly
          setMolData(moleculeInput);
        }
        
        setIsLoading(false);
      } catch (err) {
        console.error('Error fetching molecule data:', err);
        setError(err.message || 'Failed to load molecule data');
        setIsLoading(false);
      }
    };
    
    fetchMoleculeData();
  }, [moleculeInput, format]);
  
  // Initialize 3Dmol viewer when we have molecule data
  useEffect(() => {
    if (!viewerRef.current || !molData) return;
    
    // Clean up any existing viewer
    if (viewer) {
      viewer.removeAllModels();
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    }
    
    // Initialize viewer with a white background
    const newViewer = $3Dmol.createViewer(viewerRef.current, {
      backgroundColor: 'white'
    });
    
    try {
      // Add the model based on format
      let model;
      if (format === 'pdb') {
        model = newViewer.addModel(molData, 'pdb');
      } else if (format === 'smiles' || format === 'sdf' || format === 'mol') {
        model = newViewer.addModel(molData, 'mol');
      } else {
        model = newViewer.addModel(molData, format);
      }
      
      // Apply style based on initial setting
      if (modelStyle === 'stick') {
        newViewer.setStyle({}, {stick: {}}); 
      } else if (modelStyle === 'sphere') {
        newViewer.setStyle({}, {sphere: {}}); 
      } else if (modelStyle === 'line') {
        newViewer.setStyle({}, {line: {}}); 
      } else if (modelStyle === 'cartoon') {
        newViewer.setStyle({}, {cartoon: {}}); 
      }
      
      // Add a surface with 50% opacity
      newViewer.addSurface($3Dmol.SurfaceType.VDW, {
        opacity: 0.5,
        color: 'white'
      });
      
      // Center the view and render
      newViewer.zoomTo();
      newViewer.render();
      setViewer(newViewer);
      
      // Start a gentle rotation for better visualization
      const animate = () => {
        newViewer.rotate(0.5, {x: 0, y: 1, z: 0});
        newViewer.render();
        animationRef.current = requestAnimationFrame(animate);
      };
      animate();
    } catch (err) {
      console.error('Error rendering molecule:', err);
      setError(`Failed to render molecule: ${err.message}`);
    }
    
    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [molData, format, modelStyle]);
  
  // Handle fullscreen toggle
  const toggleFullscreen = () => {
    setIsFullscreen(!isFullscreen);
    
    // Need to resize the viewer after changing container size
    setTimeout(() => {
      if (viewer) {
        viewer.resize();
        viewer.render();
      }
    }, 100);
  };
  
  // Refresh the viewer
  const handleRefresh = () => {
    if (viewer) {
      viewer.removeAllModels();
      viewer.removeAllSurfaces();
      
      // Re-add the model based on format
      let model;
      if (format === 'pdb') {
        model = viewer.addModel(molData, 'pdb');
      } else if (format === 'smiles' || format === 'sdf' || format === 'mol') {
        model = viewer.addModel(molData, 'mol');
      } else {
        model = viewer.addModel(molData, format);
      }
      
      // Re-apply style
      if (modelStyle === 'stick') {
        viewer.setStyle({}, {stick: {}}); 
      } else if (modelStyle === 'sphere') {
        viewer.setStyle({}, {sphere: {}}); 
      } else if (modelStyle === 'line') {
        viewer.setStyle({}, {line: {}}); 
      } else if (modelStyle === 'cartoon') {
        viewer.setStyle({}, {cartoon: {}}); 
      }
      
      // Re-add surface
      viewer.addSurface($3Dmol.SurfaceType.VDW, {
        opacity: 0.5,
        color: 'white'
      });
      
      viewer.zoomTo();
      viewer.render();
    }
  };
  
  // Style for container based on fullscreen state
  const containerStyle = isFullscreen
    ? classes.fullscreenContainer
    : classes.viewerContainer;
  
  return (
    <div 
      style={{
        height: isFullscreen ? '100%' : height,
        width: isFullscreen ? '100%' : width,
        position: 'relative'
      }}
    >
      <div className={containerStyle} ref={viewerRef}>
        {/* Loading overlay */}
        {isLoading && (
          <div className={classes.loadingOverlay}>
            <CircularProgress />
          </div>
        )}
        
        {/* Error overlay */}
        {error && (
          <div className={classes.errorOverlay}>
            <Typography color="error">
              {error}
            </Typography>
          </div>
        )}
        
        {/* Viewer controls */}
        <div className={classes.viewerControls}>
          <Tooltip title={isFullscreen ? "Exit Fullscreen" : "Fullscreen"}>
            <IconButton size="small" onClick={toggleFullscreen}>
              {isFullscreen ? <FullscreenExitOutlined /> : <FullscreenOutlined />}
            </IconButton>
          </Tooltip>
          <Tooltip title="Refresh Viewer">
            <IconButton size="small" onClick={handleRefresh}>
              <RefreshOutlined />
            </IconButton>
          </Tooltip>
        </div>
      </div>
    </div>
  );
};

MoleculeViewer3DImproved.propTypes = {
  smiles: PropTypes.string,
  moleculeData: PropTypes.string,
  format: PropTypes.oneOf(['smiles', 'pdb', 'mol2', 'sdf', 'cif', 'xyz']),
  initialStyle: PropTypes.oneOf(['stick', 'sphere', 'line', 'cartoon']),
  height: PropTypes.oneOfType([PropTypes.number, PropTypes.string]),
  width: PropTypes.oneOfType([PropTypes.number, PropTypes.string])
};

export default MoleculeViewer3DImproved;
