import React, { useState, useEffect, useRef } from 'react';
import { 
  Paper, 
  Typography, 
  CircularProgress,
  Divider,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Card,
  CardContent,
  Box,
  Chip,
  Button,
  makeStyles 
} from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import BuildIcon from '@material-ui/icons/Build'; // Using Build icon instead of Science
import FormatQuoteIcon from '@material-ui/icons/FormatQuote';
import LocalPharmacyIcon from '@material-ui/icons/LocalPharmacy';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import WarningIcon from '@material-ui/icons/Warning';
import { Alert, AlertTitle } from '@material-ui/lab';
import MoleculeViewer3D from './MoleculeViewer3D';
import { claudeAPI } from '../services/api';

const useStyles = makeStyles((theme) => ({
  root: {
    width: '100%',
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2),
  },
  paper: {
    padding: theme.spacing(3),
  },
  loadingContainer: {
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center',
    padding: theme.spacing(4),
  },
  header: {
    display: 'flex',
    alignItems: 'center',
    marginBottom: theme.spacing(2),
  },
  headerIcon: {
    marginRight: theme.spacing(1),
    color: theme.palette.primary.main,
  },
  divider: {
    margin: theme.spacing(2, 0),
  },
  sectionTitle: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(1),
    fontWeight: 500,
  },
  accordion: {
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1),
  },
  tag: {
    marginRight: theme.spacing(1),
    marginBottom: theme.spacing(1),
  },
  timestamp: {
    marginTop: theme.spacing(1),
    color: theme.palette.text.secondary,
    fontSize: '0.9rem',
  },
  viewerCard: {
    marginBottom: theme.spacing(2),
  },
  moleculeSection: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2),
    padding: theme.spacing(2),
    backgroundColor: theme.palette.background.default,
    borderRadius: theme.shape.borderRadius,
  },
  propertyItem: {
    display: 'flex',
    justifyContent: 'space-between',
    padding: theme.spacing(0.5, 0),
  },
  statusChip: {
    marginLeft: theme.spacing(1),
  },
  codeBlock: {
    fontFamily: 'monospace',
    backgroundColor: theme.palette.background.default,
    padding: theme.spacing(2),
    borderRadius: theme.shape.borderRadius,
    overflowX: 'auto',
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1),
  },
  promptBlock: {
    backgroundColor: '#f0f7ff',
    padding: theme.spacing(2),
    borderRadius: theme.shape.borderRadius,
    marginBottom: theme.spacing(2),
    borderLeft: `4px solid ${theme.palette.primary.main}`,
  },
}));

const extractSections = (thinking) => {
  // This is a simplistic approach; in a real implementation, you'd want
  // more sophisticated parsing based on Claude's actual output structure
  const sections = [];
  
  if (!thinking) return sections;
  
  // Try to identify major sections using headers
  const lines = thinking.split('\n');
  let currentSection = null;
  let currentContent = [];
  
  lines.forEach(line => {
    const trimmedLine = line.trim();
    // Check if this line seems like a header (various formats)
    const isHeader = 
      /^(#+|Candidate|Step|Phase|Molecule \d+:|Analysis of|Design Rationale|Candidate \w+:|Synthesis|Regulatory|Safety|Pharmacokinetics|ADMET|Binding|Manufacturing|References)/i.test(trimmedLine) && 
      trimmedLine.length < 100;
    
    if (isHeader) {
      // Save previous section if it exists
      if (currentSection) {
        sections.push({
          title: currentSection,
          content: currentContent.join('\n')
        });
      }
      
      // Start new section
      currentSection = trimmedLine;
      currentContent = [];
    } else if (currentSection) {
      currentContent.push(line);
    } else {
      // First lines before any clear section header
      if (!currentSection) {
        currentSection = "Introduction";
      }
      currentContent.push(line);
    }
  });
  
  // Add the final section
  if (currentSection) {
    sections.push({
      title: currentSection,
      content: currentContent.join('\n')
    });
  }
  
  return sections;
};

const ThinkingProcess = ({ requestId, onSelectMolecule }) => {
  const classes = useStyles();
  const [thinkingData, setThinkingData] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [expandedSection, setExpandedSection] = useState(null);
  const [sections, setSections] = useState([]);
  
  useEffect(() => {
    if (!requestId) {
      setError('No request ID provided');
      setLoading(false);
      return;
    }
    
    const fetchThinkingProcess = async () => {
      try {
        setLoading(true);
        setError(null);
        
        const response = await claudeAPI.getMoleculeThinking(requestId);
        setThinkingData(response.data);
        
        // Parse the thinking process into sections
        const parsedSections = extractSections(response.data.thinking);
        setSections(parsedSections);
        
        // Expand the first section by default
        if (parsedSections.length > 0) {
          setExpandedSection(0);
        }
        
        setLoading(false);
      } catch (err) {
        console.error('Error fetching thinking process:', err);
        setError(err.response?.data?.error || 'Failed to load thinking process');
        setLoading(false);
      }
    };
    
    fetchThinkingProcess();
  }, [requestId]);
  
  const handleAccordionChange = (index) => (event, isExpanded) => {
    setExpandedSection(isExpanded ? index : null);
  };
  
  // Format timestamp
  const formatDate = (isoString) => {
    if (!isoString) return '';
    const date = new Date(isoString);
    return date.toLocaleString();
  };
  
  // Function to select a molecule when clicked
  const handleSelectMolecule = (smiles) => {
    if (onSelectMolecule) {
      onSelectMolecule(smiles);
    }
  };
  
  // Extract SMILES strings from a section
  const extractSmiles = (content) => {
    const smilesMatches = content.match(/SMILES:?\s*([^\s]+)/g);
    if (!smilesMatches) return [];
    
    return smilesMatches.map(match => {
      return match.replace(/SMILES:?\s*/, '').replace(/[`'"]/g, '');
    });
  };
  
  if (loading) {
    return (
      <Paper className={classes.loadingContainer}>
        <CircularProgress />
        <Typography variant="body1" style={{ marginTop: 16 }}>
          Loading thinking process...
        </Typography>
      </Paper>
    );
  }
  
  if (error) {
    return (
      <Alert severity="error" className={classes.root}>
        <AlertTitle>Error</AlertTitle>
        {error}
      </Alert>
    );
  }
  
  if (!thinkingData) {
    return (
      <Alert severity="warning" className={classes.root}>
        <AlertTitle>No Data</AlertTitle>
        No thinking process data is available.
      </Alert>
    );
  }
  
  return (
    <Paper className={`${classes.root} ${classes.paper}`}>
      <div className={classes.header}>
        <BuildIcon className={classes.headerIcon} fontSize="large" />
        <Typography variant="h5">
          Molecule Design Thinking Process
          <Chip 
            label={thinkingData.status} 
            color={thinkingData.status === 'completed' ? 'primary' : 'default'}
            icon={thinkingData.status === 'completed' ? <CheckCircleIcon /> : <WarningIcon />}
            size="small"
            className={classes.statusChip}
          />
        </Typography>
      </div>
      
      <Typography variant="body2" className={classes.timestamp}>
        Generated: {formatDate(thinkingData.timestamp)}
      </Typography>
      
      <Divider className={classes.divider} />
      
      <Typography variant="h6" className={classes.sectionTitle}>
        Prompt
      </Typography>
      
      <div className={classes.promptBlock}>
        <Typography variant="body1">
          {thinkingData.prompt}
        </Typography>
      </div>
      
      {thinkingData.result && thinkingData.result.smiles && thinkingData.result.smiles.length > 0 && (
        <>
          <Typography variant="h6" className={classes.sectionTitle}>
            Generated Molecules
          </Typography>
          
          <Box display="flex" flexWrap="wrap">
            {thinkingData.result.smiles.map((smiles, index) => (
              <Card key={index} className={classes.viewerCard} style={{ width: 300, margin: 8 }}>
                <CardContent>
                  <Typography variant="subtitle1" gutterBottom>
                    Molecule {index + 1}
                  </Typography>
                  <MoleculeViewer3D 
                    smiles={smiles} 
                    height={200} 
                  />
                  <Typography variant="body2" className={classes.codeBlock}>
                    {smiles}
                  </Typography>
                  <Button 
                    variant="outlined" 
                    color="primary" 
                    fullWidth
                    onClick={() => handleSelectMolecule(smiles)}
                    size="small"
                    startIcon={<LocalPharmacyIcon />}
                  >
                    Select This Molecule
                  </Button>
                </CardContent>
              </Card>
            ))}
          </Box>
        </>
      )}
      
      <Typography variant="h6" className={classes.sectionTitle}>
        Detailed Thinking Process
      </Typography>
      
      {sections.length > 0 ? (
        sections.map((section, index) => {
          const sectionSmiles = extractSmiles(section.content);
          
          return (
            <Accordion 
              key={index} 
              expanded={expandedSection === index} 
              onChange={handleAccordionChange(index)}
              className={classes.accordion}
            >
              <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                <Typography variant="subtitle1">{section.title}</Typography>
              </AccordionSummary>
              <AccordionDetails>
                <Box width="100%">
                  <Typography variant="body1" component="div">
                    {section.content.split('\n').map((line, i) => (
                      <React.Fragment key={i}>
                        {line}
                        <br />
                      </React.Fragment>
                    ))}
                  </Typography>
                  
                  {sectionSmiles.length > 0 && (
                    <div className={classes.moleculeSection}>
                      <Typography variant="subtitle2" gutterBottom>
                        Molecules in this section:
                      </Typography>
                      <Box display="flex" flexWrap="wrap">
                        {sectionSmiles.map((smiles, i) => (
                          <Chip 
                            key={i} 
                            label={`Molecule ${i+1}`} 
                            onClick={() => handleSelectMolecule(smiles)}
                            className={classes.tag}
                            icon={<LocalPharmacyIcon />}
                            color="primary"
                            variant="outlined"
                          />
                        ))}
                      </Box>
                    </div>
                  )}
                </Box>
              </AccordionDetails>
            </Accordion>
          );
        })
      ) : (
        <Typography variant="body1" color="textSecondary">
          No sections could be extracted from the thinking process.
        </Typography>
      )}
      
      <Divider className={classes.divider} />
      
      <Typography variant="h6" className={classes.sectionTitle}>
        Raw Thinking Output
      </Typography>
      
      <Accordion className={classes.accordion}>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography variant="subtitle1">
            Complete Raw Output
          </Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Typography variant="body2" component="div" style={{ whiteSpace: 'pre-wrap', width: '100%' }}>
            {thinkingData.thinking}
          </Typography>
        </AccordionDetails>
      </Accordion>
    </Paper>
  );
};

export default ThinkingProcess; 