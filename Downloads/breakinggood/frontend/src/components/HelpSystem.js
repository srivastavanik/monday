import React, { useState } from 'react';
import {
  Typography,
  Paper,
  Grid,
  Divider,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Collapse,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  IconButton,
  makeStyles,
  Tooltip,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Card,
  CardContent,
  CardHeader,
  TextField,
  Snackbar
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import SearchIcon from '@material-ui/icons/Search';
import HelpIcon from '@material-ui/icons/Help';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import InfoIcon from '@material-ui/icons/Info';
import ScienceIcon from '@material-ui/icons/Science';
import MenuBookIcon from '@material-ui/icons/MenuBook';
import MemoryIcon from '@material-ui/icons/Memory';
import ShowChartIcon from '@material-ui/icons/ShowChart';
import BuildIcon from '@material-ui/icons/Build';
import BookmarkIcon from '@material-ui/icons/Bookmark';
import SearchOutlinedIcon from '@material-ui/icons/SearchOutlined';
import HelpOutlineIcon from '@material-ui/icons/HelpOutline';
import ImportContactsIcon from '@material-ui/icons/ImportContacts';
import VideocamIcon from '@material-ui/icons/Videocam';
import ForumIcon from '@material-ui/icons/Forum';

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
  divider: {
    margin: theme.spacing(2, 0),
  },
  searchContainer: {
    marginBottom: theme.spacing(3),
  },
  helpCategoryIcon: {
    color: theme.palette.primary.main,
    fontSize: '2rem',
  },
  helpCategory: {
    marginBottom: theme.spacing(2),
  },
  accordion: {
    marginBottom: theme.spacing(1),
  },
  accordionSummary: {
    fontWeight: 500,
  },
  accordionDetails: {
    display: 'block',
  },
  videoCard: {
    marginBottom: theme.spacing(2),
  },
  videoPlaceholder: {
    backgroundColor: '#f5f5f5',
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center',
    height: 180,
    marginBottom: theme.spacing(1),
  },
  infoCard: {
    marginBottom: theme.spacing(2),
  },
  tooltip: {
    maxWidth: 320,
    fontSize: '0.875rem',
  },
  helpQuestion: {
    borderLeft: `3px solid ${theme.palette.primary.main}`,
    paddingLeft: theme.spacing(2),
    marginBottom: theme.spacing(2),
  },
}));

// Tooltip Component that displays detailed help information
export const HelpTooltip = ({ title, children }) => {
  const classes = useStyles();
  
  return (
    <Tooltip
      title={title}
      arrow
      classes={{ tooltip: classes.tooltip }}
      enterTouchDelay={50}
      leaveTouchDelay={3000}
    >
      {children || <HelpOutlineIcon color="action" fontSize="small" />}
    </Tooltip>
  );
};

// Context-sensitive help dialog that provides detailed information
export const ContextHelpDialog = ({ open, onClose, title, content }) => {
  return (
    <Dialog 
      open={open} 
      onClose={onClose}
      maxWidth="md"
      fullWidth
    >
      <DialogTitle>{title || 'Help Information'}</DialogTitle>
      <DialogContent dividers>
        {content}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose} color="primary">
          Close
        </Button>
      </DialogActions>
    </Dialog>
  );
};

// Main Help System component
const HelpSystem = () => {
  const classes = useStyles();
  const [searchQuery, setSearchQuery] = useState('');
  const [searchResults, setSearchResults] = useState(null);
  const [expanded, setExpanded] = useState(false);
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [snackbarMessage, setSnackbarMessage] = useState('');
  
  const handleSearchChange = (event) => {
    setSearchQuery(event.target.value);
  };
  
  const handleSearch = () => {
    if (!searchQuery.trim()) {
      setSearchResults(null);
      return;
    }
    
    // Simulated search functionality
    const query = searchQuery.toLowerCase();
    const results = helpContent.filter(category => {
      // Search in category title
      if (category.title.toLowerCase().includes(query)) return true;
      
      // Search in category FAQs
      if (category.faqs && category.faqs.some(faq => 
        faq.question.toLowerCase().includes(query) || 
        faq.answer.toLowerCase().includes(query)
      )) return true;
      
      // Search in tutorials
      if (category.tutorials && category.tutorials.some(tutorial => 
        tutorial.title.toLowerCase().includes(query) || 
        tutorial.description.toLowerCase().includes(query)
      )) return true;
      
      return false;
    });
    
    setSearchResults(results);
    
    if (results.length === 0) {
      setSnackbarMessage('No results found. Try different keywords.');
      setSnackbarOpen(true);
    }
  };
  
  const handleAccordionChange = (panel) => (event, isExpanded) => {
    setExpanded(isExpanded ? panel : false);
  };
  
  const handleSnackbarClose = () => {
    setSnackbarOpen(false);
  };
  
  const helpContent = [
    {
      id: 'general',
      title: 'Getting Started',
      icon: <InfoIcon className={classes.helpCategoryIcon} />,
      faqs: [
        {
          question: 'What is Breaking Good?',
          answer: 'Breaking Good is a comprehensive drug discovery platform that integrates AI-powered molecule design, simulation tools, and database access to streamline the discovery of novel therapeutic compounds. The platform is designed for researchers in academia and industry to collaborate on designing new pharmaceutical compounds with emphasis on ADHD and related disorders.'
        },
        {
          question: 'How do I navigate the interface?',
          answer: 'The Breaking Good interface is organized into several key areas: the Molecule Designer for creating new compounds, Simulations for testing interaction with biological targets, the ChEMBL and PubMed integrations for literature research, and the Data Visualization tools for analyzing results. Use the sidebar navigation to access different features of the platform.'
        },
        {
          question: 'Where is my data stored?',
          answer: 'In the current version, data is stored locally in your browser using localStorage. This means your molecules and simulations persist between browser sessions but are limited to your current device. Future updates will include cloud storage options for cross-device access and team collaboration.'
        }
      ],
      tutorials: [
        {
          title: 'Platform Overview',
          description: 'A comprehensive introduction to the Breaking Good platform and its key features.',
          videoUrl: '#'
        },
        {
          title: 'Your First Molecule',
          description: 'A step-by-step guide to designing your first molecule using the Molecule Designer.',
          videoUrl: '#'
        }
      ]
    },
    {
      id: 'molecule-designer',
      title: 'Molecule Designer',
      icon: <ScienceIcon className={classes.helpCategoryIcon} />,
      faqs: [
        {
          question: 'How do I create a new molecule?',
          answer: 'To create a new molecule, navigate to the Molecule Designer tab. You can start from scratch by using the structure editor, or provide a SMILES string of an existing molecule. The AI generator can also suggest modifications based on your target receptor selection. Once you have designed your molecule, you can save it to your local database.'
        },
        {
          question: 'What does the "Generate Molecules" function do?',
          answer: 'The "Generate Molecules" function uses AI to create novel molecules based on your inputs. You can specify a target receptor and/or a starting structure, and the system will generate candidate molecules optimized for binding to that target. The generator considers drug-like properties, synthetic feasibility, and potential activity against the specified target.'
        },
        {
          question: 'Can I import molecules from external sources?',
          answer: 'Yes! You can import molecules in several formats: SMILES strings, MOL files, or structured JSON data. Use the Import function in the Molecule Storage component to upload your molecules. You can also fetch molecules directly from the ChEMBL database using the integrated search feature.'
        }
      ],
      tutorials: [
        {
          title: 'Advanced Molecule Design',
          description: 'Learn advanced techniques for molecule design including fragment-based approaches and scaffold hopping.',
          videoUrl: '#'
        },
        {
          title: 'Working with SMILES Notation',
          description: 'Understanding SMILES notation and how to use it effectively in molecule design.',
          videoUrl: '#'
        }
      ]
    },
    {
      id: 'simulations',
      title: 'Simulations',
      icon: <MemoryIcon className={classes.helpCategoryIcon} />,
      faqs: [
        {
          question: 'What types of simulations are available?',
          answer: 'Breaking Good offers three main types of simulations: Molecular Docking (predicting binding pose and affinity), Molecular Dynamics (simulating molecular motion over time), and ADMET Prediction (estimating absorption, distribution, metabolism, excretion, and toxicity properties).'
        },
        {
          question: 'How do I interpret simulation results?',
          answer: 'Simulation results include both numerical data and visualizations. For molecular docking, look at the binding affinity (lower values indicate stronger binding) and interaction sites. For dynamics simulations, the stability percentage indicates how consistently the molecule maintains its binding pose. ADMET predictions provide estimated values for key pharmacokinetic properties with confidence intervals.'
        },
        {
          question: 'Can I compare simulation results between molecules?',
          answer: 'Yes, the platform includes a comparison tool that allows you to select multiple molecules and compare their simulation results side by side. This is particularly useful for evaluating different modifications of a base compound to identify the most promising candidates.'
        }
      ],
      tutorials: [
        {
          title: 'Running Your First Simulation',
          description: 'A guided walkthrough of setting up and running your first molecular docking simulation.',
          videoUrl: '#'
        },
        {
          title: 'Advanced Simulation Parameters',
          description: 'Understanding the advanced parameters for fine-tuning simulations for specific use cases.',
          videoUrl: '#'
        }
      ]
    },
    {
      id: 'database',
      title: 'Database Access',
      icon: <MenuBookIcon className={classes.helpCategoryIcon} />,
      faqs: [
        {
          question: 'How do I search the ChEMBL database?',
          answer: 'Navigate to the ChEMBL tab and use the search interface to query the database. You can search by molecule name, SMILES, InChI, or other identifiers. The search supports various filter types like exact match, contains, starts with, etc. Results include molecule structures, properties, and biological activity data when available.'
        },
        {
          question: 'How do I access PubMed articles?',
          answer: 'Use the PubMed/BioC tab to access scientific articles. You can search by PubMed ID (PMID) or PubMed Central ID (PMCID). The platform retrieves article content in BioC format, which preserves the structure of the document including title, abstract, body sections, and references. Note that only open access articles are available for full-text retrieval.'
        },
        {
          question: 'Can I save references for later use?',
          answer: 'Currently, you can export article information as JSON or CSV files. In future updates, we will add a reference management system that allows you to save articles, link them to molecules, and organize them into collections for specific research projects.'
        }
      ],
      tutorials: [
        {
          title: 'Effective Database Searching',
          description: 'Tips and techniques for efficient searching across chemical and literature databases.',
          videoUrl: '#'
        },
        {
          title: 'Integrating Literature with Molecule Design',
          description: 'How to use literature findings to inform your molecule design process.',
          videoUrl: '#'
        }
      ]
    },
    {
      id: 'visualization',
      title: 'Visualization Tools',
      icon: <ShowChartIcon className={classes.helpCategoryIcon} />,
      faqs: [
        {
          question: 'What visualization options are available?',
          answer: 'The platform offers both 2D and 3D molecular visualizations. In 3D mode, you can choose between ball-and-stick, stick, and surface representations. You can customize colors, show/hide hydrogens and labels, and adjust rotation speed. Visualizations can be exported as image files for use in reports or publications.'
        },
        {
          question: 'How do I view protein-ligand interactions?',
          answer: 'After running a docking simulation, you can view the protein-ligand complex in the 3D visualization tool. Key interaction sites are highlighted, and you can toggle visibility of different molecular components. The visualization tool also supports measurement of distances between atoms to analyze hydrogen bonds and other interactions.'
        },
        {
          question: 'Can I create custom visualizations of my data?',
          answer: 'Yes, the platform includes a data visualization module that allows you to create custom charts and graphs from your simulation results. You can plot properties across multiple molecules, create correlation plots, and generate heatmaps of binding affinities across different targets.'
        }
      ],
      tutorials: [
        {
          title: 'Advanced 3D Visualization',
          description: 'Mastering the 3D visualization tools for analyzing molecular structures and interactions.',
          videoUrl: '#'
        },
        {
          title: 'Creating Publication-Quality Images',
          description: 'How to create and export high-quality molecular visualizations for publications.',
          videoUrl: '#'
        }
      ]
    },
    {
      id: 'best-practices',
      title: 'Best Practices',
      icon: <BuildIcon className={classes.helpCategoryIcon} />,
      faqs: [
        {
          question: 'What are some tips for effective molecule design?',
          answer: 'Focus on small, incremental changes to molecules with known activity. Consider Lipinski\'s Rule of Five for drug-like properties. Use the AI suggestions but review them critically. Keep track of your design rationale for each modification. Regularly use the comparison tools to identify the most promising candidates.'
        },
        {
          question: 'How should I organize my workflow?',
          answer: 'Start with literature research to identify promising scaffolds or active compounds. Create a series of structural modifications using the Molecule Designer. Run docking simulations to predict binding affinity. Further assess promising candidates with ADMET predictions. Save and export your most promising compounds for experimental testing.'
        },
        {
          question: 'How can I ensure reproducibility of my results?',
          answer: 'Keep detailed records of all parameters used in simulations. Export your molecules and results regularly. Document your design rationale and hypotheses. Use version control by creating named snapshots of important molecules. When sharing results, include all relevant parameters and the exact SMILES strings of the molecules studied.'
        }
      ]
    }
  ];
  
  const renderSearchResults = () => {
    if (!searchResults) return null;
    
    if (searchResults.length === 0) {
      return (
        <Paper className={classes.paper}>
          <Typography variant="body1">
            No results found for "{searchQuery}". Please try different keywords.
          </Typography>
        </Paper>
      );
    }
    
    return (
      <div>
        <Typography variant="h6" gutterBottom>
          Search Results for "{searchQuery}"
        </Typography>
        
        {searchResults.map((category) => (
          <div key={category.id} className={classes.helpCategory}>
            <Typography variant="subtitle1" gutterBottom>
              {category.icon} {category.title}
            </Typography>
            
            {category.faqs && category.faqs.some(faq => 
              faq.question.toLowerCase().includes(searchQuery.toLowerCase()) ||
              faq.answer.toLowerCase().includes(searchQuery.toLowerCase())
            ) && (
              <div>
                <Typography variant="subtitle2" gutterBottom>
                  FAQs:
                </Typography>
                
                {category.faqs.filter(faq => 
                  faq.question.toLowerCase().includes(searchQuery.toLowerCase()) ||
                  faq.answer.toLowerCase().includes(searchQuery.toLowerCase())
                ).map((faq, index) => (
                  <div key={index} className={classes.helpQuestion}>
                    <Typography variant="body1" gutterBottom>
                      <strong>Q: {faq.question}</strong>
                    </Typography>
                    <Typography variant="body2">
                      {faq.answer}
                    </Typography>
                  </div>
                ))}
              </div>
            )}
            
            {category.tutorials && category.tutorials.some(tutorial => 
              tutorial.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
              tutorial.description.toLowerCase().includes(searchQuery.toLowerCase())
            ) && (
              <div>
                <Typography variant="subtitle2" gutterBottom>
                  Tutorials:
                </Typography>
                
                {category.tutorials.filter(tutorial => 
                  tutorial.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
                  tutorial.description.toLowerCase().includes(searchQuery.toLowerCase())
                ).map((tutorial, index) => (
                  <Card key={index} className={classes.videoCard}>
                    <CardContent>
                      <Typography variant="subtitle1">
                        {tutorial.title}
                      </Typography>
                      <Typography variant="body2" color="textSecondary">
                        {tutorial.description}
                      </Typography>
                    </CardContent>
                  </Card>
                ))}
              </div>
            )}
            
            <Divider className={classes.divider} />
          </div>
        ))}
      </div>
    );
  };
  
  return (
    <div className={classes.root}>
      <Typography variant="h5" className={classes.title}>
        Help & Documentation
      </Typography>
      
      <Paper className={classes.searchContainer}>
        <Grid container spacing={2} alignItems="center" style={{ padding: 16 }}>
          <Grid item xs>
            <TextField
              fullWidth
              variant="outlined"
              placeholder="Search for help topics..."
              value={searchQuery}
              onChange={handleSearchChange}
              InputProps={{
                startAdornment: <SearchOutlinedIcon color="action" style={{ marginRight: 8 }} />,
              }}
              onKeyPress={(e) => {
                if (e.key === 'Enter') {
                  handleSearch();
                }
              }}
            />
          </Grid>
          <Grid item>
            <Button
              variant="contained"
              color="primary"
              startIcon={<SearchIcon />}
              onClick={handleSearch}
            >
              Search
            </Button>
          </Grid>
        </Grid>
      </Paper>
      
      {searchResults ? (
        renderSearchResults()
      ) : (
        <Grid container spacing={3}>
          <Grid item xs={12} md={8}>
            <Paper className={classes.paper}>
              <Typography variant="h6" gutterBottom>
                Frequently Asked Questions
              </Typography>
              
              {helpContent.map((category) => (
                <Accordion 
                  key={category.id}
                  expanded={expanded === category.id}
                  onChange={handleAccordionChange(category.id)}
                  className={classes.accordion}
                >
                  <AccordionSummary
                    expandIcon={<ExpandMoreIcon />}
                    aria-controls={`${category.id}-content`}
                    id={`${category.id}-header`}
                  >
                    <Grid container spacing={1} alignItems="center">
                      <Grid item>
                        {category.icon}
                      </Grid>
                      <Grid item>
                        <Typography className={classes.accordionSummary}>
                          {category.title}
                        </Typography>
                      </Grid>
                    </Grid>
                  </AccordionSummary>
                  <AccordionDetails className={classes.accordionDetails}>
                    {category.faqs && category.faqs.map((faq, index) => (
                      <div key={index} className={classes.helpQuestion}>
                        <Typography variant="body1" gutterBottom>
                          <strong>Q: {faq.question}</strong>
                        </Typography>
                        <Typography variant="body2">
                          {faq.answer}
                        </Typography>
                      </div>
                    ))}
                  </AccordionDetails>
                </Accordion>
              ))}
            </Paper>
          </Grid>
          
          <Grid item xs={12} md={4}>
            <Paper className={classes.paper}>
              <Typography variant="h6" gutterBottom>
                Video Tutorials
              </Typography>
              
              {helpContent.flatMap(category => 
                category.tutorials || []
              ).slice(0, 3).map((tutorial, index) => (
                <Card key={index} className={classes.videoCard}>
                  <div className={classes.videoPlaceholder}>
                    <VideocamIcon style={{ fontSize: 48, color: '#757575', marginBottom: 8 }} />
                    <Typography variant="body2" color="textSecondary">
                      Video Tutorial
                    </Typography>
                  </div>
                  <CardContent>
                    <Typography variant="subtitle1">
                      {tutorial.title}
                    </Typography>
                    <Typography variant="body2" color="textSecondary">
                      {tutorial.description}
                    </Typography>
                  </CardContent>
                </Card>
              ))}
              
              <Button
                variant="outlined"
                color="primary"
                startIcon={<ImportContactsIcon />}
                fullWidth
                style={{ marginTop: 16 }}
              >
                View All Tutorials
              </Button>
            </Paper>
            
            <Paper className={classes.paper} style={{ marginTop: 24 }}>
              <Typography variant="h6" gutterBottom>
                Additional Resources
              </Typography>
              
              <List>
                <ListItem button>
                  <ListItemIcon>
                    <ImportContactsIcon />
                  </ListItemIcon>
                  <ListItemText 
                    primary="User Guide" 
                    secondary="Complete documentation of all features"
                  />
                </ListItem>
                <ListItem button>
                  <ListItemIcon>
                    <BookmarkIcon />
                  </ListItemIcon>
                  <ListItemText 
                    primary="API Reference" 
                    secondary="For developers and integrations"
                  />
                </ListItem>
                <ListItem button>
                  <ListItemIcon>
                    <ForumIcon />
                  </ListItemIcon>
                  <ListItemText 
                    primary="Community Forum" 
                    secondary="Ask questions and share insights"
                  />
                </ListItem>
              </List>
            </Paper>
          </Grid>
        </Grid>
      )}
      
      <Snackbar
        open={snackbarOpen}
        autoHideDuration={6000}
        onClose={handleSnackbarClose}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
      >
        <Alert onClose={handleSnackbarClose} severity="info">
          {snackbarMessage}
        </Alert>
      </Snackbar>
    </div>
  );
};

export default HelpSystem; 