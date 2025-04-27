import React, { useState } from 'react';
import { 
  Typography, 
  Grid, 
  Paper, 
  TextField, 
  Button, 
  FormControl, 
  InputLabel, 
  Select, 
  MenuItem, 
  Divider,
  CircularProgress,
  FormControlLabel,
  Switch,
  Chip,
  Tabs,
  Tab,
  Box,
  makeStyles 
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import TimelineIcon from '@material-ui/icons/Timeline';
import AssignmentIcon from '@material-ui/icons/Assignment';
import BuildIcon from '@material-ui/icons/Build';
import TrendingUpIcon from '@material-ui/icons/TrendingUp';
import { regulatoryAPI } from '../services/api';
import ReactMarkdown from 'react-markdown';

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
    margin: theme.spacing(1, 0),
    minWidth: '100%',
  },
  sectionTitle: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(2),
  },
  divider: {
    margin: theme.spacing(3, 0),
  },
  progress: {
    display: 'flex',
    justifyContent: 'center',
    margin: theme.spacing(4, 0),
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
  switchFormControl: {
    margin: theme.spacing(1, 0),
  },
  alertMargin: {
    marginBottom: theme.spacing(2),
  },
  reportContainer: {
    marginTop: theme.spacing(2),
    padding: theme.spacing(2),
    backgroundColor: theme.palette.grey[50],
    borderRadius: theme.shape.borderRadius,
    border: `1px solid ${theme.palette.divider}`,
    maxHeight: '70vh',
    overflowY: 'auto',
  },
  markdownContent: {
    '& h1, & h2, & h3, & h4': {
        marginTop: theme.spacing(2.5),
        marginBottom: theme.spacing(1),
        fontWeight: 500,
    },
    '& p': {
        marginBottom: theme.spacing(1.5),
        lineHeight: 1.6,
    },
    '& ul, & ol': {
        marginBottom: theme.spacing(1.5),
        paddingLeft: theme.spacing(3),
    },
    '& li': {
        marginBottom: theme.spacing(0.5),
    },
    '& pre': {
      backgroundColor: theme.palette.grey[200],
      padding: theme.spacing(1),
      borderRadius: theme.shape.borderRadius,
      overflowX: 'auto',
      whiteSpace: 'pre-wrap',
      wordBreak: 'break-all',
    },
    '& code': {
      fontFamily: 'monospace',
      backgroundColor: theme.palette.grey[100],
      padding: theme.spacing(0.2, 0.5),
      borderRadius: 3,
    },
  },
}));

const RegulatoryAnalysis = () => {
  const classes = useStyles();
  const [loading, setLoading] = useState(false);
  const [parameters, setParameters] = useState({
    smiles: '',
    drugClass: 'CNS stimulant',
    novelMechanism: false,
    orphanDrug: false,
    fastTrack: false,
    targetIndication: 'ADHD',
    primaryMechanism: 'dopamine/norepinephrine reuptake inhibition',
  });
  const [analysisReportText, setAnalysisReportText] = useState(null);
  const [error, setError] = useState(null);
  
  const molecules = [
    { id: 1, name: 'Methylphenidate', smiles: 'CN(C)C(C1=CC=CC=C1)C(C)OC(=O)C' },
    { id: 2, name: 'Amphetamine', smiles: 'CC(N)CC1=CC=CC=C1' },
    { id: 3, name: 'Atomoxetine', smiles: 'CC(C)NCC1=CC=CC(OC2=CC=CC=C2)=C1' },
    { id: 4, name: 'Novel Amphetamine Derivative', smiles: 'CC(CC1=CC=C(C=C1)O)NC' }
  ];
  
  const drugClasses = [
    'CNS stimulant',
    'Non-stimulant',
    'Novel mechanism',
    'Prodrug'
  ];
  
  const targetIndications = [
    'ADHD',
    'ADHD with comorbid anxiety',
    'Adult ADHD',
    'Pediatric ADHD'
  ];
  
  const mechanisms = [
    'dopamine/norepinephrine reuptake inhibition',
    'norepinephrine reuptake inhibition',
    'dopamine/norepinephrine release',
    'alpha-2A adrenergic receptor agonism',
    'novel dual mechanism'
  ];
  
  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setParameters({
      ...parameters,
      [name]: value
    });
  };
  
  const handleSwitchChange = (e) => {
    const { name, checked } = e.target;
    setParameters({
      ...parameters,
      [name]: checked
    });
  };
  
  const handleReset = () => {
    setParameters({
      smiles: '',
      drugClass: 'CNS stimulant',
      novelMechanism: false,
      orphanDrug: false,
      fastTrack: false,
      targetIndication: 'ADHD',
      primaryMechanism: 'dopamine/norepinephrine reuptake inhibition',
    });
    setAnalysisReportText(null);
    setError(null);
  };
  
  const handleRunAnalysis = async () => {
    if (!parameters.smiles) {
      setError('Please select a molecule (SMILES required)');
      return;
    }
    
    setLoading(true);
    setError(null);
    setAnalysisReportText(null);
    
    try {
      const response = await regulatoryAPI.generateRegulatoryReport(parameters); 
      
      if (response.data && response.data.regulatoryReport) {
        setAnalysisReportText(response.data.regulatoryReport);
      } else {
        throw new Error('Invalid response from server or missing report data');
      }

    } catch (err) {
      console.error('Error running regulatory analysis:', err);
      setError(err.response?.data?.error || 'An error occurred during analysis. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  const renderMarkdown = (text) => {
      if (!text) return null;
      if (text.startsWith('Error:')){
          return <Alert severity="warning">{text}</Alert>;
      }
      return <ReactMarkdown className={classes.markdownContent}>{text}</ReactMarkdown>;
  };

  return (
    <div className={classes.root}>
      <Typography variant="h4" className={classes.title}>
        Regulatory Analysis Report Generator
      </Typography>
      
      <Grid container spacing={3}>
        <Grid item xs={12} md={4}>
          <Paper className={classes.paper}>
            <Typography variant="h6" gutterBottom>
              Analysis Parameters
            </Typography>
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel id="molecule-select-label">Select Molecule</InputLabel>
              <Select
                labelId="molecule-select-label"
                name="smiles"
                value={parameters.smiles}
                onChange={handleInputChange}
                label="Select Molecule"
              >
                <MenuItem value="">
                  <em>Select a molecule</em>
                </MenuItem>
                {molecules.map((molecule) => (
                  <MenuItem key={molecule.id} value={molecule.smiles}>
                    {molecule.name} ({molecule.smiles.substring(0, 20)}...)
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel>Drug Class</InputLabel>
              <Select name="drugClass" value={parameters.drugClass} onChange={handleInputChange} label="Drug Class">
                {drugClasses.map((cls) => <MenuItem key={cls} value={cls}>{cls}</MenuItem>)}
              </Select>
            </FormControl>
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel>Target Indication</InputLabel>
              <Select name="targetIndication" value={parameters.targetIndication} onChange={handleInputChange} label="Target Indication">
                {targetIndications.map((ind) => <MenuItem key={ind} value={ind}>{ind}</MenuItem>)}
              </Select>
            </FormControl>
            
            <FormControl variant="outlined" className={classes.formControl}>
              <InputLabel>Primary Mechanism</InputLabel>
              <Select name="primaryMechanism" value={parameters.primaryMechanism} onChange={handleInputChange} label="Primary Mechanism">
                {mechanisms.map((mech) => <MenuItem key={mech} value={mech}>{mech}</MenuItem>)}
              </Select>
            </FormControl>
            
            <Typography variant="subtitle1" className={classes.sectionTitle}>
              Special Designations
            </Typography>
            
            <FormControlLabel control={<Switch checked={parameters.novelMechanism} onChange={handleSwitchChange} name="novelMechanism" color="primary" />} label="Novel Mechanism" className={classes.switchFormControl}/>
            
            <FormControlLabel control={<Switch checked={parameters.orphanDrug} onChange={handleSwitchChange} name="orphanDrug" color="primary" />} label="Orphan Drug" className={classes.switchFormControl} />
            
            <FormControlLabel control={<Switch checked={parameters.fastTrack} onChange={handleSwitchChange} name="fastTrack" color="primary" />} label="Fast Track" className={classes.switchFormControl} />
            
            {error && (
              <Alert severity="error" className={classes.alertMargin}>{error}</Alert>
            )}
            
            <div className={classes.buttonWrapper}>
              <Button variant="contained" color="primary" fullWidth onClick={handleRunAnalysis} disabled={loading || !parameters.smiles}>
                Generate Report
              </Button>
              {loading && <CircularProgress size={24} className={classes.buttonProgress} />}
            </div>
            
            {analysisReportText && (
              <Button variant="outlined" color="secondary" fullWidth onClick={handleReset} style={{ marginTop: 16 }}>
                Reset / New Analysis
              </Button>
            )}
          </Paper>
        </Grid>
        
        <Grid item xs={12} md={8}>
          {loading ? (
             <Paper className={classes.paper} style={{ display: 'flex', alignItems: 'center', justifyContent: 'center'}}>
                <CircularProgress />
             </Paper>
          ) : analysisReportText ? (
            <Paper className={`${classes.paper} ${classes.reportContainer}`}>
              {renderMarkdown(analysisReportText)}
            </Paper>
          ) : (
            <Paper className={classes.paper} style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', minHeight: 500, flexDirection: 'column' }}>
              <AssignmentIcon className={classes.emptyStateIcon} />
              <Typography variant="h6" color="textSecondary">
                No regulatory analysis generated
              </Typography>
              <Typography variant="body2" color="textSecondary">
                Select a molecule and parameters, then click "Generate Report"
              </Typography>
            </Paper>
          )}
        </Grid>
      </Grid>
    </div>
  );
};

export default RegulatoryAnalysis; 