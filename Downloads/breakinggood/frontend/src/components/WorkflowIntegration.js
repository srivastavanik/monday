import React, { useState, useEffect } from 'react';
import {
  Typography,
  Paper,
  Grid,
  Card,
  CardContent,
  CardHeader,
  CardActions,
  Button,
  Stepper,
  Step,
  StepLabel,
  StepContent,
  Divider,
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Chip,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  CircularProgress,
  makeStyles,
  IconButton
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';
import SearchIcon from '@material-ui/icons/Search';
import ArrowForwardIcon from '@material-ui/icons/ArrowForward';
import ScienceIcon from '@material-ui/icons/Science';
import MemoryIcon from '@material-ui/icons/Memory';
import MenuBookIcon from '@material-ui/icons/MenuBook';
import VisibilityIcon from '@material-ui/icons/Visibility';
import SaveIcon from '@material-ui/icons/Save';
import PlayArrowIcon from '@material-ui/icons/PlayArrow';
import EditIcon from '@material-ui/icons/Edit';
import AddIcon from '@material-ui/icons/Add';
import SettingsIcon from '@material-ui/icons/Settings';

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
  card: {
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    position: 'relative',
  },
  cardContent: {
    flexGrow: 1,
  },
  workflowStepper: {
    backgroundColor: 'transparent',
  },
  workflowStep: {
    cursor: 'pointer',
    '&:hover': {
      backgroundColor: theme.palette.action.hover,
    },
  },
  workflowStepActive: {
    backgroundColor: theme.palette.action.selected,
  },
  activeTaskIndicator: {
    position: 'absolute',
    top: 0,
    left: 0,
    width: '100%',
    height: '4px',
    backgroundColor: theme.palette.primary.main,
  },
  chip: {
    margin: theme.spacing(0.5),
  },
  taskIcon: {
    fontSize: '3rem',
    marginBottom: theme.spacing(1),
    color: theme.palette.primary.main,
  },
  taskIconSecondary: {
    fontSize: '3rem',
    marginBottom: theme.spacing(1),
    color: theme.palette.secondary.main,
  },
  arrowIcon: {
    fontSize: '2rem',
    color: theme.palette.text.secondary,
    alignSelf: 'center',
  },
  workflowInfo: {
    marginBottom: theme.spacing(2),
  },
  workflowCard: {
    marginBottom: theme.spacing(2),
  },
  taskList: {
    marginTop: theme.spacing(2),
  },
  taskStatusChip: {
    marginLeft: theme.spacing(1),
  },
  button: {
    margin: theme.spacing(1),
  },
  closeButton: {
    position: 'absolute',
    right: theme.spacing(1),
    top: theme.spacing(1),
  },
  summaryItem: {
    display: 'flex',
    justifyContent: 'space-between',
    padding: theme.spacing(1, 0),
    borderBottom: `1px solid ${theme.palette.divider}`,
  },
  summaryLabel: {
    fontWeight: 500,
  },
}));

// This component provides workflow integration between different parts of the application
const WorkflowIntegration = ({ onLaunchComponent, currentUser }) => {
  const classes = useStyles();
  const [activeWorkflow, setActiveWorkflow] = useState(null);
  const [workflowDialogOpen, setWorkflowDialogOpen] = useState(false);
  const [newWorkflowName, setNewWorkflowName] = useState('');
  const [newWorkflowDescription, setNewWorkflowDescription] = useState('');
  const [workflows, setWorkflows] = useState([]);
  const [loading, setLoading] = useState(true);
  const [activeStep, setActiveStep] = useState(0);
  const [error, setError] = useState(null);
  const [taskDialogOpen, setTaskDialogOpen] = useState(false);
  const [selectedTask, setSelectedTask] = useState(null);

  // Load workflows from localStorage
  useEffect(() => {
    const loadWorkflows = () => {
      try {
        const storedWorkflows = localStorage.getItem('workflows');
        
        setTimeout(() => {
          if (storedWorkflows) {
            setWorkflows(JSON.parse(storedWorkflows));
          } else {
            // Create example workflow if none exists
            const exampleWorkflows = [
              {
                id: 'workflow-1',
                name: 'ADHD Drug Discovery',
                description: 'Workflow for discovering novel ADHD treatments targeting the dopamine transporter',
                created: '2023-05-15',
                lastModified: '2023-06-01',
                status: 'in_progress',
                currentStep: 2,
                steps: [
                  {
                    id: 'step-1',
                    name: 'Literature Review',
                    description: 'Research existing ADHD treatments and their mechanisms',
                    component: 'PubMedBioC',
                    status: 'completed',
                    results: {
                      articles: 5,
                      keyFindings: [
                        'Methylphenidate blocks dopamine reuptake',
                        'Amphetamines increase dopamine release',
                        'Non-stimulants affect norepinephrine pathways'
                      ]
                    }
                  },
                  {
                    id: 'step-2',
                    name: 'Structural Analysis',
                    description: 'Analyze structures of known active compounds',
                    component: 'ChEMBL',
                    status: 'completed',
                    results: {
                      molecules: 3,
                      structures: ['Methylphenidate', 'Amphetamine', 'Atomoxetine']
                    }
                  },
                  {
                    id: 'step-3',
                    name: 'Molecule Design',
                    description: 'Design novel compounds based on known actives',
                    component: 'MoleculeDesigner',
                    status: 'in_progress',
                    results: {
                      designedMolecules: 2,
                      pending: true
                    }
                  },
                  {
                    id: 'step-4',
                    name: 'Docking Simulation',
                    description: 'Simulate binding to dopamine transporter',
                    component: 'Simulations',
                    status: 'not_started',
                    results: null
                  },
                  {
                    id: 'step-5',
                    name: 'Results Analysis',
                    description: 'Analyze and visualize simulation results',
                    component: 'Visualization',
                    status: 'not_started',
                    results: null
                  }
                ],
                molecules: [
                  {
                    id: 'mol-1',
                    name: 'Modified Methylphenidate',
                    smiles: 'CN(C)C(C1=CC=CC=C1)C(C)OC(=O)C',
                    status: 'design_complete'
                  },
                  {
                    id: 'mol-2',
                    name: 'Novel DAT Inhibitor',
                    smiles: 'CC1=CC(=C(C=C1)NC(=O)C)NC2=CC=CC=C2',
                    status: 'design_complete'
                  }
                ]
              },
              {
                id: 'workflow-2',
                name: 'Novel Stimulant Development',
                description: 'Research and design of alternative stimulant compounds with reduced side effects',
                created: '2023-04-10',
                lastModified: '2023-04-25',
                status: 'paused',
                currentStep: 1,
                steps: [
                  {
                    id: 'step-1',
                    name: 'Target Selection',
                    description: 'Identify and prioritize molecular targets',
                    component: 'PubMedBioC',
                    status: 'completed',
                    results: {
                      targets: ['Dopamine Transporter', 'Dopamine Receptors', 'Norepinephrine Transporters']
                    }
                  },
                  {
                    id: 'step-2',
                    name: 'Compound Screening',
                    description: 'Screen existing compounds for activity',
                    component: 'ChEMBL',
                    status: 'in_progress',
                    results: {
                      screened: 25,
                      hits: 3
                    }
                  }
                ],
                molecules: []
              }
            ];
            
            setWorkflows(exampleWorkflows);
            localStorage.setItem('workflows', JSON.stringify(exampleWorkflows));
          }
          
          setLoading(false);
        }, 1000);
      } catch (error) {
        console.error('Error loading workflows:', error);
        setError('Error loading workflows. Please try again.');
        setLoading(false);
      }
    };
    
    loadWorkflows();
  }, []);
  
  // Set active workflow on initial load
  useEffect(() => {
    if (workflows.length > 0 && !activeWorkflow) {
      const inProgressWorkflows = workflows.filter(w => w.status === 'in_progress');
      if (inProgressWorkflows.length > 0) {
        setActiveWorkflow(inProgressWorkflows[0]);
        setActiveStep(inProgressWorkflows[0].currentStep);
      } else {
        setActiveWorkflow(workflows[0]);
        setActiveStep(workflows[0].currentStep || 0);
      }
    }
  }, [workflows, activeWorkflow]);
  
  const handleWorkflowDialogOpen = () => {
    setWorkflowDialogOpen(true);
  };
  
  const handleWorkflowDialogClose = () => {
    setWorkflowDialogOpen(false);
    setNewWorkflowName('');
    setNewWorkflowDescription('');
  };
  
  const handleCreateWorkflow = () => {
    if (!newWorkflowName.trim()) {
      setError('Please enter a workflow name');
      return;
    }
    
    const newWorkflow = {
      id: `workflow-${Date.now()}`,
      name: newWorkflowName,
      description: newWorkflowDescription || 'No description provided',
      created: new Date().toISOString().split('T')[0],
      lastModified: new Date().toISOString().split('T')[0],
      status: 'in_progress',
      currentStep: 0,
      steps: [],
      molecules: []
    };
    
    const updatedWorkflows = [...workflows, newWorkflow];
    setWorkflows(updatedWorkflows);
    setActiveWorkflow(newWorkflow);
    setActiveStep(0);
    
    // Save to localStorage
    localStorage.setItem('workflows', JSON.stringify(updatedWorkflows));
    
    handleWorkflowDialogClose();
  };
  
  const handleSelectWorkflow = (workflow) => {
    setActiveWorkflow(workflow);
    setActiveStep(workflow.currentStep || 0);
  };
  
  const handleNextStep = () => {
    if (!activeWorkflow || activeStep >= activeWorkflow.steps.length - 1) return;
    
    const nextStep = activeStep + 1;
    setActiveStep(nextStep);
    
    // Update workflow
    const updatedWorkflow = { ...activeWorkflow, currentStep: nextStep };
    const updatedWorkflows = workflows.map(wf => 
      wf.id === activeWorkflow.id ? updatedWorkflow : wf
    );
    
    setActiveWorkflow(updatedWorkflow);
    setWorkflows(updatedWorkflows);
    
    // Save to localStorage
    localStorage.setItem('workflows', JSON.stringify(updatedWorkflows));
  };
  
  const handlePreviousStep = () => {
    if (!activeWorkflow || activeStep <= 0) return;
    
    const prevStep = activeStep - 1;
    setActiveStep(prevStep);
    
    // Update workflow
    const updatedWorkflow = { ...activeWorkflow, currentStep: prevStep };
    const updatedWorkflows = workflows.map(wf => 
      wf.id === activeWorkflow.id ? updatedWorkflow : wf
    );
    
    setActiveWorkflow(updatedWorkflow);
    setWorkflows(updatedWorkflows);
    
    // Save to localStorage
    localStorage.setItem('workflows', JSON.stringify(updatedWorkflows));
  };
  
  const handleLaunchComponent = (componentName) => {
    if (onLaunchComponent) {
      onLaunchComponent(componentName);
    } else {
      console.log(`Launch component requested: ${componentName}`);
    }
  };
  
  const handleTaskClick = (task) => {
    setSelectedTask(task);
    setTaskDialogOpen(true);
  };
  
  const handleTaskDialogClose = () => {
    setTaskDialogOpen(false);
    setSelectedTask(null);
  };
  
  const getStatusChip = (status) => {
    let color, label;
    
    switch (status) {
      case 'completed':
        color = 'primary';
        label = 'Completed';
        break;
      case 'in_progress':
        color = 'secondary';
        label = 'In Progress';
        break;
      case 'not_started':
        color = 'default';
        label = 'Not Started';
        break;
      case 'paused':
        color = 'default';
        label = 'Paused';
        break;
      default:
        color = 'default';
        label = status.replace('_', ' ');
    }
    
    return <Chip color={color} size="small" label={label} />;
  };
  
  const getComponentIcon = (componentName) => {
    switch (componentName) {
      case 'MoleculeDesigner':
        return <ScienceIcon className={classes.taskIcon} />;
      case 'Simulations':
        return <MemoryIcon className={classes.taskIcon} />;
      case 'ChEMBL':
      case 'PubMedBioC':
        return <MenuBookIcon className={classes.taskIcon} />;
      case 'Visualization':
        return <VisibilityIcon className={classes.taskIcon} />;
      default:
        return <ScienceIcon className={classes.taskIcon} />;
    }
  };
  
  // Renders active workflow step details
  const renderActiveStepDetails = () => {
    if (!activeWorkflow || !activeWorkflow.steps || activeWorkflow.steps.length === 0) {
      return (
        <Card className={classes.card}>
          <CardContent>
            <Typography variant="h6" gutterBottom>
              No Steps Defined
            </Typography>
            <Typography variant="body2" color="textSecondary">
              This workflow doesn't have any steps defined yet. Add steps to get started.
            </Typography>
          </CardContent>
          <CardActions>
            <Button
              size="small"
              color="primary"
              startIcon={<AddIcon />}
            >
              Add First Step
            </Button>
          </CardActions>
        </Card>
      );
    }
    
    const currentStep = activeWorkflow.steps[activeStep];
    
    if (!currentStep) return null;
    
    return (
      <Card className={classes.card}>
        <div className={classes.activeTaskIndicator} />
        <CardHeader
          title={currentStep.name}
          subheader={`Step ${activeStep + 1} of ${activeWorkflow.steps.length}`}
          action={
            <Chip 
              size="small" 
              label={currentStep.status.replace('_', ' ')} 
              color={currentStep.status === 'completed' ? 'primary' : (currentStep.status === 'in_progress' ? 'secondary' : 'default')}
            />
          }
        />
        <CardContent className={classes.cardContent}>
          <Typography variant="body1" paragraph>
            {currentStep.description}
          </Typography>
          
          <div style={{ textAlign: 'center', margin: '16px 0' }}>
            {getComponentIcon(currentStep.component)}
            <Typography variant="subtitle1" gutterBottom>
              {currentStep.component}
            </Typography>
          </div>
          
          {currentStep.results && (
            <div>
              <Typography variant="subtitle2" gutterBottom>
                Results:
              </Typography>
              
              {Object.entries(currentStep.results).map(([key, value]) => (
                key !== 'pending' && (
                  <div key={key} className={classes.summaryItem}>
                    <Typography variant="body2" className={classes.summaryLabel}>
                      {key.charAt(0).toUpperCase() + key.slice(1).replace(/([A-Z])/g, ' $1')}:
                    </Typography>
                    <Typography variant="body2">
                      {Array.isArray(value) ? value.join(', ') : value.toString()}
                    </Typography>
                  </div>
                )
              ))}
            </div>
          )}
        </CardContent>
        <CardActions>
          <Button
            size="small"
            color="primary"
            startIcon={<PlayArrowIcon />}
            onClick={() => handleLaunchComponent(currentStep.component)}
          >
            Open {currentStep.component}
          </Button>
          
          <Button
            size="small"
            startIcon={<EditIcon />}
          >
            Edit Step
          </Button>
        </CardActions>
      </Card>
    );
  };
  
  // Dialog content for task details view
  const renderTaskDetails = () => {
    if (!selectedTask) return null;
    
    return (
      <div>
        <Typography variant="subtitle1" gutterBottom>
          {selectedTask.description}
        </Typography>
        
        <Divider className={classes.divider} />
        
        <Typography variant="subtitle2" gutterBottom>
          Status:
        </Typography>
        <div style={{ margin: '8px 0' }}>
          {getStatusChip(selectedTask.status)}
        </div>
        
        <Typography variant="subtitle2" gutterBottom>
          Component:
        </Typography>
        <Typography variant="body2" paragraph>
          {selectedTask.component}
        </Typography>
        
        {selectedTask.results && (
          <div>
            <Typography variant="subtitle2" gutterBottom>
              Results:
            </Typography>
            
            {Object.entries(selectedTask.results).map(([key, value]) => (
              key !== 'pending' && (
                <div key={key} className={classes.summaryItem}>
                  <Typography variant="body2" className={classes.summaryLabel}>
                    {key.charAt(0).toUpperCase() + key.slice(1).replace(/([A-Z])/g, ' $1')}:
                  </Typography>
                  <Typography variant="body2">
                    {Array.isArray(value) ? value.join(', ') : value.toString()}
                  </Typography>
                </div>
              )
            ))}
          </div>
        )}
      </div>
    );
  };
  
  // Renders the workflow stepper
  const renderWorkflowStepper = () => {
    if (!activeWorkflow || !activeWorkflow.steps || activeWorkflow.steps.length === 0) {
      return (
        <Typography variant="body2" color="textSecondary" style={{ padding: 16 }}>
          No steps defined in this workflow.
        </Typography>
      );
    }
    
    return (
      <Stepper 
        activeStep={activeStep} 
        orientation="vertical" 
        className={classes.workflowStepper}
        nonLinear
      >
        {activeWorkflow.steps.map((step, index) => (
          <Step key={step.id}>
            <StepLabel 
              className={index === activeStep ? classes.workflowStepActive : classes.workflowStep}
              onClick={() => setActiveStep(index)}
            >
              <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                <div>{step.name}</div>
                {getStatusChip(step.status)}
              </div>
            </StepLabel>
            <StepContent>
              <Typography variant="body2" color="textSecondary">
                {step.description}
              </Typography>
              <div style={{ marginTop: 8 }}>
                <Button
                  size="small"
                  color="primary"
                  onClick={() => handleLaunchComponent(step.component)}
                >
                  Open {step.component}
                </Button>
              </div>
            </StepContent>
          </Step>
        ))}
      </Stepper>
    );
  };
  
  return (
    <div className={classes.root}>
      <Typography variant="h5" className={classes.title}>
        Workflow Integration
      </Typography>
      
      {error && (
        <Alert severity="error" style={{ marginBottom: 16 }}>
          {error}
        </Alert>
      )}
      
      {loading ? (
        <div style={{ textAlign: 'center', padding: 40 }}>
          <CircularProgress />
          <Typography style={{ marginTop: 16 }}>
            Loading workflows...
          </Typography>
        </div>
      ) : (
        <Grid container spacing={3}>
          <Grid item xs={12} md={4}>
            <Paper className={classes.paper}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 16 }}>
                <Typography variant="h6">
                  Workflows
                </Typography>
                <Button
                  variant="outlined"
                  color="primary"
                  size="small"
                  startIcon={<AddIcon />}
                  onClick={handleWorkflowDialogOpen}
                >
                  New
                </Button>
              </div>
              
              {workflows.length === 0 ? (
                <Typography variant="body2" color="textSecondary">
                  No workflows created yet. Create a new workflow to get started.
                </Typography>
              ) : (
                <List>
                  {workflows.map((workflow) => (
                    <ListItem
                      key={workflow.id}
                      button
                      selected={activeWorkflow && activeWorkflow.id === workflow.id}
                      onClick={() => handleSelectWorkflow(workflow)}
                    >
                      <ListItemText
                        primary={workflow.name}
                        secondary={
                          <React.Fragment>
                            <Typography variant="body2" color="textSecondary" noWrap>
                              {workflow.description}
                            </Typography>
                            <Typography variant="caption" color="textSecondary">
                              Last modified: {workflow.lastModified}
                            </Typography>
                          </React.Fragment>
                        }
                      />
                      <Chip
                        size="small"
                        label={workflow.status.replace('_', ' ')}
                        color={workflow.status === 'completed' ? 'primary' : (workflow.status === 'in_progress' ? 'secondary' : 'default')}
                        className={classes.taskStatusChip}
                      />
                    </ListItem>
                  ))}
                </List>
              )}
            </Paper>
          </Grid>
          
          <Grid item xs={12} md={8}>
            {activeWorkflow ? (
              <div>
                <Paper className={classes.workflowInfo}>
                  <div style={{ padding: 16 }}>
                    <Typography variant="h6" gutterBottom>
                      {activeWorkflow.name}
                    </Typography>
                    <Typography variant="body2" color="textSecondary" paragraph>
                      {activeWorkflow.description}
                    </Typography>
                    
                    <Grid container spacing={2}>
                      <Grid item xs={6}>
                        <Typography variant="body2">
                          <strong>Created:</strong> {activeWorkflow.created}
                        </Typography>
                      </Grid>
                      <Grid item xs={6}>
                        <Typography variant="body2">
                          <strong>Last Modified:</strong> {activeWorkflow.lastModified}
                        </Typography>
                      </Grid>
                      <Grid item xs={6}>
                        <Typography variant="body2">
                          <strong>Status:</strong> {activeWorkflow.status.replace('_', ' ')}
                        </Typography>
                      </Grid>
                      <Grid item xs={6}>
                        <Typography variant="body2">
                          <strong>Molecules:</strong> {activeWorkflow.molecules?.length || 0}
                        </Typography>
                      </Grid>
                    </Grid>
                  </div>
                </Paper>
                
                <Grid container spacing={3}>
                  <Grid item xs={12} md={4}>
                    <Paper className={classes.paper}>
                      <Typography variant="subtitle1" gutterBottom>
                        Workflow Steps
                      </Typography>
                      {renderWorkflowStepper()}
                      
                      <div style={{ 
                        display: 'flex', 
                        justifyContent: 'space-between', 
                        marginTop: 16, 
                        paddingTop: 16,
                        borderTop: '1px solid rgba(0, 0, 0, 0.12)'
                      }}>
                        <Button
                          variant="outlined"
                          size="small"
                          onClick={handlePreviousStep}
                          disabled={activeStep === 0}
                        >
                          Previous
                        </Button>
                        <Button
                          variant="outlined"
                          size="small"
                          onClick={handleNextStep}
                          disabled={!activeWorkflow.steps || activeStep >= activeWorkflow.steps.length - 1}
                        >
                          Next
                        </Button>
                      </div>
                    </Paper>
                  </Grid>
                  
                  <Grid item xs={12} md={8}>
                    {renderActiveStepDetails()}
                    
                    {activeWorkflow.molecules && activeWorkflow.molecules.length > 0 && (
                      <Card className={classes.workflowCard} style={{ marginTop: 16 }}>
                        <CardHeader title="Molecules in Workflow" />
                        <CardContent>
                          <List dense>
                            {activeWorkflow.molecules.map((molecule) => (
                              <ListItem key={molecule.id}>
                                <ListItemIcon>
                                  <ScienceIcon color="primary" />
                                </ListItemIcon>
                                <ListItemText
                                  primary={molecule.name}
                                  secondary={molecule.smiles}
                                />
                                <Chip
                                  size="small"
                                  label={molecule.status.replace('_', ' ')}
                                  color={molecule.status.includes('complete') ? 'primary' : 'default'}
                                />
                              </ListItem>
                            ))}
                          </List>
                        </CardContent>
                        <CardActions>
                          <Button
                            size="small"
                            color="primary"
                            startIcon={<VisibilityIcon />}
                            onClick={() => handleLaunchComponent('MoleculeStorage')}
                          >
                            View All Molecules
                          </Button>
                        </CardActions>
                      </Card>
                    )}
                  </Grid>
                </Grid>
              </div>
            ) : (
              <Paper className={classes.paper} style={{ textAlign: 'center', padding: 40 }}>
                <Typography variant="h6" gutterBottom>
                  No Workflow Selected
                </Typography>
                <Typography variant="body2" color="textSecondary" paragraph>
                  Select a workflow from the list or create a new one to get started.
                </Typography>
                <Button
                  variant="contained"
                  color="primary"
                  startIcon={<AddIcon />}
                  onClick={handleWorkflowDialogOpen}
                >
                  Create New Workflow
                </Button>
              </Paper>
            )}
          </Grid>
          
          <Grid item xs={12}>
            <Paper className={classes.paper}>
              <Typography variant="h6" gutterBottom>
                Component Connections
              </Typography>
              <Grid container spacing={2} alignItems="center" justifyContent="center">
                <Grid item xs={12} sm={2}>
                  <div style={{ textAlign: 'center' }}>
                    <MenuBookIcon className={classes.taskIcon} />
                    <Typography variant="subtitle1">Literature</Typography>
                    <Button
                      size="small"
                      variant="outlined"
                      onClick={() => handleLaunchComponent('PubMedBioC')}
                      style={{ marginTop: 8 }}
                    >
                      Open
                    </Button>
                  </div>
                </Grid>
                
                <Grid item xs={12} sm={1}>
                  <ArrowForwardIcon className={classes.arrowIcon} />
                </Grid>
                
                <Grid item xs={12} sm={2}>
                  <div style={{ textAlign: 'center' }}>
                    <MenuBookIcon className={classes.taskIcon} />
                    <Typography variant="subtitle1">ChEMBL</Typography>
                    <Button
                      size="small"
                      variant="outlined"
                      onClick={() => handleLaunchComponent('ChEMBL')}
                      style={{ marginTop: 8 }}
                    >
                      Open
                    </Button>
                  </div>
                </Grid>
                
                <Grid item xs={12} sm={1}>
                  <ArrowForwardIcon className={classes.arrowIcon} />
                </Grid>
                
                <Grid item xs={12} sm={2}>
                  <div style={{ textAlign: 'center' }}>
                    <ScienceIcon className={classes.taskIconSecondary} />
                    <Typography variant="subtitle1">Designer</Typography>
                    <Button
                      size="small"
                      variant="outlined"
                      color="secondary"
                      onClick={() => handleLaunchComponent('MoleculeDesigner')}
                      style={{ marginTop: 8 }}
                    >
                      Open
                    </Button>
                  </div>
                </Grid>
                
                <Grid item xs={12} sm={1}>
                  <ArrowForwardIcon className={classes.arrowIcon} />
                </Grid>
                
                <Grid item xs={12} sm={2}>
                  <div style={{ textAlign: 'center' }}>
                    <MemoryIcon className={classes.taskIcon} />
                    <Typography variant="subtitle1">Simulations</Typography>
                    <Button
                      size="small"
                      variant="outlined"
                      onClick={() => handleLaunchComponent('Simulations')}
                      style={{ marginTop: 8 }}
                    >
                      Open
                    </Button>
                  </div>
                </Grid>
              </Grid>
              
              <Divider style={{ margin: '24px 0 16px' }} />
              
              <Typography variant="subtitle2" gutterBottom>
                Quick Actions
              </Typography>
              
              <Grid container spacing={2}>
                <Grid item xs={12} sm={6} md={3}>
                  <Button
                    variant="outlined"
                    fullWidth
                    size="small"
                    startIcon={<SearchIcon />}
                    onClick={() => handleLaunchComponent('MoleculeStorage')}
                    className={classes.button}
                  >
                    Stored Molecules
                  </Button>
                </Grid>
                
                <Grid item xs={12} sm={6} md={3}>
                  <Button
                    variant="outlined"
                    fullWidth
                    size="small"
                    startIcon={<VisibilityIcon />}
                    onClick={() => handleLaunchComponent('Visualization')}
                    className={classes.button}
                  >
                    Visualize Data
                  </Button>
                </Grid>
                
                <Grid item xs={12} sm={6} md={3}>
                  <Button
                    variant="outlined"
                    fullWidth
                    size="small"
                    startIcon={<SaveIcon />}
                    className={classes.button}
                  >
                    Export Results
                  </Button>
                </Grid>
                
                <Grid item xs={12} sm={6} md={3}>
                  <Button
                    variant="outlined"
                    fullWidth
                    size="small"
                    startIcon={<SettingsIcon />}
                    className={classes.button}
                  >
                    Workflow Settings
                  </Button>
                </Grid>
              </Grid>
            </Paper>
          </Grid>
        </Grid>
      )}
      
      {/* New Workflow Dialog */}
      <Dialog
        open={workflowDialogOpen}
        onClose={handleWorkflowDialogClose}
        aria-labelledby="new-workflow-dialog-title"
      >
        <DialogTitle id="new-workflow-dialog-title">Create New Workflow</DialogTitle>
        <DialogContent>
          <TextField
            autoFocus
            margin="dense"
            label="Workflow Name"
            fullWidth
            value={newWorkflowName}
            onChange={(e) => setNewWorkflowName(e.target.value)}
            required
          />
          <TextField
            margin="dense"
            label="Description"
            fullWidth
            multiline
            rows={3}
            value={newWorkflowDescription}
            onChange={(e) => setNewWorkflowDescription(e.target.value)}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={handleWorkflowDialogClose} color="primary">
            Cancel
          </Button>
          <Button onClick={handleCreateWorkflow} color="primary" variant="contained">
            Create
          </Button>
        </DialogActions>
      </Dialog>
      
      {/* Task Details Dialog */}
      <Dialog
        open={taskDialogOpen}
        onClose={handleTaskDialogClose}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>
          {selectedTask?.name || 'Task Details'}
        </DialogTitle>
        <DialogContent>
          {renderTaskDetails()}
        </DialogContent>
        <DialogActions>
          <Button onClick={handleTaskDialogClose} color="primary">
            Close
          </Button>
          <Button 
            color="primary" 
            variant="contained"
            onClick={() => {
              handleLaunchComponent(selectedTask?.component);
              handleTaskDialogClose();
            }}
            disabled={!selectedTask?.component}
          >
            Open Component
          </Button>
        </DialogActions>
      </Dialog>
    </div>
  );
};

export default WorkflowIntegration; 