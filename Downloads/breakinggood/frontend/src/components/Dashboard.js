import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import { 
  Typography, 
  Grid, 
  Paper, 
  Card, 
  CardContent, 
  CardActions, 
  Button, 
  Divider,
  makeStyles,
  LinearProgress
} from '@material-ui/core';
import TrendingUpIcon from '@material-ui/icons/TrendingUp';
import TrendingDownIcon from '@material-ui/icons/TrendingDown';
import FiberNewIcon from '@material-ui/icons/FiberNew';
import ScienceIcon from '@material-ui/icons/Science';
import MenuBookIcon from '@material-ui/icons/MenuBook';
import AssessmentIcon from '@material-ui/icons/Assessment';

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
    padding: theme.spacing(2),
    height: '100%',
  },
  card: {
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    transition: 'transform 0.3s',
    '&:hover': {
      transform: 'translateY(-5px)',
      boxShadow: '0 4px 10px rgba(0,0,0,0.1)',
    },
  },
  cardContent: {
    flexGrow: 1,
  },
  statCard: {
    padding: theme.spacing(2),
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    backgroundColor: '#f5f5f5',
    height: '100%',
  },
  statValue: {
    fontWeight: 500,
    fontSize: '1.5rem',
  },
  statLabel: {
    color: theme.palette.text.secondary,
    fontSize: '0.875rem',
  },
  trendIcon: {
    marginLeft: theme.spacing(1),
    fontSize: '1rem',
  },
  trendUp: {
    color: '#4caf50',
  },
  trendDown: {
    color: '#f44336',
  },
  icon: {
    fontSize: '2.5rem',
    color: theme.palette.grey[500],
  },
  divider: {
    margin: theme.spacing(2, 0),
  },
  progress: {
    height: 8,
    borderRadius: 4,
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1),
  },
  recentItem: {
    padding: theme.spacing(1.5),
    borderLeft: `3px solid ${theme.palette.primary.main}`,
    marginBottom: theme.spacing(1),
    backgroundColor: '#f9f9f9',
  },
  recentTitle: {
    fontWeight: 500,
  },
  recentMeta: {
    fontSize: '0.75rem',
    color: theme.palette.text.secondary,
  },
  newIcon: {
    color: theme.palette.secondary.main,
    fontSize: '0.875rem',
    marginLeft: theme.spacing(1),
    verticalAlign: 'middle',
  },
  cardAction: {
    justifyContent: 'flex-end',
  },
  sectionTitle: {
    marginBottom: theme.spacing(2),
    marginTop: theme.spacing(4),
  },
}));

const Dashboard = () => {
  const classes = useStyles();
  const [loading, setLoading] = useState(true);
  const [stats, setStats] = useState({
    pendingMolecules: 0,
    completedDesigns: 0,
    literatureCount: 0,
    simulationRuns: 0
  });
  
  useEffect(() => {
    // Simulate loading data
    const timer = setTimeout(() => {
      setStats({
        pendingMolecules: 4,
        completedDesigns: 12,
        literatureCount: 347,
        simulationRuns: 28
      });
      setLoading(false);
    }, 1000);
    
    return () => clearTimeout(timer);
  }, []);
  
  const recentMolecules = [
    {
      id: 'mol-1234',
      name: 'Modified Methylphenidate Analog',
      date: '2023-06-10',
      status: 'In Simulation',
      isNew: true
    },
    {
      id: 'mol-1233',
      name: 'Amphetamine Derivative XJ-42',
      date: '2023-06-08',
      status: 'FDA Pathway Analysis',
      isNew: false
    },
    {
      id: 'mol-1232',
      name: 'Novel Dopamine Reuptake Inhibitor',
      date: '2023-06-05',
      status: 'Completed',
      isNew: false
    }
  ];
  
  const projectProgress = [
    {
      name: 'Project Alpha: Methylphenidate Alternative',
      progress: 78,
      phase: 'Production Feasibility'
    },
    {
      name: 'Project Beta: Non-stimulant ADHD Treatment',
      progress: 45,
      phase: 'Simulation'
    },
    {
      name: 'Project Gamma: Extended Release Formulation',
      progress: 20,
      phase: 'Initial Design'
    }
  ];

  return (
    <div className={classes.root}>
      <Typography variant="h4" className={classes.title}>
        Dashboard
      </Typography>
      
      {loading ? (
        <LinearProgress />
      ) : (
        <>
          <Grid container spacing={3}>
            <Grid item xs={12} sm={6} md={3}>
              <Paper className={classes.statCard}>
                <div>
                  <Typography className={classes.statLabel}>
                    Pending Molecules
                  </Typography>
                  <Typography className={classes.statValue}>
                    {stats.pendingMolecules}
                    <TrendingUpIcon className={`${classes.trendIcon} ${classes.trendUp}`} />
                  </Typography>
                </div>
                <ScienceIcon className={classes.icon} />
              </Paper>
            </Grid>
            
            <Grid item xs={12} sm={6} md={3}>
              <Paper className={classes.statCard}>
                <div>
                  <Typography className={classes.statLabel}>
                    Completed Designs
                  </Typography>
                  <Typography className={classes.statValue}>
                    {stats.completedDesigns}
                    <TrendingUpIcon className={`${classes.trendIcon} ${classes.trendUp}`} />
                  </Typography>
                </div>
                <AssessmentIcon className={classes.icon} />
              </Paper>
            </Grid>
            
            <Grid item xs={12} sm={6} md={3}>
              <Paper className={classes.statCard}>
                <div>
                  <Typography className={classes.statLabel}>
                    Literature References
                  </Typography>
                  <Typography className={classes.statValue}>
                    {stats.literatureCount}
                    <TrendingUpIcon className={`${classes.trendIcon} ${classes.trendUp}`} />
                  </Typography>
                </div>
                <MenuBookIcon className={classes.icon} />
              </Paper>
            </Grid>
            
            <Grid item xs={12} sm={6} md={3}>
              <Paper className={classes.statCard}>
                <div>
                  <Typography className={classes.statLabel}>
                    Simulation Runs
                  </Typography>
                  <Typography className={classes.statValue}>
                    {stats.simulationRuns}
                    <TrendingDownIcon className={`${classes.trendIcon} ${classes.trendDown}`} />
                  </Typography>
                </div>
                <AssessmentIcon className={classes.icon} />
              </Paper>
            </Grid>
          </Grid>
          
          <Typography variant="h5" className={classes.sectionTitle}>
            Project Progress
          </Typography>
          
          <Grid container spacing={3}>
            {projectProgress.map((project, index) => (
              <Grid item xs={12} md={4} key={index}>
                <Paper className={classes.paper}>
                  <Typography variant="h6">{project.name}</Typography>
                  <Typography variant="body2" color="textSecondary">
                    Current Phase: {project.phase}
                  </Typography>
                  <LinearProgress 
                    variant="determinate" 
                    value={project.progress} 
                    className={classes.progress}
                  />
                  <Typography variant="body2" align="right">
                    {project.progress}% Complete
                  </Typography>
                </Paper>
              </Grid>
            ))}
          </Grid>
          
          <Grid container spacing={3} style={{ marginTop: '16px' }}>
            <Grid item xs={12} md={8}>
              <Paper className={classes.paper}>
                <Typography variant="h6">Recent Molecules</Typography>
                <Divider className={classes.divider} />
                
                {recentMolecules.map((molecule) => (
                  <div className={classes.recentItem} key={molecule.id}>
                    <Typography className={classes.recentTitle}>
                      {molecule.name}
                      {molecule.isNew && <FiberNewIcon className={classes.newIcon} />}
                    </Typography>
                    <Typography className={classes.recentMeta}>
                      ID: {molecule.id} • Date: {molecule.date} • Status: {molecule.status}
                    </Typography>
                  </div>
                ))}
                
                <Button 
                  component={Link} 
                  to="/molecule-designer" 
                  color="primary" 
                  size="small"
                >
                  View All Molecules
                </Button>
              </Paper>
            </Grid>
            
            <Grid item xs={12} md={4}>
              <Card className={classes.card}>
                <CardContent className={classes.cardContent}>
                  <ScienceIcon style={{ fontSize: 40, marginBottom: 16 }} />
                  <Typography variant="h6" gutterBottom>
                    Create New Molecule
                  </Typography>
                  <Typography variant="body2" color="textSecondary">
                    Start designing a new molecule with the help of our AI-powered tools.
                    Access literature references, simulation tools, and regulatory analysis.
                  </Typography>
                </CardContent>
                <CardActions className={classes.cardAction}>
                  <Button 
                    component={Link} 
                    to="/molecule-designer" 
                    color="primary" 
                    variant="contained"
                  >
                    Get Started
                  </Button>
                </CardActions>
              </Card>
            </Grid>
          </Grid>
        </>
      )}
    </div>
  );
};

export default Dashboard; 