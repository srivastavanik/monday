import React from 'react';
import { Link, useLocation } from 'react-router-dom';
import { 
  Drawer, 
  List, 
  ListItem, 
  ListItemIcon, 
  ListItemText, 
  Divider, 
  makeStyles 
} from '@material-ui/core';
import DashboardIcon from '@material-ui/icons/Dashboard';
import ScienceIcon from '@material-ui/icons/Science';
import MenuBookIcon from '@material-ui/icons/MenuBook';
import BarChartIcon from '@material-ui/icons/BarChart';
import AssignmentIcon from '@material-ui/icons/Assignment';
import CompareArrowsIcon from '@material-ui/icons/CompareArrows';
import SettingsIcon from '@material-ui/icons/Settings';
import HelpIcon from '@material-ui/icons/Help';

const drawerWidth = 240;

const useStyles = makeStyles((theme) => ({
  drawer: {
    width: drawerWidth,
    flexShrink: 0,
    [theme.breakpoints.down('sm')]: {
      display: 'none',
    },
  },
  drawerPaper: {
    width: drawerWidth,
    borderRight: '1px solid #e0e0e0',
    backgroundColor: '#f5f5f5',
  },
  toolbar: theme.mixins.toolbar,
  listItem: {
    borderLeft: '3px solid transparent',
    '&:hover': {
      backgroundColor: '#eeeeee',
    },
  },
  activeListItem: {
    borderLeft: '3px solid #000000',
    backgroundColor: '#eeeeee',
    '&:hover': {
      backgroundColor: '#e0e0e0',
    },
  },
  listItemIcon: {
    minWidth: 40,
    color: '#757575',
  },
  activeIcon: {
    color: '#000000',
  },
  subheader: {
    textTransform: 'uppercase',
    fontWeight: 500,
    fontSize: '0.75rem',
    color: '#9e9e9e',
    padding: theme.spacing(2, 2, 1, 2),
    letterSpacing: '1px',
  },
}));

const Sidebar = () => {
  const classes = useStyles();
  const location = useLocation();

  const isActive = (path) => {
    return location.pathname === path;
  };

  return (
    <Drawer
      className={classes.drawer}
      variant="permanent"
      classes={{
        paper: classes.drawerPaper,
      }}
    >
      <div className={classes.toolbar} />
      
      <List component="nav" aria-label="main navigation">
        <ListItem 
          button 
          component={Link} 
          to="/" 
          className={isActive('/') ? classes.activeListItem : classes.listItem}
        >
          <ListItemIcon className={`${classes.listItemIcon} ${isActive('/') ? classes.activeIcon : ''}`}>
            <DashboardIcon />
          </ListItemIcon>
          <ListItemText primary="Dashboard" />
        </ListItem>
        
        <ListItem 
          button 
          component={Link} 
          to="/molecule-designer" 
          className={isActive('/molecule-designer') ? classes.activeListItem : classes.listItem}
        >
          <ListItemIcon className={`${classes.listItemIcon} ${isActive('/molecule-designer') ? classes.activeIcon : ''}`}>
            <ScienceIcon />
          </ListItemIcon>
          <ListItemText primary="Molecule Designer" />
        </ListItem>
        
        <ListItem 
          button 
          component={Link} 
          to="/literature" 
          className={isActive('/literature') ? classes.activeListItem : classes.listItem}
        >
          <ListItemIcon className={`${classes.listItemIcon} ${isActive('/literature') ? classes.activeIcon : ''}`}>
            <MenuBookIcon />
          </ListItemIcon>
          <ListItemText primary="Literature Explorer" />
        </ListItem>
        
        <ListItem 
          button 
          component={Link} 
          to="/simulation" 
          className={isActive('/simulation') ? classes.activeListItem : classes.listItem}
        >
          <ListItemIcon className={`${classes.listItemIcon} ${isActive('/simulation') ? classes.activeIcon : ''}`}>
            <BarChartIcon />
          </ListItemIcon>
          <ListItemText primary="Simulation Panel" />
        </ListItem>
        
        <ListItem 
          button 
          component={Link} 
          to="/regulatory" 
          className={isActive('/regulatory') ? classes.activeListItem : classes.listItem}
        >
          <ListItemIcon className={`${classes.listItemIcon} ${isActive('/regulatory') ? classes.activeIcon : ''}`}>
            <AssignmentIcon />
          </ListItemIcon>
          <ListItemText primary="Regulatory Analysis" />
        </ListItem>
        
        <ListItem 
          button 
          component={Link} 
          to="/comparison" 
          className={isActive('/comparison') ? classes.activeListItem : classes.listItem}
        >
          <ListItemIcon className={`${classes.listItemIcon} ${isActive('/comparison') ? classes.activeIcon : ''}`}>
            <CompareArrowsIcon />
          </ListItemIcon>
          <ListItemText primary="Comparison Tool" />
        </ListItem>
      </List>
      
      <Divider />
      
      <div className={classes.subheader}>Settings</div>
      
      <List component="nav" aria-label="settings">
        <ListItem button className={classes.listItem}>
          <ListItemIcon className={classes.listItemIcon}>
            <SettingsIcon />
          </ListItemIcon>
          <ListItemText primary="Preferences" />
        </ListItem>
        
        <ListItem button className={classes.listItem}>
          <ListItemIcon className={classes.listItemIcon}>
            <HelpIcon />
          </ListItemIcon>
          <ListItemText primary="Help & Documentation" />
        </ListItem>
      </List>
    </Drawer>
  );
};

export default Sidebar; 