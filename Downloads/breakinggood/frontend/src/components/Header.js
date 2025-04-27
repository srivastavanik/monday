import React from 'react';
import { AppBar, Toolbar, Typography, Button, IconButton, makeStyles } from '@material-ui/core';
import SearchIcon from '@material-ui/icons/Search';
import HelpOutlineIcon from '@material-ui/icons/HelpOutline';
import NotificationsNoneIcon from '@material-ui/icons/NotificationsNone';
import MenuIcon from '@material-ui/icons/Menu';

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
  },
  appBar: {
    zIndex: theme.zIndex.drawer + 1,
    backgroundColor: '#000',
    boxShadow: 'none',
    borderBottom: '1px solid #e0e0e0',
  },
  menuButton: {
    marginRight: theme.spacing(2),
    [theme.breakpoints.up('md')]: {
      display: 'none',
    },
  },
  title: {
    flexGrow: 1,
    fontWeight: 500,
    letterSpacing: '1px',
    color: '#fff',
  },
  logo: {
    marginRight: theme.spacing(2),
    height: '32px',
  },
  searchButton: {
    marginRight: theme.spacing(2),
    color: '#fff',
  },
  iconButton: {
    color: '#fff',
    marginLeft: theme.spacing(1),
  },
  userName: {
    marginLeft: theme.spacing(2),
    marginRight: theme.spacing(1),
    color: '#fff',
  },
}));

const Header = () => {
  const classes = useStyles();

  return (
    <div className={classes.root}>
      <AppBar position="static" className={classes.appBar}>
        <Toolbar>
          <IconButton
            edge="start"
            className={classes.menuButton}
            color="inherit"
            aria-label="menu"
          >
            <MenuIcon />
          </IconButton>
          
          <Typography variant="h6" className={classes.title}>
            BREAKING GOOD
          </Typography>
          
          <Button
            startIcon={<SearchIcon />}
            className={classes.searchButton}
            aria-label="search"
          >
            Search
          </Button>
          
          <IconButton className={classes.iconButton} aria-label="help">
            <HelpOutlineIcon />
          </IconButton>
          
          <IconButton className={classes.iconButton} aria-label="notifications">
            <NotificationsNoneIcon />
          </IconButton>
          
          <Typography variant="body2" className={classes.userName}>
            Dr. Jane Researcher
          </Typography>
        </Toolbar>
      </AppBar>
    </div>
  );
};

export default Header; 