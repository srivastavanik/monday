import React, { useState } from 'react';
import { 
  ThemeProvider, 
  createTheme, 
  makeStyles,
  CssBaseline,
  AppBar,
  Toolbar,
  Typography,
  Tabs,
  Tab,
  Box,
  Container,
  IconButton,
  Button,
  Drawer,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Divider
} from '@material-ui/core';
import MenuIcon from '@material-ui/icons/Menu';
import HomeIcon from '@material-ui/icons/Home';
import BuildIcon from '@material-ui/icons/Build'; // Using Build icon instead of Science
import DescriptionIcon from '@material-ui/icons/Description';
import CompareArrowsIcon from '@material-ui/icons/CompareArrows';
import TimelineIcon from '@material-ui/icons/Timeline';
import AccountCircleIcon from '@material-ui/icons/AccountCircle';

// Import all the components
import LandingPage from './components/LandingPage';
import MoleculeDesigner from './components/MoleculeDesigner';
import SimulationPanel from './components/SimulationPanel';
import LiteratureExplorer from './components/LiteratureExplorer';
import ComparisonTool from './components/ComparisonTool';
import RegulatoryAnalysis from './components/RegulatoryAnalysis';

// Create TabPanel component
function TabPanel(props) {
  const { children, value, index, ...other } = props;

  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`simple-tabpanel-${index}`}
      aria-labelledby={`simple-tab-${index}`}
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

// Define Anthropic/Claude-inspired theme (Light)
const landingTheme = createTheme({
  palette: {
    type: 'light',
    primary: {
      main: '#E27B58', // Terracotta Orange
      light: '#F29E7A', // Lighter accent for hover
      dark: '#C4644A', // Darker accent for hover
      contrastText: '#FFFFFF',
    },
    secondary: {
      main: '#4F4F4F', // Warm Gray (used for secondary text)
      contrastText: '#FFFFFF',
    },
    background: {
      default: '#FBF9F7', // Off-white
      paper: '#FFFFFF',
    },
    text: {
      primary: '#111111', // Charcoal Black
      secondary: '#4F4F4F', // Warm Gray
      disabled: '#BDBDBD',
    },
    divider: '#EAE8E4', // Very light gray
  },
  typography: {
    // Base font size is controlled by html global style (112.5% = 18px)
    // Define variants relative to this base
    fontFamily: [
      'Tiempos', 'ui-serif', 'Georgia', 'Cambria', '"Times New Roman"', 'Times', 'serif',
      'Styrene', 'Inter', 'Helvetica Neue', 'Arial', 'sans-serif'
    ].join(','),
    h1: {
      fontFamily: 'ui-serif, Georgia, Cambria, "Times New Roman", Times, serif',
      // Adjust rem based on new root size (18px)
      // 4rem on desktop (4 * 18 = 72px), 3rem on mobile (3*18=54px), 4.5rem on larger screens (4.5*18=81px)
      fontSize: '3rem', // Start with mobile size
      fontWeight: 400,
      lineHeight: 1.1, // Tightened line height
      color: '#111111',
      // Responsive styles applied via makeStyles or media queries in component
    },
    // H2-H6 use Styrene fallback
    h2: { fontFamily: 'Styrene, sans-serif', fontSize: '2.2rem', fontWeight: 400, color: '#111111' }, // Adjusted rem
    h3: { fontFamily: 'Styrene, sans-serif', fontSize: '1.8rem', fontWeight: 400, color: '#111111' }, // Adjusted rem
    h4: { fontFamily: 'Styrene, sans-serif', fontSize: '1.6rem', fontWeight: 400, color: '#111111' }, // Adjusted rem
    h5: { fontFamily: 'Styrene, sans-serif', fontSize: '1.4rem', fontWeight: 400, color: '#111111' }, // Adjusted rem
    h6: { fontFamily: 'Styrene, sans-serif', fontSize: '1.2rem', fontWeight: 400, color: '#111111' }, // Adjusted rem
    // Body/Tagline uses sans-serif
    body1: { // Used for tagline
      fontFamily: 'Inter, Helvetica Neue, Arial, sans-serif',
      fontSize: '1rem', // Start with mobile size (1rem * 18px = 18px, spec is 16px? using 1rem for baseline)
      lineHeight: 1.5,
      color: '#4F4F4F', 
      letterSpacing: '0.5px',
    },
    body2: { // For general paragraph text
      fontFamily: 'Tiempos, ui-serif, Georgia, serif', // Use serif for body copy
      fontSize: '1rem', // 18px body
      lineHeight: 1.6, 
      color: '#3D3929', // Use primary text color
    },
    // Nav/Button Text
    button: {
      fontFamily: 'ui-serif, Georgia, Cambria, "Times New Roman", Times, serif',
      fontSize: '1rem', // 18px based on new root
      fontWeight: 400,
      lineHeight: 1.4,
      textTransform: 'none',
      letterSpacing: '0em',
    },
    caption: {
      fontFamily: 'Inter, Helvetica Neue, Arial, sans-serif',
      fontSize: '0.77rem', // ~14px relative to 18px base
      lineHeight: 1.4,
      color: '#6A6A6A',
    }
  },
  shape: {
    borderRadius: 6, // Button radius
  },
  overrides: {
    MuiCssBaseline: {
        '@global': {
            html: {
                fontSize: '112.5%', // 18px base font size
                scrollBehavior: 'smooth',
            },
            body: {
                 backgroundColor: '#FBF9F7', // Use theme default background
            },
            '@keyframes fadeIn': {
                'from': { opacity: 0, transform: 'translateY(10px)' },
                'to': { opacity: 1, transform: 'translateY(0)' },
            },
             '@keyframes bounce': {
                'from': { transform: 'translateY(0)' },
                'to': { transform: 'translateY(8px)' },
            }
        }
    },
    MuiAppBar: {
      colorPrimary: {
        backgroundColor: 'transparent', 
        color: '#111111', 
        boxShadow: 'none', 
        borderBottom: '1px solid #EAE8E4',
      },
    },
    MuiButton: {
      root: {
        borderRadius: 6,
        padding: '12px 32px', // Updated padding: 0.75rem -> 12px, 2rem -> 32px
      },
      containedPrimary: {
        backgroundColor: '#E27B58', // Terracotta Orange
        color: '#FFFFFF',
        boxShadow: 'none',
        '&:hover': {
          backgroundColor: '#C4644A', // Darken accent
          boxShadow: 'none',
        },
      },
      outlinedPrimary: { // This corresponds to the secondary button style
          borderColor: '#111111',
          color: '#111111',
          padding: '11px 31px', // Account for border
          '&:hover': {
            backgroundColor: 'rgba(17, 17, 17, 0.04)'
          }
      },
      textPrimary: {
          color: '#E27B58',
         '&:hover': {
           backgroundColor: 'rgba(226, 123, 88, 0.04)'
         }
      }
    },
    MuiPaper: {
      rounded: {
        borderRadius: 12, // Card radius
      },
      elevation1: { boxShadow: '0px 1px 4px rgba(0,0,0,0.08)' }, // Subtle shadow
    },
    MuiOutlinedInput: {
        root: {
            borderRadius: 16, // Input radius
            backgroundColor: '#FFFFFF', 
            '&$notchedOutline': { // Correct way to target the notch
                borderColor: '#EAE8E4',
            },
            '&:hover $notchedOutline': {
                borderColor: '#C4C1BB',
            },
            '&.Mui-focused $notchedOutline': {
                borderColor: '#E27B58', 
                borderWidth: '1px',
            },
        },
        input: {
            padding: '16px', 
            fontSize: '1rem', // Input font size (1 * 18px = 18px)
        }
    },
     MuiInputLabel: {
         outlined: {
             fontSize: '1rem', // Match input field size
            transform: 'translate(16px, 18px) scale(1)',
            '&.MuiInputLabel-shrink': {
                transform: 'translate(14px, -6px) scale(0.75)',
            }
         }
     }
  },
});

// Define styles AFTER the theme
const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    display: 'flex',
    flexDirection: 'column',
    minHeight: '100vh',
  },
  appBar: {
    zIndex: theme.zIndex.drawer + 1,
    backgroundColor: theme.palette.background.paper, // Use paper for AppBar BG
    color: theme.palette.text.primary, // Use primary text for AppBar
    boxShadow: theme.shadows[1], // Use subtle shadow
  },
  title: {
    flexGrow: 1,
  },
  drawer: {
    width: 240,
    flexShrink: 0,
  },
  drawerPaper: {
    width: 240,
    backgroundColor: theme.palette.background.default, // Off-white drawer
    borderRight: 'none',
  },
  drawerContainer: {
    overflow: 'auto',
  },
  content: {
    flexGrow: 1,
    padding: 0,
    backgroundColor: theme.palette.background.default, // Main background
  },
  toolbar: theme.mixins.toolbar,
  tabsContainer: {
    backgroundColor: theme.palette.background.paper, // Paper for tabs container
    boxShadow: theme.shadows[1], 
    borderBottom: `1px solid ${theme.palette.divider}`, // Apply border here
  },
  tab: {
    minHeight: 64, // Adjust height slightly
    color: theme.palette.text.secondary,
    '&.Mui-selected': {
        color: theme.palette.text.primary, 
    },
  },
  tabPanel: {
    padding: theme.spacing(3), // Add padding around tab content
  },
}));

function App() {
  const classes = useStyles();
  const [tabValue, setTabValue] = useState(0);
  const [drawerOpen, setDrawerOpen] = useState(false);
  const [showLandingPage, setShowLandingPage] = useState(true);

  const handleTabChange = (event, newValue) => {
    setTabValue(newValue);
  };

  const toggleDrawer = () => {
    setDrawerOpen(!drawerOpen);
  };

  const navigateToApp = (initialTabIndex = 0) => {
    setTabValue(initialTabIndex);
    setShowLandingPage(false);
  };

  const menuItems = [
    { text: 'Dashboard', icon: <HomeIcon />, index: 0 },
    { text: 'Molecule Designer', icon: <BuildIcon />, index: 1 },
    { text: 'Simulations', icon: <TimelineIcon />, index: 2 },
    { text: 'Literature Explorer', icon: <DescriptionIcon />, index: 3 },
    { text: 'Comparison Tool', icon: <CompareArrowsIcon />, index: 4 },
    { text: 'Regulatory Analysis', icon: <TimelineIcon />, index: 5 },
  ];

  return (
    <ThemeProvider theme={landingTheme}>
      <CssBaseline />
      {showLandingPage ? (
        <LandingPage onNavigate={navigateToApp} />
      ) : (
        <div className={classes.root}>
          <AppBar position="fixed" className={classes.appBar}>
            <Toolbar>
              <IconButton
                edge="start"
                color="inherit"
                aria-label="menu"
                onClick={toggleDrawer}
              >
                <MenuIcon />
              </IconButton>
              <Typography variant="h6" className={classes.title}>
                Breaking Good | {menuItems[tabValue]?.text || 'Platform'}
              </Typography>
              <Button color="inherit" onClick={() => setShowLandingPage(true)}>Home</Button>
            </Toolbar>
          </AppBar>

          <Drawer
            className={classes.drawer}
            variant="temporary"
            open={drawerOpen}
            onClose={toggleDrawer}
            classes={{
              paper: classes.drawerPaper,
            }}
          >
            <div className={classes.toolbar} />
            <div className={classes.drawerContainer}>
              <List>
                {menuItems.map((item) => (
                  <ListItem 
                    button 
                    key={item.text}
                    selected={item.index === tabValue}
                    onClick={() => {
                      setTabValue(item.index);
                      setDrawerOpen(false);
                    }}
                  >
                    <ListItemIcon>{item.icon}</ListItemIcon>
                    <ListItemText primary={item.text} />
                  </ListItem>
                ))}
              </List>
              <Divider />
              <List>
                <ListItem button>
                  <ListItemIcon><DescriptionIcon /></ListItemIcon>
                  <ListItemText primary="Documentation" />
                </ListItem>
              </List>
            </div>
          </Drawer>

          <main className={classes.content}>
            <div className={classes.toolbar} />
            
            <Box className={classes.tabsContainer}>
              <Container>
                <Tabs
                  value={tabValue}
                  onChange={handleTabChange}
                  variant="scrollable"
                  scrollButtons="auto"
                  indicatorColor="primary"
                  textColor="primary"
                  aria-label="app navigation tabs"
                >
                  {menuItems.map((item) => (
                    <Tab 
                      key={item.text}
                      className={classes.tab}
                      icon={item.icon}
                      label={item.text}
                    />
                  ))}
                </Tabs>
              </Container>
            </Box>

            <Container className={classes.tabPanel}>
              <TabPanel value={tabValue} index={0}>
                <Typography variant="h4" gutterBottom>
                  Dashboard
                </Typography>
                <Typography paragraph>
                  Welcome to Breaking Good - an innovative platform for neuropharmacology drug discovery focused on designing safer alternatives to Adderall for ADHD treatment.
                </Typography>
                <Typography paragraph>
                  Navigate through the tabs to access different tools and functionalities of the platform. You can design molecules, run simulations, explore scientific literature, compare different molecules, and analyze regulatory pathways.
                </Typography>
              </TabPanel>
              
              <TabPanel value={tabValue} index={1}>
                <MoleculeDesigner />
              </TabPanel>
              
              <TabPanel value={tabValue} index={2}>
                <SimulationPanel />
              </TabPanel>
              
              <TabPanel value={tabValue} index={3}>
                <LiteratureExplorer />
              </TabPanel>
              
              <TabPanel value={tabValue} index={4}>
                <ComparisonTool />
              </TabPanel>
              
              <TabPanel value={tabValue} index={5}>
                <RegulatoryAnalysis />
              </TabPanel>
            </Container>
          </main>
        </div>
      )}
    </ThemeProvider>
  );
}

export default App; 