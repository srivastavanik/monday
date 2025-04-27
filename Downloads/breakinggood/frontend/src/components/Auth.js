import React, { useState, useEffect } from 'react';
import { 
  Paper, 
  Typography, 
  TextField, 
  Button, 
  Grid, 
  Link, 
  FormControlLabel, 
  Checkbox,
  Box,
  Tabs,
  Tab,
  makeStyles,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  CircularProgress,
  Snackbar
} from '@material-ui/core';
import Alert from '@material-ui/lab/Alert';

const useStyles = makeStyles((theme) => ({
  root: {
    padding: theme.spacing(3),
    maxWidth: 500,
    margin: '0 auto',
    marginTop: theme.spacing(10)
  },
  form: {
    width: '100%',
    marginTop: theme.spacing(1),
  },
  submit: {
    margin: theme.spacing(3, 0, 2),
  },
  tabPanel: {
    padding: theme.spacing(3, 0),
  },
  title: {
    marginBottom: theme.spacing(3),
    textAlign: 'center',
    fontWeight: 500
  },
  formControl: {
    marginBottom: theme.spacing(2),
    minWidth: '100%',
  },
}));

function TabPanel(props) {
  const { children, value, index, ...other } = props;

  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`auth-tabpanel-${index}`}
      aria-labelledby={`auth-tab-${index}`}
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

const Auth = ({ onLogin }) => {
  const classes = useStyles();
  const [tabValue, setTabValue] = useState(0);
  const [loading, setLoading] = useState(false);
  const [success, setSuccess] = useState(false);
  const [error, setError] = useState(null);
  
  // Login form state
  const [loginEmail, setLoginEmail] = useState('');
  const [loginPassword, setLoginPassword] = useState('');
  const [rememberMe, setRememberMe] = useState(false);
  
  // Register form state
  const [registerEmail, setRegisterEmail] = useState('');
  const [registerPassword, setRegisterPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [firstName, setFirstName] = useState('');
  const [lastName, setLastName] = useState('');
  const [organization, setOrganization] = useState('');
  const [role, setRole] = useState('');
  const [researchInterests, setResearchInterests] = useState('');
  
  // Check for stored user info
  useEffect(() => {
    const storedUser = localStorage.getItem('user');
    if (storedUser) {
      try {
        const user = JSON.parse(storedUser);
        if (user && user.email) {
          setLoginEmail(user.email);
          setRememberMe(true);
        }
      } catch (e) {
        console.error('Error parsing stored user data', e);
      }
    }
  }, []);
  
  const handleTabChange = (event, newValue) => {
    setTabValue(newValue);
    setError(null);
  };
  
  const handleLogin = (e) => {
    e.preventDefault();
    if (!loginEmail || !loginPassword) {
      setError('Please fill in all required fields');
      return;
    }
    
    setLoading(true);
    setError(null);
    
    // Simulate API call for login
    setTimeout(() => {
      setLoading(false);
      
      // For demonstration, accept any email with valid format and password length >= 6
      if (loginEmail.includes('@') && loginPassword.length >= 6) {
        // Store user info if remember me is checked
        if (rememberMe) {
          localStorage.setItem('user', JSON.stringify({ email: loginEmail }));
        } else {
          localStorage.removeItem('user');
        }
        
        // Create a user object with mock data
        const user = {
          id: `user-${Date.now()}`,
          email: loginEmail,
          firstName: 'Jane',
          lastName: 'Researcher',
          role: 'Principal Investigator',
          organization: 'University Research Lab',
          token: `token-${Date.now()}`
        };
        
        // Store current user in sessionStorage
        sessionStorage.setItem('currentUser', JSON.stringify(user));
        
        // Notify parent component about successful login
        if (onLogin) {
          onLogin(user);
        }
        
        setSuccess(true);
      } else {
        setError('Invalid email or password. Password must be at least 6 characters.');
      }
    }, 1000);
  };
  
  const handleRegister = (e) => {
    e.preventDefault();
    
    // Validation
    if (!registerEmail || !registerPassword || !confirmPassword || !firstName || !lastName || !organization || !role) {
      setError('Please fill in all required fields');
      return;
    }
    
    if (registerPassword !== confirmPassword) {
      setError('Passwords do not match');
      return;
    }
    
    if (registerPassword.length < 6) {
      setError('Password must be at least 6 characters');
      return;
    }
    
    if (!registerEmail.includes('@')) {
      setError('Please enter a valid email');
      return;
    }
    
    setLoading(true);
    setError(null);
    
    // Simulate API call for registration
    setTimeout(() => {
      setLoading(false);
      
      // Create a user object with registration data
      const user = {
        id: `user-${Date.now()}`,
        email: registerEmail,
        firstName,
        lastName,
        role,
        organization,
        researchInterests,
        token: `token-${Date.now()}`
      };
      
      // For demonstration, always succeed
      setSuccess(true);
      
      // Switch to login tab after successful registration
      setTabValue(0);
      setLoginEmail(registerEmail);
      
      // Clear register form
      setRegisterPassword('');
      setConfirmPassword('');
    }, 1000);
  };
  
  const roles = [
    'Principal Investigator',
    'Researcher',
    'Student',
    'Industry Professional',
    'Regulatory Affairs Specialist',
    'Other'
  ];
  
  const handleSnackbarClose = () => {
    setSuccess(false);
  };
  
  return (
    <Paper className={classes.root} elevation={3}>
      <Typography variant="h5" className={classes.title}>
        Breaking Good
      </Typography>
      
      <Tabs
        value={tabValue}
        onChange={handleTabChange}
        indicatorColor="primary"
        textColor="primary"
        variant="fullWidth"
      >
        <Tab label="Login" />
        <Tab label="Register" />
      </Tabs>
      
      <TabPanel value={tabValue} index={0} className={classes.tabPanel}>
        <form className={classes.form} onSubmit={handleLogin}>
          <TextField
            variant="outlined"
            margin="normal"
            required
            fullWidth
            id="email"
            label="Email Address"
            name="email"
            autoComplete="email"
            autoFocus
            value={loginEmail}
            onChange={(e) => setLoginEmail(e.target.value)}
            disabled={loading}
          />
          <TextField
            variant="outlined"
            margin="normal"
            required
            fullWidth
            name="password"
            label="Password"
            type="password"
            id="password"
            autoComplete="current-password"
            value={loginPassword}
            onChange={(e) => setLoginPassword(e.target.value)}
            disabled={loading}
          />
          <FormControlLabel
            control={
              <Checkbox 
                value="remember" 
                color="primary" 
                checked={rememberMe}
                onChange={(e) => setRememberMe(e.target.checked)}
                disabled={loading}
              />
            }
            label="Remember me"
          />
          
          {error && (
            <Alert severity="error" style={{ marginTop: 16 }}>
              {error}
            </Alert>
          )}
          
          <Button
            type="submit"
            fullWidth
            variant="contained"
            color="primary"
            className={classes.submit}
            disabled={loading}
          >
            {loading ? <CircularProgress size={24} /> : "Sign In"}
          </Button>
          
          <Grid container>
            <Grid item xs>
              <Link href="#" variant="body2">
                Forgot password?
              </Link>
            </Grid>
            <Grid item>
              <Link href="#" variant="body2" onClick={(e) => {
                e.preventDefault();
                setTabValue(1);
              }}>
                {"Don't have an account? Sign Up"}
              </Link>
            </Grid>
          </Grid>
        </form>
      </TabPanel>
      
      <TabPanel value={tabValue} index={1} className={classes.tabPanel}>
        <form className={classes.form} onSubmit={handleRegister}>
          <Grid container spacing={2}>
            <Grid item xs={12} sm={6}>
              <TextField
                autoComplete="fname"
                name="firstName"
                variant="outlined"
                required
                fullWidth
                id="firstName"
                label="First Name"
                autoFocus
                value={firstName}
                onChange={(e) => setFirstName(e.target.value)}
                disabled={loading}
              />
            </Grid>
            <Grid item xs={12} sm={6}>
              <TextField
                variant="outlined"
                required
                fullWidth
                id="lastName"
                label="Last Name"
                name="lastName"
                autoComplete="lname"
                value={lastName}
                onChange={(e) => setLastName(e.target.value)}
                disabled={loading}
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                variant="outlined"
                required
                fullWidth
                id="email"
                label="Email Address"
                name="email"
                autoComplete="email"
                value={registerEmail}
                onChange={(e) => setRegisterEmail(e.target.value)}
                disabled={loading}
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                variant="outlined"
                required
                fullWidth
                name="password"
                label="Password"
                type="password"
                id="password"
                autoComplete="new-password"
                value={registerPassword}
                onChange={(e) => setRegisterPassword(e.target.value)}
                disabled={loading}
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                variant="outlined"
                required
                fullWidth
                name="confirmPassword"
                label="Confirm Password"
                type="password"
                id="confirmPassword"
                value={confirmPassword}
                onChange={(e) => setConfirmPassword(e.target.value)}
                disabled={loading}
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                variant="outlined"
                required
                fullWidth
                name="organization"
                label="Organization"
                id="organization"
                value={organization}
                onChange={(e) => setOrganization(e.target.value)}
                disabled={loading}
              />
            </Grid>
            <Grid item xs={12}>
              <FormControl variant="outlined" className={classes.formControl}>
                <InputLabel id="role-label">Role</InputLabel>
                <Select
                  labelId="role-label"
                  id="role"
                  value={role}
                  onChange={(e) => setRole(e.target.value)}
                  label="Role"
                  required
                  disabled={loading}
                >
                  <MenuItem value=""><em>Select a role</em></MenuItem>
                  {roles.map((roleOption) => (
                    <MenuItem key={roleOption} value={roleOption}>
                      {roleOption}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            </Grid>
            <Grid item xs={12}>
              <TextField
                variant="outlined"
                fullWidth
                name="researchInterests"
                label="Research Interests"
                id="researchInterests"
                multiline
                rows={3}
                value={researchInterests}
                onChange={(e) => setResearchInterests(e.target.value)}
                helperText="Optional: Share your research interests to help us personalize your experience"
                disabled={loading}
              />
            </Grid>
          </Grid>
          
          {error && (
            <Alert severity="error" style={{ marginTop: 16 }}>
              {error}
            </Alert>
          )}
          
          <Button
            type="submit"
            fullWidth
            variant="contained"
            color="primary"
            className={classes.submit}
            disabled={loading}
          >
            {loading ? <CircularProgress size={24} /> : "Sign Up"}
          </Button>
          
          <Grid container justifyContent="flex-end">
            <Grid item>
              <Link href="#" variant="body2" onClick={(e) => {
                e.preventDefault();
                setTabValue(0);
              }}>
                Already have an account? Sign in
              </Link>
            </Grid>
          </Grid>
        </form>
      </TabPanel>
      
      <Snackbar 
        open={success} 
        autoHideDuration={6000} 
        onClose={handleSnackbarClose}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
      >
        <Alert onClose={handleSnackbarClose} severity="success">
          {tabValue === 0 ? "Login successful!" : "Registration successful! You can now log in."}
        </Alert>
      </Snackbar>
    </Paper>
  );
};

export default Auth; 