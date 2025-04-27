import React from 'react';
import { 
    Box, 
    Container, 
    Typography, 
    Button, 
    AppBar, 
    Toolbar, 
    IconButton, 
    useTheme, 
    useMediaQuery,
    makeStyles 
} from '@material-ui/core';
import MenuIcon from '@material-ui/icons/Menu';
import ArrowDownwardIcon from '@material-ui/icons/ArrowDownward'; // For scroll indicator

const useStyles = makeStyles((theme) => ({
    root: {
        display: 'flex',
        flexDirection: 'column',
        minHeight: '100vh',
        backgroundColor: theme.palette.background.default, // #FBF9F7
    },
    appBar: {
        backgroundColor: 'transparent', // As per spec
        boxShadow: 'none',
        borderBottom: `1px solid ${theme.palette.divider}`, // #EAE8E4
    },
    toolbar: {
        display: 'flex',
        justifyContent: 'space-between',
    },
    navLinks: {
        display: 'flex',
        gap: theme.spacing(4),
    },
    navLink: {
        fontFamily: 'ui-serif, Georgia, Cambria, "Times New Roman", Times, serif',
        fontSize: '1rem', // Updated size (1rem * 18px = 18px)
        fontWeight: 400,
        lineHeight: 1.4,
        color: theme.palette.text.primary,
        textDecoration: 'none',
        letterSpacing: '0em',
        position: 'relative', // For animated underline
        paddingBottom: theme.spacing(0.5),
        '&::after': { // Underline element
            content: '""',
            position: 'absolute',
            left: 0,
            bottom: 0,
            width: '0%',
            height: '1px',
            backgroundColor: theme.palette.primary.main,
            transition: 'width 0.2s ease-out',
        },
        '&:hover::after': {
            width: '100%', // Animate width on hover
        },
        '&:hover': {
            color: theme.palette.primary.main, 
            backgroundColor: 'transparent', // Prevent default hover background
        }
    },
    heroContainer: {
        flexGrow: 1,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        textAlign: 'center',
        paddingTop: theme.spacing(10), // ~80px
        paddingBottom: theme.spacing(10),
        paddingLeft: theme.spacing(3), // ~24px mobile
        paddingRight: theme.spacing(3),
        [theme.breakpoints.up('md')]: {
            paddingLeft: theme.spacing(20), // ~160px desktop
            paddingRight: theme.spacing(20),
        },
        // Apply fade-in animation
        animation: '$fadeIn 0.8s ease-out forwards',
        opacity: 0, // Start hidden
    },
    headlineBreaking: {
        color: theme.palette.text.primary, 
    },
    headlineGood: {
        color: theme.palette.primary.main, 
    },
    // Apply responsive H1 size using media queries
    headline: {
        fontSize: '3rem', // Mobile first (54px)
        [theme.breakpoints.up('sm')]: {
           fontSize: '4rem', // Tablet+ (72px)
        },
         [theme.breakpoints.up('lg')]: {
           fontSize: '4.5rem', // Large desktop (81px)
        },
    },
    tagline: {
        marginTop: theme.spacing(2), // 1rem relative to 18px base = 18px (Spec asked for 16px, using 1rem)
        fontFamily: 'Inter, Helvetica Neue, Arial, sans-serif',
        fontSize: '1rem', // Mobile first (1rem = 18px)
        lineHeight: 1.5,
        color: theme.palette.text.secondary, 
        maxWidth: '600px',
        marginLeft: 'auto',
        marginRight: 'auto',
        letterSpacing: '0.5px', // Subtle tracking
        // Responsive size
        [theme.breakpoints.up('sm')]: {
            fontSize: '1.125rem', // 20px on tablet+
        },
        [theme.breakpoints.up('lg')]: {
            fontSize: '1.375rem', // 22px on large desktop
        },
    },
    loginButton: {
        // Uses button overrides from theme
    },
    ctaButton: {
        marginTop: theme.spacing(4), // 2rem relative to 18px = 36px
        display: 'inline-block',
        padding: '12px 32px', // Already set in theme overrides for containedPrimary
    },
    scrollIndicator: {
        marginTop: theme.spacing(6), // 3rem relative to 18px = 54px
        fontSize: '2rem', // Approx 36px
        color: theme.palette.text.secondary, // Use warm gray
        animation: '$bounce 1.5s infinite alternate ease-in-out',
    },
    // Keyframes for animations
    '@keyframes fadeIn': {
        'from': { opacity: 0, transform: 'translateY(10px)' },
        'to': { opacity: 1, transform: 'translateY(0)' },
    },
    '@keyframes bounce': {
        'from': { transform: 'translateY(0)' },
        'to': { transform: 'translateY(8px)' },
    }
}));

const LandingPage = ({ onNavigate }) => {
    const classes = useStyles();
    const theme = useTheme();
    const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

    // Placeholder navigation items - replace with actual app sections
    const navItems = [
        { label: 'Dashboard', index: 0 },
        { label: 'Molecule Designer', index: 1 },
        { label: 'Simulations', index: 2 },
        { label: 'Literature', index: 3 },
        { label: 'Compare', index: 4 },
        { label: 'Regulatory', index: 5 },
    ];

    return (
        <div className={classes.root}>
            <AppBar position="static" className={classes.appBar} elevation={0}>
                <Container maxWidth="lg">
                    <Toolbar disableGutters className={classes.toolbar}>
                        {/* Left: Hamburger (Mobile) or Logo (Desktop) */}
                        <Box sx={{ flexGrow: { xs: 1, md: 0 } }}>
                            {isMobile ? (
                                <IconButton edge="start" color="inherit" aria-label="menu">
                                    <MenuIcon />
                                </IconButton>
                            ) : (
                                <Typography variant="h6" component="div" sx={{ fontFamily: 'Copernicus, ui-serif'}}> 
                                    {/* Replace with actual logo if available */}
                                    Anthropic / Claude
                                </Typography>
                            )}
                        </Box>

                        {/* Center: Nav Links (Desktop) */}
                        {!isMobile && (
                            <Box className={classes.navLinks} sx={{ flexGrow: 1, justifyContent: 'center' }}>
                                {navItems.map((item) => (
                                    <Button 
                                        key={item.label}
                                        className={classes.navLink}
                                        onClick={() => onNavigate(item.index)} // Use callback to switch view
                                    >
                                        {item.label}
                                    </Button>
                                ))}
                            </Box>
                        )}

                        {/* Right: Login Button */}
                        <Box sx={{ flexGrow: { xs: 1, md: 0 }, textAlign: 'right' }}>
                             {/* Using outlined as secondary style */} 
                            <Button variant="outlined" color="primary" className={classes.loginButton}>
                                Login
                            </Button>
                        </Box>
                    </Toolbar>
                </Container>
            </AppBar>

            <Container maxWidth="md" className={classes.heroContainer}>
                <Box>
                    {/* Apply responsive H1 style */}
                    <Typography variant="h1" className={classes.headline}>
                        <span className={classes.headlineBreaking}>Breaking</span>{' '}
                        <span className={classes.headlineGood}>Good</span>
                    </Typography>
                    <Typography variant="body1" className={classes.tagline}>
                        The future of drug-discovery research, analysis, and production. Powered by Claude.
                    </Typography>
                    
                    {/* Primary CTA Button */}
                     <Button 
                         variant="contained" 
                         color="primary" 
                         className={classes.ctaButton}
                         onClick={() => onNavigate(0)} // Navigate to Dashboard (index 0)
                     >
                         Get Started
                     </Button> 
                     
                     {/* Scroll Down Indicator */}
                     <div className={classes.scrollIndicator}>
                         <ArrowDownwardIcon fontSize="inherit" />
                     </div>
                </Box>
            </Container>
        </div>
    );
};

export default LandingPage; 