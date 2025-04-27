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
  TextField,
  IconButton,
  Grid,
  Avatar,
  Container,
  makeStyles,
  Collapse,
  useTheme
} from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import BuildIcon from '@material-ui/icons/Build';
import FormatQuoteIcon from '@material-ui/icons/FormatQuote';
import LocalPharmacyIcon from '@material-ui/icons/LocalPharmacy';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import WarningIcon from '@material-ui/icons/Warning';
import SendIcon from '@material-ui/icons/Send';
import ScienceIcon from '@material-ui/icons/EmojiObjects'; // Using EmojiObjects instead of Science
import ChatIcon from '@material-ui/icons/Chat';
import { Alert, AlertTitle } from '@material-ui/lab';
import MoleculeViewer3D from './MoleculeViewer3DImproved';
import ReactMarkdown from 'react-markdown';
import { claudeAPI } from '../services/api';

const useStyles = makeStyles((theme) => ({
  root: {
    width: '100%',
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2),
    display: 'flex',
    flexDirection: 'column',
    height: '100%'
  },
  header: {
    display: 'flex',
    alignItems: 'center',
    marginBottom: theme.spacing(2),
  },
  headerTitle: {
    fontWeight: 500
  },
  messageContainer: {
    flexGrow: 1,
    overflowY: 'auto',
    padding: theme.spacing(2, 0),
    marginBottom: theme.spacing(2)
  },
  messageRow: {
    display: 'flex',
    marginBottom: theme.spacing(3),
    alignItems: 'flex-start'
  },
  avatar: {
    backgroundColor: theme.palette.primary.main,
    marginRight: theme.spacing(2),
    width: 40,
    height: 40
  },
  userAvatar: {
    backgroundColor: theme.palette.secondary.main
  },
  messageContent: {
    flexGrow: 1,
    backgroundColor: theme.palette.background.paper,
    borderRadius: theme.shape.borderRadius,
    padding: theme.spacing(2),
    boxShadow: theme.shadows[1],
  },
  userMessageContent: {
    backgroundColor: theme.palette.primary.light,
    color: theme.palette.primary.contrastText,
  },
  inputContainer: {
    display: 'flex',
    alignItems: 'center',
    borderTop: `1px solid ${theme.palette.divider}`,
    borderRadius: theme.shape.borderRadius,
    padding: theme.spacing(1),
    marginTop: 'auto'
  },
  input: {
    flexGrow: 1,
    borderRadius: theme.shape.borderRadius,
  },
  sendButton: {
    marginLeft: theme.spacing(1)
  },
  loadingContainer: {
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center',
    padding: theme.spacing(4),
  },
  sectionTitle: {
    marginTop: theme.spacing(3),
    marginBottom: theme.spacing(1),
    fontWeight: 500
  },
  accordionRoot: {
    marginBottom: theme.spacing(1),
    border: `1px solid ${theme.palette.divider}`,
    boxShadow: 'none',
    '&:before': {
      display: 'none',
    },
    '&.Mui-expanded': {
      margin: theme.spacing(1, 0),
    },
  },
  accordionSummary: {
    backgroundColor: 'rgba(0, 0, 0, .03)',
    borderBottom: '1px solid rgba(0, 0, 0, .125)',
    marginBottom: -1,
    minHeight: 56,
    '&.Mui-expanded': {
      minHeight: 56,
    },
  },
  accordionDetails: {
    flexDirection: 'column',
    padding: theme.spacing(2),
  },
  statusChip: {
    marginLeft: theme.spacing(1),
  },
  codeBlock: {
    fontFamily: 'monospace',
    backgroundColor: theme.palette.grey[200],
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    overflowX: 'auto',
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1),
    color: theme.palette.text.primary
  },
  moleculesGrid: {
    display: 'flex',
    flexWrap: 'wrap',
    gap: theme.spacing(2),
    marginTop: theme.spacing(2)
  },
  moleculeCard: {
    width: 280,
    margin: theme.spacing(1),
    overflow: 'hidden'
  },
  viewerContainer: {
    height: 200,
    width: '100%',
    backgroundColor: '#fff',
    borderRadius: `${theme.shape.borderRadius}px ${theme.shape.borderRadius}px 0 0`
  },
  cardContent: {
    padding: theme.spacing(2),
  },
  moleculeName: {
    fontSize: '0.875rem',
    fontWeight: 500,
    marginBottom: theme.spacing(1)
  },
  smilesCode: {
    fontSize: '0.75rem',
    wordBreak: 'break-all',
    marginBottom: theme.spacing(1),
    fontFamily: 'monospace',
    backgroundColor: theme.palette.grey[200],
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius
  },
  cardActions: {
    justifyContent: 'space-between',
    padding: theme.spacing(1, 2)
  },
  markdownContent: {
    '& h1, & h2, & h3, & h4, & h5, & h6': {
      marginTop: theme.spacing(2),
      marginBottom: theme.spacing(1),
    },
    '& a': {
      color: theme.palette.primary.main
    },
    '& ul, & ol': {
      paddingLeft: theme.spacing(3)
    },
    '& pre': {
      backgroundColor: theme.palette.grey[200],
      padding: theme.spacing(1),
      borderRadius: theme.shape.borderRadius,
      overflowX: 'auto',
      fontFamily: 'monospace',
      fontSize: '0.875rem'
    },
    '& code': {
      fontFamily: 'monospace',
      backgroundColor: theme.palette.grey[100],
      padding: theme.spacing(0.2, 0.5),
      borderRadius: 3,
    },
    '& table': {
      borderCollapse: 'collapse',
      width: '100%',
      marginBottom: theme.spacing(2)
    },
    '& th, & td': {
      border: `1px solid ${theme.palette.divider}`,
      padding: theme.spacing(1)
    },
    '& th': {
      backgroundColor: theme.palette.grey[100],
    }
  },
  toggleButton: {
    marginTop: theme.spacing(1),
  },
}));

// This function parses the molecule data from the thinking process
const extractMoleculesFromThinking = (thinking) => {
  if (!thinking) return [];
  
  const molecules = [];
  
  // Find SMILES strings in the thinking process
  const smilesRegex = /SMILES:\s*([^\s\n]+)/g;
  let match;
  
  while ((match = smilesRegex.exec(thinking)) !== null) {
    const smiles = match[1].replace(/[`'"]/g, '');
    
    // Try to find a name for this molecule near the SMILES string
    const nameRegex = /\*\*Chemical Name:\*\*\s*([^\n]+)/;
    const nameMatch = thinking.substring(Math.max(0, match.index - 200), match.index + 300).match(nameRegex);
    
    molecules.push({
      smiles,
      name: nameMatch ? nameMatch[1].trim() : `Molecule ${molecules.length + 1}`
    });
  }
  
  return molecules;
};

// Extract sections from the thinking process for structured display
const extractSections = (thinking) => {
  if (!thinking) return [];
  
  const sections = [];
  const lines = thinking.split('\n');
  let currentSection = null;
  let currentContent = [];
  
  lines.forEach(line => {
    const trimmedLine = line.trim();
    // Check if this line is a header (using # markdown syntax)
    const isHeader = /^#{1,3}\s+(.+)$/.test(trimmedLine);
    
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

const ThinkingProcessImproved = ({ requestId, onSelectMolecule, onSaveMolecule }) => {
  const classes = useStyles();
  const theme = useTheme();
  const [thinkingData, setThinkingData] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [sections, setSections] = useState([]);
  const [expandedSection, setExpandedSection] = useState(null);
  const [molecules, setMolecules] = useState([]);
  const [chatMessages, setChatMessages] = useState([]);
  const [userInput, setUserInput] = useState('');
  const [sendingMessage, setSendingMessage] = useState(false);
  const messageContainerRef = useRef(null);
  const [isExpanded, setIsExpanded] = useState(false);
  
  // Fetch thinking process data when requestId changes
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
        
        // Extract molecules from the thinking process
        const extractedMolecules = extractMoleculesFromThinking(response.data.thinking);
        setMolecules(extractedMolecules);
        
        // Parse the thinking process into sections
        const parsedSections = extractSections(response.data.thinking);
        setSections(parsedSections);
        
        // Initialize chat messages with the thinking process
        setChatMessages([
          {
            id: 'system-thinking',
            type: 'claude',
            content: response.data.thinking,
            timestamp: response.data.timestamp
          }
        ]);
        
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
  
  // Scroll to bottom of messages when new messages are added
  useEffect(() => {
    if (messageContainerRef.current) {
      messageContainerRef.current.scrollTop = messageContainerRef.current.scrollHeight;
    }
  }, [chatMessages]);
  
  // Handle sending user message to continue the conversation
  const handleSendMessage = async () => {
    if (!userInput.trim() || sendingMessage) return;
    
    // Add user message to chat
    const userMessage = {
      id: `user-${Date.now()}`,
      type: 'user',
      content: userInput,
      timestamp: new Date().toISOString()
    };
    
    setChatMessages(prev => [...prev, userMessage]);
    setUserInput('');
    setSendingMessage(true);
    
    try {
      // Send user message to Claude API
      const response = await claudeAPI.askQuestion(userInput, { requestId, context: 'molecule-design' });
      
      // Add Claude's response to chat
      const claudeMessage = {
        id: `claude-${Date.now()}`,
        type: 'claude',
        content: response.data.response,
        timestamp: new Date().toISOString()
      };
      
      setChatMessages(prev => [...prev, claudeMessage]);
      
      // Check if the response contains new molecules
      const newMolecules = extractMoleculesFromThinking(response.data.response);
      if (newMolecules.length > 0) {
        setMolecules(prev => [...prev, ...newMolecules]);
      }
      
      setSendingMessage(false);
    } catch (err) {
      console.error('Error sending message:', err);
      
      // Add error message to chat
      const errorMessage = {
        id: `error-${Date.now()}`,
        type: 'system',
        content: 'Sorry, there was an error processing your request. Please try again.',
        timestamp: new Date().toISOString()
      };
      
      setChatMessages(prev => [...prev, errorMessage]);
      setSendingMessage(false);
    }
  };
  
  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSendMessage();
    }
  };
  
  const handleAccordionChange = (index) => (event, isExpanded) => {
    setExpandedSection(isExpanded ? index : null);
  };
  
  const handleSelectMolecule = (molecule) => {
    if (onSelectMolecule) {
      onSelectMolecule(molecule.smiles);
    }
  };
  
  const handleSaveMolecule = (molecule) => {
    if (onSaveMolecule) {
      onSaveMolecule(molecule);
    }
  };
  
  // Format timestamp
  const formatDate = (isoString) => {
    if (!isoString) return '';
    const date = new Date(isoString);
    return date.toLocaleString();
  };
  
  const handleToggleExpand = () => {
    setIsExpanded(!isExpanded);
  };
  
  if (loading) {
    return (
      <Paper className={classes.claudeContainer}>
        <div className={classes.loadingContainer}>
          <CircularProgress color="primary" />
          <Typography variant="body1" style={{ marginTop: 16, color: '#fff' }}>
            Loading thinking process...
          </Typography>
        </div>
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
  
  // Extract the main thinking text
  let thinkingText = 'No thinking process text available.';
  if (thinkingData && thinkingData.thinking && Array.isArray(thinkingData.thinking.content) && thinkingData.thinking.content.length > 0) {
    thinkingText = thinkingData.thinking.content[0]?.text || 'No text content found.';
  } else if (typeof thinkingData?.thinking === 'string') {
    // Fallback if thinking is just a string (older format?)
    thinkingText = thinkingData.thinking;
  }
  
  const truncatedText = thinkingText.substring(0, 1000); // Show first 1000 characters
  const needsTruncation = thinkingText.length > 1000;
  
  return (
    <div className={classes.root}>
      <div className={classes.header}>
        <Avatar className={classes.avatar}>
          <ScienceIcon />
        </Avatar>
        <Typography variant="h6" className={classes.headerTitle}>
          AI Thinking Process
          {thinkingData.status && (
            <Chip 
              label={thinkingData.status} 
              color="primary"
              size="small"
              icon={thinkingData.status === 'completed' ? <CheckCircleIcon /> : <CircularProgress size={16} />}
              className={classes.statusChip}
            />
          )}
        </Typography>
      </div>
      
      <Typography variant="body2" style={{ color: 'rgba(255,255,255,0.7)', marginBottom: 16 }}>
        Generated: {formatDate(thinkingData.timestamp)}
      </Typography>
      
      {/* Molecule Display Section */}
      {molecules.length > 0 && (
        <>
          <Typography variant="h6" className={classes.sectionTitle}>
            Generated Molecules
          </Typography>
          
          <div className={classes.moleculesGrid}>
            {molecules.map((molecule, index) => (
              <Card key={index} className={classes.moleculeCard}>
                <div className={classes.viewerContainer}>
                  <MoleculeViewer3D 
                    moleculeData={molecule.smiles}
                    format="smiles"
                    height={200}
                    width="100%"
                  />
                </div>
                <CardContent className={classes.cardContent}>
                  <Typography className={classes.moleculeName}>
                    {molecule.name}
                  </Typography>
                  <div className={classes.smilesCode}>
                    {molecule.smiles}
                  </div>
                  <Box display="flex" justifyContent="space-between" mt={1}>
                    <Button 
                      size="small" 
                      variant="outlined" 
                      color="primary"
                      startIcon={<LocalPharmacyIcon />}
                      onClick={() => handleSelectMolecule(molecule)}
                    >
                      Select
                    </Button>
                    <Button 
                      size="small" 
                      variant="contained" 
                      color="primary"
                      onClick={() => handleSaveMolecule(molecule)}
                    >
                      Save
                    </Button>
                  </Box>
                </CardContent>
              </Card>
            ))}
          </div>
        </>
      )}
      
      {/* Chat Messages */}
      <Typography variant="h6" className={classes.sectionTitle}>
        Detailed Analysis
      </Typography>
      
      <div className={classes.messageContainer} ref={messageContainerRef}>
        {chatMessages.map((message) => (
          <div key={message.id} className={classes.messageRow}>
            <Avatar className={`${classes.avatar} ${message.type === 'user' ? classes.userAvatar : ''}`}>
              {message.type === 'user' ? 'U' : message.type === 'claude' ? <ScienceIcon fontSize="small"/> : '!'}
            </Avatar>
            <div className={`${classes.messageContent} ${message.type === 'user' ? classes.userMessageContent : ''}`}>
              <div className={classes.markdownContent}>
                <ReactMarkdown>
                  {message.content}
                </ReactMarkdown>
              </div>
              <Typography variant="caption" style={{ color: theme.palette.text.secondary, marginTop: 8, display: 'block' }}>
                {formatDate(message.timestamp)}
              </Typography>
            </div>
          </div>
        ))}
        
        {sendingMessage && (
          <div className={classes.messageRow}>
            <Avatar className={classes.avatar}>
              <ScienceIcon fontSize="small"/>
            </Avatar>
            <div className={classes.messageContent} style={{ display: 'flex', alignItems: 'center' }}>
              <CircularProgress size={20} style={{ marginRight: 8 }} />
              <Typography>Claude is thinking...</Typography>
            </div>
          </div>
        )}
      </div>
      
      {/* Input Field for Continued Conversation */}
      <div className={classes.inputContainer}>
        <TextField
          className={classes.input}
          variant="outlined"
          placeholder="Ask a follow-up question or request refinements..."
          fullWidth
          multiline
          rows={2}
          value={userInput}
          onChange={(e) => setUserInput(e.target.value)}
          onKeyPress={handleKeyPress}
          disabled={sendingMessage}
        />
        <IconButton 
          className={classes.sendButton} 
          onClick={handleSendMessage}
          disabled={!userInput.trim() || sendingMessage}
          color="primary"
        >
          <SendIcon />
        </IconButton>
      </div>
      
      {/* Structured Thinking Process (Accordions) */}
      <Typography variant="h6" className={classes.sectionTitle} style={{ marginTop: 16 }}>
        Raw Thinking Output (Sections)
      </Typography>
      
      {sections.length > 0 ? (
        sections.map((section, index) => (
          <Accordion
            key={index}
            expanded={expandedSection === index}
            onChange={handleAccordionChange(index)}
            className={classes.accordionRoot}
            variant="outlined"
          >
            <AccordionSummary 
              expandIcon={<ExpandMoreIcon />}
              className={classes.accordionSummary}
            >
              <Typography>{section.title.replace(/^[#]+\s*/, '')}</Typography>
            </AccordionSummary>
            <AccordionDetails className={classes.accordionDetails}>
              <div className={classes.markdownContent}>
                <Collapse in={expandedSection === index} collapsedSize="150px">
                <ReactMarkdown>
                  {section.content}
                </ReactMarkdown>
                </Collapse>
              </div>
            </AccordionDetails>
          </Accordion>
        ))
      ) : (
        <Typography variant="body2" color="textSecondary">
          No structured sections found in the thinking process.
        </Typography>
      )}
    </div>
  );
};

export default ThinkingProcessImproved;
