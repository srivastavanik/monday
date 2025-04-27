import React, { useState, useEffect, useRef } from 'react';
import {
  Paper,
  TextField,
  Button,
  List,
  ListItem,
  Typography,
  makeStyles,
  CircularProgress,
  Box,
  IconButton,
  Avatar,
  Container, // Added for potential centering
  useTheme,
  InputAdornment
} from '@material-ui/core';
import SendIcon from '@material-ui/icons/Send';
import PersonIcon from '@material-ui/icons/Person';
import AssistantIcon from '@material-ui/icons/EmojiObjects'; // Using EmojiObjects for AI
import Alert from '@material-ui/lab/Alert';
import { claudeAPI } from '../services/api';

// TODO: Define styles based on specifications
const useStyles = makeStyles((theme) => ({
  // --- Root Container --- 
  chatPanel: {
    background: '#FCF9F2', // Warm off-white
    display: 'flex',
    flexDirection: 'column',
    height: 'calc(100vh - 180px)', // Example height, adjust as needed for sticky positioning
    maxHeight: '700px', // Max height before internal scroll
    width: '100%',
    maxWidth: '640px', // Centered, focused container
    margin: '0 auto', // Center horizontally
    borderRadius: '16px',
    boxShadow: '0px 4px 12px rgba(0,0,0,0.08)', // Softer shadow
    overflow: 'hidden',
    position: 'relative', // Needed for absolute positioned elements if any
  },
  stickyContainer: { // Apply this class to the PARENT Grid item in MoleculeDesigner
      position: 'sticky',
      top: '80px', // Adjust based on your header height
      alignSelf: 'flex-start', // Prevent stretching if parent is flex
      height: 'fit-content', // Important for sticky
  },
  // --- Message List --- 
  messageList: {
    flexGrow: 1,
    overflowY: 'auto',
    padding: theme.spacing(2, 3), // More horizontal padding
    '&::-webkit-scrollbar': {
        width: '6px',
    },
    '&::-webkit-scrollbar-thumb': {
        backgroundColor: 'rgba(0,0,0,0.2)',
        borderRadius: '3px',
    },
    '&::-webkit-scrollbar-track': {
        backgroundColor: 'transparent',
    },
  },
  // --- Message Bubbles & Items --- 
  messageItem: {
    display: 'flex',
    marginBottom: theme.spacing(2),
    padding: 0, // Remove default padding
  },
  userMessageItem: {
    justifyContent: 'flex-end',
  },
  assistantMessageItem: {
    justifyContent: 'flex-start',
  },
  messageBubble: {
    padding: theme.spacing(1.5, 2), // Generous padding
    borderRadius: '12px', // Rounded geometry
    maxWidth: '85%',
    lineHeight: 1.5,
    boxShadow: '0px 2px 4px rgba(0,0,0,0.05)', // Subtle shadow
    wordWrap: 'break-word',
  },
  userBubble: {
    backgroundColor: '#DF8A65', // Soft terracotta accent
    color: '#FFFFFF', // White text
    marginLeft: 'auto',
  },
  assistantBubble: {
    backgroundColor: '#FFFFFF', // Surface card white
    color: '#2B2B2B', // Dark charcoal text
    marginRight: 'auto',
    border: `1px solid rgba(0,0,0,0.08)`,
  },
  errorBubble: {
    backgroundColor: theme.palette.error.light,
    color: theme.palette.error.contrastText,
    marginRight: 'auto',
    border: `1px solid ${theme.palette.error.dark}`,
  },
  // --- Avatars --- 
  avatar: {
    width: 36,
    height: 36,
    boxShadow: '0px 1px 3px rgba(0,0,0,0.1)',
  },
  userAvatar: {
    backgroundColor: '#DF8A65', // Accent color
    color: '#FFFFFF',
    marginLeft: theme.spacing(1.5),
    order: 2, // Place avatar after bubble for user
  },
  assistantAvatar: {
    backgroundColor: theme.palette.grey[200],
    color: theme.palette.grey[600],
    marginRight: theme.spacing(1.5),
  },
  // --- Input Area --- 
  messageInputContainer: {
    display: 'flex',
    padding: theme.spacing(1.5, 2), // Spacing
    borderTop: `1px solid rgba(0,0,0,0.08)`,
    backgroundColor: 'rgba(255, 255, 255, 0.7)', // Slight transparency
    backdropFilter: 'blur(5px)', // Optional blur effect
    alignItems: 'center',
  },
  inputField: {
    flexGrow: 1,
    marginRight: theme.spacing(1),
    '& .MuiOutlinedInput-root': {
      borderRadius: '12px', // Rounded field
      backgroundColor: '#FFFFFF', // White input background
      fontSize: '16px',
      // Adjust padding if needed
    },
    '& .MuiOutlinedInput-notchedOutline': {
        borderColor: 'rgba(0, 0, 0, 0.1)', // Lighter border
    },
    '& .MuiOutlinedInput-root.Mui-focused .MuiOutlinedInput-notchedOutline': {
        borderColor: '#DF8A65', // Accent focus border
    },
  },
  sendButton: {
    backgroundColor: '#DF8A65',
    color: '#FFFFFF',
    borderRadius: '12px',
    padding: theme.spacing(1.2), // Adjust padding for size
    boxShadow: '0px 2px 5px rgba(223, 138, 101, 0.3)',
    '&:hover': {
        backgroundColor: '#ca7a58', // Darker terracotta
    },
  },
  // --- Typography --- 
  messageText: {
      fontFamily: 'Inter, Roboto, sans-serif', // Sans-serif body
      fontSize: '16px',
      fontWeight: 400,
      lineHeight: 1.5, // Readability
      color: 'inherit', // Inherit color from bubble
  },
  errorText: {
      fontFamily: 'Inter, Roboto, sans-serif',
      fontSize: '14px',
      fontWeight: 500,
  },
  greetingHeader: { // Style for the potential greeting/status area
      fontFamily: 'Merriweather, serif', // Serif heading
      fontSize: '20px',
      fontWeight: 500,
      color: '#2B2B2B',
      padding: theme.spacing(2, 3, 1, 3),
      borderBottom: `1px solid rgba(0,0,0,0.08)`,
  },
  contextText: { // Small text for context messages
      fontSize: '12px', 
      color: theme.palette.text.secondary,
      textAlign: 'center',
      padding: theme.spacing(0.5, 0)
  }
}));

const AIChatInterfaceRedesigned = ({ initialContext = 'General drug design', onMoleculeMentioned = () => {}, selectedMolecule = null }) => {
  const classes = useStyles();
  const theme = useTheme(); // Get theme object
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const messagesEndRef = useRef(null);
  const [currentContext, setCurrentContext] = useState(initialContext);

  // Scroll to bottom effect
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  // Context update effect
  useEffect(() => {
    let newContext = initialContext;
    if (selectedMolecule) {
      newContext = `Discussing ${selectedMolecule.name || 'selected molecule'} (SMILES: ${selectedMolecule.smiles || 'N/A'})`;
    }
    // Only update context state, don't add visible system messages for now
    if (newContext !== currentContext) {
        setCurrentContext(newContext);
    }
  }, [initialContext, selectedMolecule, currentContext]);

  // Initial greeting message effect
  useEffect(() => {
    setMessages([
      {
        role: 'assistant',
        content: [{ type: 'text', text: `Hello! How can I help with drug design today?` }]
      }
    ]);
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [initialContext]);

  const handleInputChange = (event) => {
    setInput(event.target.value);
  };

  const handleSendMessage = async (event) => {
     event.preventDefault();
     const userMessageContent = input.trim();
     if (!userMessageContent || isLoading) return;

     const newUserMessage = { role: 'user', content: [{ type: 'text', text: userMessageContent }] };
     
     // Include current visible context message + new user message for API call
     let apiMessages = messages
         .filter(msg => msg.role !== 'system' && !msg.isError) // Exclude system/error messages
         .map(msg => ({ role: msg.role, content: msg.content })); // Map to API format
     apiMessages.push(newUserMessage); // Add the new user message

     // Add implicit focus context if a molecule is selected
     if (selectedMolecule) {
         const contextMessageForApi = {
             role: 'user', // Prepend as user instruction
             content: [{ type: 'text', text: `[Focus on Molecule: ${selectedMolecule.name}, SMILES: ${selectedMolecule.smiles}]` }]
         };
         // Insert context message before the last user message in the API payload
         apiMessages.splice(apiMessages.length - 1, 0, contextMessageForApi);
     }

     setMessages(prev => [...prev, newUserMessage]); // Add user message to UI
     setInput('');
     setIsLoading(true);
     setError(null);

     try {
         const response = await claudeAPI.continueChat(apiMessages); // Send curated messages

         if (response.data && response.data.content && Array.isArray(response.data.content)) {
             const assistantMessage = {
                 role: 'assistant',
                 content: response.data.content
             };
             setMessages(prevMessages => [...prevMessages, assistantMessage]);

             // SMILES extraction logic
             const responseText = assistantMessage.content[0]?.type === 'text' ? assistantMessage.content[0].text : '';
             const smilesRegex = /\b([A-Za-z0-9@+\-\[\]\(\)\\\/%=#$!.~{},*]+)\b/g;
             const potentialSmiles = responseText.match(smilesRegex);
             if (potentialSmiles) {
                 potentialSmiles.forEach(smiles => {
                     if (smiles.length > 5 && smiles.includes('C') && (smiles.includes('(') || smiles.includes('='))) {
                         onMoleculeMentioned(smiles);
                     }
                 });
             }
         } else {
             console.error('Invalid response structure from /chat API:', response.data);
             throw new Error('Invalid response format from API');
         }
     } catch (err) {
         console.error('Error sending chat message:', err);
         const errorText = err.response?.data?.details?.error?.message || err.response?.data?.error || err.message || 'Failed to get response from AI';
         setError(errorText);
         // Add error message bubble to the chat
         setMessages(prevMessages => [...prevMessages, { role: 'assistant', content: [{ type: 'text', text: `Error: ${errorText}` }], isError: true }]);
     } finally {
         setIsLoading(false);
     }
  };

  return (
    // --- Redesigned Structure --- 
    <Paper className={classes.chatPanel} elevation={3}>
        {/* Optional Header/Greeting Area */}
         <Box className={classes.greetingHeader}>
             {/* Example: Dynamic greeting could go here */}
             <Typography variant="h6" component="h2" style={{ fontFamily: 'Merriweather, serif' }}>
                 Claude Assistant
             </Typography>
         </Box>
         
        <List className={classes.messageList}>
            {messages.map((msg, index) => (
                <ListItem 
                    key={index} 
                    className={`${classes.messageItem} ${msg.role === 'user' ? classes.userMessageItem : classes.assistantMessageItem}`}
                >
                    {msg.role === 'assistant' && (
                        <Avatar className={`${classes.avatar} ${classes.assistantAvatar}`}>
                            <AssistantIcon fontSize="small"/>
                        </Avatar>
                    )}
                    <Box 
                        className={`${classes.messageBubble} ${msg.role === 'user' ? classes.userBubble : msg.isError ? classes.errorBubble : classes.assistantBubble}`}
                    >
                        {Array.isArray(msg.content) ? 
                        msg.content.map((part, partIndex) => (
                            part.type === 'text' ? 
                            <Typography key={partIndex} variant="body1" className={classes.messageText}>
                                {part.text}
                            </Typography> : null
                        ))
                        : <Typography variant="body1" className={classes.messageText}>{typeof msg.content === 'string' ? msg.content : 'Invalid format'}</Typography>
                        }
                    </Box>
                    {msg.role === 'user' && (
                        <Avatar className={`${classes.avatar} ${classes.userAvatar}`}>
                            <PersonIcon fontSize="small"/>
                        </Avatar>
                    )}
                </ListItem>
            ))}
            {isLoading && (
                <ListItem className={`${classes.messageItem} ${classes.assistantMessageItem}`}>
                    <Avatar className={`${classes.avatar} ${classes.assistantAvatar}`}>
                        <AssistantIcon fontSize="small"/>
                    </Avatar>
                    <Box className={`${classes.messageBubble} ${classes.assistantBubble}`}>
                        <CircularProgress size={20} color="inherit" />
                    </Box>
                </ListItem>
            )}
            {/* Error Alert at the bottom - might be redundant with error bubbles */}
            {/* error && (
                 <ListItem>
                     <Alert severity="error" style={{ width: '100%', borderRadius: '8px' }}>{error}</Alert>
                 </ListItem>
             ) */}
            <div ref={messagesEndRef} /> { /* Scroll anchor */ }
        </List>

        <form onSubmit={handleSendMessage} className={classes.messageInputContainer}>
            <TextField
                className={classes.inputField}
                variant="outlined"
                size="small"
                placeholder={selectedMolecule ? `Ask Claude about ${selectedMolecule.name}...` : 'Ask Claude anything...'}
                value={input}
                onChange={handleInputChange}
                disabled={isLoading}
                multiline
                maxRows={4} // Allow multi-line input
                InputProps={{
                    style: { fontFamily: 'Inter, Roboto, sans-serif', fontSize: '16px' },
                }}
            />
            <IconButton type="submit" className={classes.sendButton} disabled={isLoading || !input.trim()} aria-label="Send message">
                <SendIcon fontSize="small"/>
            </IconButton>
        </form>
    </Paper>
  );
};

export default AIChatInterfaceRedesigned; 