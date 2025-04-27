import React, { useState, useEffect, useRef } from 'react';
import {
  Paper,
  TextField,
  Button,
  List,
  ListItem,
  ListItemText,
  Typography,
  makeStyles,
  CircularProgress,
  Box,
  IconButton,
  Avatar
} from '@material-ui/core';
import SendIcon from '@material-ui/icons/Send';
import PersonIcon from '@material-ui/icons/Person';
import AssistantIcon from '@material-ui/icons/EmojiObjects'; // Using EmojiObjects for AI
import Alert from '@material-ui/lab/Alert';
import { claudeAPI } from '../services/api';

const useStyles = makeStyles((theme) => ({
  root: {
    display: 'flex',
    flexDirection: 'column',
    height: '500px', // Or adjust as needed
    border: `1px solid ${theme.palette.divider}`,
    borderRadius: theme.shape.borderRadius,
    backgroundColor: theme.palette.background.paper, // Use theme paper background
  },
  messageList: {
    flexGrow: 1,
    overflowY: 'auto',
    padding: theme.spacing(2),
  },
  messageInputContainer: {
    display: 'flex',
    padding: theme.spacing(1),
    borderTop: `1px solid ${theme.palette.divider}`,
    alignItems: 'center',
  },
  inputField: {
    flexGrow: 1,
    marginRight: theme.spacing(1),
  },
  userMessage: {
    display: 'flex',
    justifyContent: 'flex-end',
    marginBottom: theme.spacing(1),
  },
  assistantMessage: {
    display: 'flex',
    justifyContent: 'flex-start',
    marginBottom: theme.spacing(1),
  },
  messageBubble: {
    padding: theme.spacing(1, 2),
    borderRadius: theme.shape.borderRadius,
    maxWidth: '80%',
    wordWrap: 'break-word',
    boxShadow: theme.shadows[1], // Use subtle shadow from theme
    lineHeight: 1.6,
  },
  userBubble: {
    backgroundColor: theme.palette.primary.dark, // Richer terracotta #BD5D3A
    color: theme.palette.primary.contrastText, // White text
    marginLeft: 'auto',
  },
  assistantBubble: {
    backgroundColor: theme.palette.background.paper, // White paper
    color: theme.palette.text.primary, // Dark charcoal text
    marginRight: 'auto',
    border: `1px solid ${theme.palette.divider}` // Subtle divider border
  },
  errorBubble: {
    backgroundColor: theme.palette.error.light,
    color: theme.palette.error.contrastText,
    marginRight: 'auto',
    border: `1px solid ${theme.palette.error.main}`,
  },
  avatar: {
    backgroundColor: theme.palette.grey[300], // Use a light grey for assistant avatar
    color: theme.palette.text.secondary,
    marginRight: theme.spacing(1),
    marginLeft: theme.spacing(1),
    width: 32, 
    height: 32, 
  },
  userAvatar: {
    backgroundColor: theme.palette.primary.main, // Terracotta for user avatar
    color: theme.palette.primary.contrastText,
    // Inherits size from avatar class
  }
}));

const AIChatInterface = ({ initialContext = 'General drug design', onMoleculeMentioned = () => {}, selectedMolecule = null }) => {
  const classes = useStyles();
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const messagesEndRef = useRef(null);
  const [currentContext, setCurrentContext] = useState(initialContext);

  // Scroll to bottom when messages update
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  // Update context when initialContext or selectedMolecule changes
  useEffect(() => {
    let newContext = initialContext;
    if (selectedMolecule) {
      newContext = `Discussing ${selectedMolecule.name || 'selected molecule'} (SMILES: ${selectedMolecule.smiles || 'N/A'})`;
    }
    setCurrentContext(newContext);
    // Add a context update message to the chat if it changed significantly
    if (newContext !== currentContext && messages.length > 0) { // Avoid adding on initial load
         setMessages(prevMessages => [...prevMessages, {
             role: 'system', // Use a system role for context messages
             content: [{ type: 'text', text: `Context updated: ${newContext}` }],
             isContextUpdate: true // Custom flag
         }]);
    }

  }, [initialContext, selectedMolecule, currentContext, messages.length]); // Added dependencies

  // Add initial assistant message based on the initial context
  useEffect(() => {
    setMessages([
      {
        role: 'assistant',
        content: [{ type: 'text', text: `Hello! I'm here to help with drug design. Context: ${currentContext}` }]
      }
    ]);
  // Depend only on initialContext to set the *very first* message
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
    let updatedMessages = [...messages, newUserMessage]; // Start with user message

    // Prepend current context implicitly for the API call *if* a molecule is selected
    // This helps Claude remember what is being discussed without cluttering the visible chat too much
    let apiMessages = [...messages, newUserMessage].map(msg => ({ role: msg.role, content: msg.content }));
    if (selectedMolecule) {
        const contextMessageForApi = {
            role: 'user', // Pretend user mentioned it for context
            content: [{ type: 'text', text: `[Current Focus: ${selectedMolecule.name}, SMILES: ${selectedMolecule.smiles}]` }]
        };
        // Add context *before* the latest user message for better flow
        apiMessages.splice(apiMessages.length - 1, 0, contextMessageForApi);
    }
    
    setMessages(updatedMessages); // Update UI immediately with user message
    setInput('');
    setIsLoading(true);
    setError(null);

    try {
      // Use the apiMessages which might include the implicit context
      const response = await claudeAPI.continueChat(apiMessages, currentContext);
      
      // --- Corrected Response Check --- 
      // Check for the actual content array from Anthropic's response
      if (response.data && response.data.content && Array.isArray(response.data.content) && response.data.content.length > 0) {
        // --- Corrected Message Construction --- 
        // Use response.data.content directly as it's the message content array
        const assistantMessage = { 
            role: 'assistant', 
            content: response.data.content 
        }; 
        setMessages(prevMessages => [...prevMessages, assistantMessage]);
        
        // --- Corrected SMILES Extraction Target --- 
        // Check response for SMILES within the first text block
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
        // Log the unexpected structure for debugging
        console.error('Invalid or unexpected response structure from /chat API:', response.data);
        throw new Error('Invalid response format from API');
      }
    } catch (err) {
      console.error('Error sending chat message:', err);
      const errorText = err.message || 'Failed to get response from AI';
      setError(errorText);
      setMessages(prevMessages => [...prevMessages, { role: 'assistant', content: [{ type: 'text', text: `Error: ${errorText}`}], isError: true }]);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <Paper className={classes.root} elevation={1}>
      <List className={classes.messageList}>
        {messages.map((msg, index) => (
          // Don't render system context update messages
          !msg.isContextUpdate && (
            <ListItem key={index} className={msg.role === 'user' ? classes.userMessage : classes.assistantMessage}>
              {msg.role === 'assistant' && (
                <Avatar className={classes.avatar}>
                  <AssistantIcon fontSize="small"/>
                </Avatar>
              )}
              <Box className={`${classes.messageBubble} ${msg.role === 'user' ? classes.userBubble : msg.isError ? classes.errorBubble : classes.assistantBubble}`}> {/* Added errorBubble style */}
                {Array.isArray(msg.content) ? 
                  msg.content.map((part, partIndex) => (
                    part.type === 'text' ? <Typography key={partIndex} variant="body1">{part.text}</Typography> : null
                  ))
                  : <Typography variant="body1">{typeof msg.content === 'string' ? msg.content : 'Invalid message format'}</Typography>
                }
              </Box>
              {msg.role === 'user' && (
                <Avatar className={`${classes.avatar} ${classes.userAvatar}`}>
                  <PersonIcon fontSize="small"/>
                </Avatar>
              )}
            </ListItem>
          )
        ))}
        {isLoading && (
          <ListItem className={classes.assistantMessage}>
            <Avatar className={classes.avatar}>
              <AssistantIcon fontSize="small"/>
            </Avatar>
            <Box className={`${classes.messageBubble} ${classes.assistantBubble}`}>
              <CircularProgress size={20} />
            </Box>
          </ListItem>
        )}
        {error && (
          <ListItem>
            <Alert severity="error" style={{ width: '100%' }}>{error}</Alert>
          </ListItem>
        )}
        <div ref={messagesEndRef} />
      </List>
      <form onSubmit={handleSendMessage} className={classes.messageInputContainer}>
        <TextField
          className={classes.inputField}
          variant="outlined"
          size="small"
          placeholder="Ask Claude about the molecule..."
          value={input}
          onChange={handleInputChange}
          disabled={isLoading}
          multiline
          maxRows={3}
        />
        <IconButton type="submit" color="primary" disabled={isLoading || !input.trim()}>
          <SendIcon />
        </IconButton>
      </form>
    </Paper>
  );
};

export default AIChatInterface; 