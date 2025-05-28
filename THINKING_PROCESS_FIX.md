# Thinking Process Display Fix

## Issue
For reasoning and deep research modes, the thinking process was being displayed in the main panel (middle) instead of the left visualization panel. Additionally, the response was being cut off.

## Root Cause
1. The backend was properly extracting and filtering thinking tags, but wasn't passing the thinking data to the frontend
2. The frontend wasn't configured to display the complete thinking process in the visualization panel

## Fix Applied

### 1. Backend Changes (`backend-server.js`)
- Added `thinking: perplexityResponse.thinking || ''` to the voice_response data
- Ensures the thinking process is sent separately from the filtered main response

### 2. Frontend Type Updates (`monday-perplexity-system.tsx`)
- Added `thinking?: string` to the WebSocketMessage interface
- Updated handleVoiceResponse to include thinking in visualizationData

### 3. Visualization Panel Updates (`adaptive-visualization-panel.tsx`)
- Added `thinking?: string` to AdaptiveVisualizationPanelProps
- Updated ProgressiveThinkingDisplay to accept and display the thinking prop
- Added a new section to display the complete thinking process

### 4. Component Integration
- Updated the AdaptiveVisualizationPanel call to pass `thinking={visualizationData?.thinking}`

## Result

### Reasoning Mode
- **Left Panel**: Shows the complete reasoning process in a formatted display
- **Middle Panel**: Shows only the final response without thinking tags
- **No cutoff**: Full response is displayed

### Deep Research Mode
- **Left Panel**: Shows the complete research process
- **Middle Panel**: Shows only the final response without research tags
- **No cutoff**: Full response is displayed

## How It Works

1. Perplexity API returns response with `<think>` or `<research>` tags
2. Backend extracts thinking process using `extractThinkingProcess()`
3. Backend filters out thinking tags using `filterThinkingTags()`
4. Backend sends:
   - `message`: Filtered response (for main panel)
   - `data.thinking`: Complete thinking process (for left panel)
5. Frontend displays thinking in left panel and response in middle panel

## Testing

You can test by saying:
- "Hey Monday, think about why the sky is blue" - See reasoning in left panel
- "Hey Monday, research into the history of coffee" - See research process in left panel

The thinking process will appear in the left panel with proper formatting, while the main response appears in the middle panel without any thinking tags or cutoffs. 