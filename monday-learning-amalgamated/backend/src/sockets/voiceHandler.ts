export const handleVoiceQuery = (socket: any, io: any) => {
  socket.on('voice_command', async (data: any) => {
    console.log('Voice command received:', data)
    
    // Process the voice command
    const response = {
      type: 'monday_response',
      content: `Received command: ${data.command}`,
      timestamp: Date.now()
    }
    
    socket.emit('monday_response', response)
  })
} 