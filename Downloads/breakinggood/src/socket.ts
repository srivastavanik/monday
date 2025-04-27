import { serve } from 'bun';

console.log('WebSocket server started...');

serve<{
  channel: string;
}> ({
  fetch(req, server) {
    const url = new URL(req.url);
    const channel = url.pathname.slice(1);

    if (
      server.upgrade(req, {
        data: {
          channel,
        },
      })
    )
      return; // do not return a Response

    return new Response('Upgrade failed :(', { status: 500 });
  },
  websocket: {
    open(ws) {
      console.log(`Client connected to channel: ${ws.data.channel}`);
      ws.subscribe(ws.data.channel);
    },
    message(ws, message) {
      console.log(`Received message on channel ${ws.data.channel}:`, message);
      // the server broadcasts messages to all clients subscribed to the same channel
      ws.publish(ws.data.channel, message);
    },
    close(ws) {
      console.log(`Client disconnected from channel: ${ws.data.channel}`);
      ws.unsubscribe(ws.data.channel);
    },
  },
  port: 3001,
  // uncomment this to allow connections in windows wsl
  hostname: "0.0.0.0",
}); 