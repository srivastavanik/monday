const router = {
  post: (path: string, handler: any) => {},
  get: (path: string, handler: any) => {}
}

router.post('/start', async (req: any, res: any) => {
  res.json({ success: true, message: 'Session started' })
})

router.post('/end', async (req: any, res: any) => {
  res.json({ success: true, message: 'Session ended' })
})

router.get('/status', async (req: any, res: any) => {
  res.json({ success: true, message: 'Session status' })
})

export default router 