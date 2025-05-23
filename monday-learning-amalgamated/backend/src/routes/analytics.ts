const router = {
  get: (path: string, handler: any) => {},
  post: (path: string, handler: any) => {}
}

router.get('/metrics', async (req: any, res: any) => {
  res.json({ success: true, message: 'Analytics metrics' })
})

router.post('/track', async (req: any, res: any) => {
  res.json({ success: true, message: 'Event tracked' })
})

export default router 