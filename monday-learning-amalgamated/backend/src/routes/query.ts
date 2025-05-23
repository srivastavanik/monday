import { Router } from 'express'

const router = Router()

router.post('/basic', async (req, res) => {
  try {
    res.json({ success: true, message: 'Basic query endpoint' })
  } catch (error) {
    res.status(500).json({ success: false, error: 'Internal server error' })
  }
})

router.post('/reasoning', async (req, res) => {
  try {
    res.json({ success: true, message: 'Reasoning query endpoint' })
  } catch (error) {
    res.status(500).json({ success: false, error: 'Internal server error' })
  }
})

router.post('/research', async (req, res) => {
  try {
    res.json({ success: true, message: 'Research query endpoint' })
  } catch (error) {
    res.status(500).json({ success: false, error: 'Internal server error' })
  }
})

export default router 