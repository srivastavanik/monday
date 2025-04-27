const express = require('express');
const router = express.Router();
const jwt = require('jsonwebtoken');
const bcrypt = require('bcrypt');
const { v4: uuidv4 } = require('uuid');
const path = require('path');
const fs = require('fs').promises;

// Constants
const JWT_SECRET = process.env.JWT_SECRET || 'your_jwt_secret_should_be_in_env';
const TOKEN_EXPIRY = '24h';
const USERS_FILE = path.join(__dirname, '../data/users.json');

// Ensure the data directory exists and create users.json if it doesn't exist
(async () => {
  try {
    const dataDir = path.join(__dirname, '../data');
    await fs.mkdir(dataDir, { recursive: true });
    
    try {
      await fs.access(USERS_FILE);
    } catch (err) {
      // File doesn't exist, create it with default admin user
      const defaultUsers = [
        {
          id: uuidv4(),
          username: 'admin',
          email: 'admin@example.com',
          passwordHash: await bcrypt.hash('admin123', 10),
          role: 'admin',
          createdAt: new Date().toISOString()
        }
      ];
      
      await fs.writeFile(USERS_FILE, JSON.stringify(defaultUsers, null, 2));
      console.log('Created users.json with default admin user');
    }
  } catch (err) {
    console.error('Error initializing auth data:', err);
  }
})();

// Helper function to read users
const getUsers = async () => {
  try {
    const data = await fs.readFile(USERS_FILE, 'utf8');
    return JSON.parse(data);
  } catch (err) {
    console.error('Error reading users:', err);
    return [];
  }
};

// Helper function to write users
const saveUsers = async (users) => {
  try {
    await fs.writeFile(USERS_FILE, JSON.stringify(users, null, 2));
    return true;
  } catch (err) {
    console.error('Error saving users:', err);
    return false;
  }
};

// Middleware to verify JWT token
const authenticateToken = async (req, res, next) => {
  const authHeader = req.headers['authorization'];
  const token = authHeader && authHeader.split(' ')[1];
  
  if (!token) {
    return res.status(401).json({ error: 'Authentication token is required' });
  }
  
  try {
    const user = jwt.verify(token, JWT_SECRET);
    
    // Check if user still exists
    const users = await getUsers();
    const userExists = users.find(u => u.id === user.id);
    
    if (!userExists) {
      return res.status(401).json({ error: 'User no longer exists' });
    }
    
    req.user = user;
    next();
  } catch (err) {
    return res.status(403).json({ error: 'Invalid or expired token' });
  }
};

// Register a new user
router.post('/register', async (req, res) => {
  try {
    const { username, email, password } = req.body;
    
    if (!username || !email || !password) {
      return res.status(400).json({ error: 'Username, email, and password are required' });
    }
    
    // Input validation
    if (username.length < 3) {
      return res.status(400).json({ error: 'Username must be at least 3 characters long' });
    }
    
    if (password.length < 6) {
      return res.status(400).json({ error: 'Password must be at least 6 characters long' });
    }
    
    if (!/^\S+@\S+\.\S+$/.test(email)) {
      return res.status(400).json({ error: 'Invalid email format' });
    }
    
    // Check if user already exists
    const users = await getUsers();
    
    if (users.some(user => user.username === username)) {
      return res.status(400).json({ error: 'Username already exists' });
    }
    
    if (users.some(user => user.email === email)) {
      return res.status(400).json({ error: 'Email already exists' });
    }
    
    // Create new user
    const newUser = {
      id: uuidv4(),
      username,
      email,
      passwordHash: await bcrypt.hash(password, 10),
      role: 'user',
      createdAt: new Date().toISOString()
    };
    
    // Save user
    users.push(newUser);
    await saveUsers(users);
    
    // Create JWT token
    const token = jwt.sign(
      { id: newUser.id, username: newUser.username, email: newUser.email, role: newUser.role },
      JWT_SECRET,
      { expiresIn: TOKEN_EXPIRY }
    );
    
    // Return user info (without password) and token
    return res.status(201).json({
      user: {
        id: newUser.id,
        username: newUser.username,
        email: newUser.email,
        role: newUser.role
      },
      token
    });
    
  } catch (error) {
    console.error('Registration error:', error);
    return res.status(500).json({ error: 'Error registering user' });
  }
});

// Login
router.post('/login', async (req, res) => {
  try {
    const { username, password } = req.body;
    
    if (!username || !password) {
      return res.status(400).json({ error: 'Username and password are required' });
    }
    
    // Find user
    const users = await getUsers();
    const user = users.find(u => u.username === username);
    
    if (!user) {
      return res.status(401).json({ error: 'Invalid username or password' });
    }
    
    // Check password
    const passwordValid = await bcrypt.compare(password, user.passwordHash);
    
    if (!passwordValid) {
      return res.status(401).json({ error: 'Invalid username or password' });
    }
    
    // Create JWT token
    const token = jwt.sign(
      { id: user.id, username: user.username, email: user.email, role: user.role },
      JWT_SECRET,
      { expiresIn: TOKEN_EXPIRY }
    );
    
    // Return user info and token
    return res.json({
      user: {
        id: user.id,
        username: user.username,
        email: user.email,
        role: user.role
      },
      token
    });
    
  } catch (error) {
    console.error('Login error:', error);
    return res.status(500).json({ error: 'Error logging in' });
  }
});

// Get current user
router.get('/me', authenticateToken, async (req, res) => {
  try {
    const users = await getUsers();
    const user = users.find(u => u.id === req.user.id);
    
    if (!user) {
      return res.status(404).json({ error: 'User not found' });
    }
    
    return res.json({
      id: user.id,
      username: user.username,
      email: user.email,
      role: user.role
    });
    
  } catch (error) {
    console.error('Get user error:', error);
    return res.status(500).json({ error: 'Error getting user information' });
  }
});

// Change password
router.post('/change-password', authenticateToken, async (req, res) => {
  try {
    const { currentPassword, newPassword } = req.body;
    
    if (!currentPassword || !newPassword) {
      return res.status(400).json({ error: 'Current password and new password are required' });
    }
    
    if (newPassword.length < 6) {
      return res.status(400).json({ error: 'New password must be at least 6 characters long' });
    }
    
    // Get users
    const users = await getUsers();
    const userIndex = users.findIndex(u => u.id === req.user.id);
    
    if (userIndex === -1) {
      return res.status(404).json({ error: 'User not found' });
    }
    
    // Verify current password
    const passwordValid = await bcrypt.compare(currentPassword, users[userIndex].passwordHash);
    
    if (!passwordValid) {
      return res.status(401).json({ error: 'Current password is incorrect' });
    }
    
    // Update password
    users[userIndex].passwordHash = await bcrypt.hash(newPassword, 10);
    users[userIndex].updatedAt = new Date().toISOString();
    
    // Save users
    await saveUsers(users);
    
    return res.json({ message: 'Password changed successfully' });
    
  } catch (error) {
    console.error('Change password error:', error);
    return res.status(500).json({ error: 'Error changing password' });
  }
});

// Update profile
router.put('/profile', authenticateToken, async (req, res) => {
  try {
    const { email } = req.body;
    
    if (!email) {
      return res.status(400).json({ error: 'Email is required' });
    }
    
    if (!/^\S+@\S+\.\S+$/.test(email)) {
      return res.status(400).json({ error: 'Invalid email format' });
    }
    
    // Get users
    const users = await getUsers();
    const userIndex = users.findIndex(u => u.id === req.user.id);
    
    if (userIndex === -1) {
      return res.status(404).json({ error: 'User not found' });
    }
    
    // Check if email is already used by another user
    if (users.some(u => u.email === email && u.id !== req.user.id)) {
      return res.status(400).json({ error: 'Email already in use' });
    }
    
    // Update profile
    users[userIndex].email = email;
    users[userIndex].updatedAt = new Date().toISOString();
    
    // Save users
    await saveUsers(users);
    
    return res.json({
      id: users[userIndex].id,
      username: users[userIndex].username,
      email: users[userIndex].email,
      role: users[userIndex].role
    });
    
  } catch (error) {
    console.error('Update profile error:', error);
    return res.status(500).json({ error: 'Error updating profile' });
  }
});

// Admin routes

// Get all users (admin only)
router.get('/users', authenticateToken, async (req, res) => {
  try {
    // Check if user is admin
    if (req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Forbidden: Admin access required' });
    }
    
    const users = await getUsers();
    
    // Return users without password hashes
    const sanitizedUsers = users.map(({ passwordHash, ...user }) => user);
    
    return res.json(sanitizedUsers);
    
  } catch (error) {
    console.error('Get all users error:', error);
    return res.status(500).json({ error: 'Error getting users' });
  }
});

// Create user (admin only)
router.post('/users', authenticateToken, async (req, res) => {
  try {
    // Check if user is admin
    if (req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Forbidden: Admin access required' });
    }
    
    const { username, email, password, role = 'user' } = req.body;
    
    if (!username || !email || !password) {
      return res.status(400).json({ error: 'Username, email, and password are required' });
    }
    
    // Input validation
    if (username.length < 3) {
      return res.status(400).json({ error: 'Username must be at least 3 characters long' });
    }
    
    if (password.length < 6) {
      return res.status(400).json({ error: 'Password must be at least 6 characters long' });
    }
    
    if (!/^\S+@\S+\.\S+$/.test(email)) {
      return res.status(400).json({ error: 'Invalid email format' });
    }
    
    if (!['user', 'admin'].includes(role)) {
      return res.status(400).json({ error: 'Invalid role' });
    }
    
    // Check if user already exists
    const users = await getUsers();
    
    if (users.some(user => user.username === username)) {
      return res.status(400).json({ error: 'Username already exists' });
    }
    
    if (users.some(user => user.email === email)) {
      return res.status(400).json({ error: 'Email already exists' });
    }
    
    // Create new user
    const newUser = {
      id: uuidv4(),
      username,
      email,
      passwordHash: await bcrypt.hash(password, 10),
      role,
      createdAt: new Date().toISOString()
    };
    
    // Save user
    users.push(newUser);
    await saveUsers(users);
    
    // Return user info (without password)
    return res.status(201).json({
      id: newUser.id,
      username: newUser.username,
      email: newUser.email,
      role: newUser.role
    });
    
  } catch (error) {
    console.error('Create user error:', error);
    return res.status(500).json({ error: 'Error creating user' });
  }
});

// Delete user (admin only)
router.delete('/users/:id', authenticateToken, async (req, res) => {
  try {
    // Check if user is admin
    if (req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Forbidden: Admin access required' });
    }
    
    const { id } = req.params;
    
    // Get users
    const users = await getUsers();
    const userIndex = users.findIndex(u => u.id === id);
    
    if (userIndex === -1) {
      return res.status(404).json({ error: 'User not found' });
    }
    
    // Prevent deleting self
    if (id === req.user.id) {
      return res.status(400).json({ error: 'Cannot delete your own account' });
    }
    
    // Remove user
    users.splice(userIndex, 1);
    
    // Save users
    await saveUsers(users);
    
    return res.json({ message: 'User deleted successfully' });
    
  } catch (error) {
    console.error('Delete user error:', error);
    return res.status(500).json({ error: 'Error deleting user' });
  }
});

// Export the authenticateToken middleware for use in other routes
module.exports = {
  router,
  authenticateToken
}; 