const express = require('express');
const router = express.Router();
const path = require('path');
const fs = require('fs').promises;
const { v4: uuidv4 } = require('uuid');
const { authenticateToken } = require('./auth');

// Constants
const DATA_DIR = path.join(__dirname, '../data');
const PROJECTS_FILE = path.join(DATA_DIR, 'projects.json');
const COMMENTS_FILE = path.join(DATA_DIR, 'comments.json');
const ACTIVITIES_FILE = path.join(DATA_DIR, 'activities.json');

// Ensure data files exist
(async () => {
  try {
    // Create data directory if it doesn't exist
    await fs.mkdir(DATA_DIR, { recursive: true });
    
    // Initialize projects.json if it doesn't exist
    try {
      await fs.access(PROJECTS_FILE);
    } catch (err) {
      await fs.writeFile(PROJECTS_FILE, JSON.stringify([], null, 2));
      console.log('Created projects.json');
    }
    
    // Initialize comments.json if it doesn't exist
    try {
      await fs.access(COMMENTS_FILE);
    } catch (err) {
      await fs.writeFile(COMMENTS_FILE, JSON.stringify([], null, 2));
      console.log('Created comments.json');
    }
    
    // Initialize activities.json if it doesn't exist
    try {
      await fs.access(ACTIVITIES_FILE);
    } catch (err) {
      await fs.writeFile(ACTIVITIES_FILE, JSON.stringify([], null, 2));
      console.log('Created activities.json');
    }
  } catch (err) {
    console.error('Error initializing collaboration data files:', err);
  }
})();

// Helper functions
const getProjects = async () => {
  try {
    const data = await fs.readFile(PROJECTS_FILE, 'utf8');
    return JSON.parse(data);
  } catch (err) {
    console.error('Error reading projects:', err);
    return [];
  }
};

const saveProjects = async (projects) => {
  try {
    await fs.writeFile(PROJECTS_FILE, JSON.stringify(projects, null, 2));
    return true;
  } catch (err) {
    console.error('Error saving projects:', err);
    return false;
  }
};

const getComments = async () => {
  try {
    const data = await fs.readFile(COMMENTS_FILE, 'utf8');
    return JSON.parse(data);
  } catch (err) {
    console.error('Error reading comments:', err);
    return [];
  }
};

const saveComments = async (comments) => {
  try {
    await fs.writeFile(COMMENTS_FILE, JSON.stringify(comments, null, 2));
    return true;
  } catch (err) {
    console.error('Error saving comments:', err);
    return false;
  }
};

const getActivities = async () => {
  try {
    const data = await fs.readFile(ACTIVITIES_FILE, 'utf8');
    return JSON.parse(data);
  } catch (err) {
    console.error('Error reading activities:', err);
    return [];
  }
};

const saveActivities = async (activities) => {
  try {
    await fs.writeFile(ACTIVITIES_FILE, JSON.stringify(activities, null, 2));
    return true;
  } catch (err) {
    console.error('Error saving activities:', err);
    return false;
  }
};

const logActivity = async (userId, username, projectId, action, details = {}) => {
  try {
    const activities = await getActivities();
    
    const activity = {
      id: uuidv4(),
      userId,
      username,
      projectId,
      action,
      details,
      timestamp: new Date().toISOString()
    };
    
    activities.push(activity);
    
    // Keep only the last 1000 activities to prevent the file from growing too large
    if (activities.length > 1000) {
      activities.sort((a, b) => new Date(b.timestamp) - new Date(a.timestamp));
      activities.splice(1000);
    }
    
    await saveActivities(activities);
    
    return activity;
  } catch (err) {
    console.error('Error logging activity:', err);
    return null;
  }
};

// Check if user is a member of a project
const isProjectMember = (project, userId) => {
  return project.members.some(member => member.userId === userId);
};

// Check if user is an admin of a project
const isProjectAdmin = (project, userId) => {
  const member = project.members.find(member => member.userId === userId);
  return member && member.role === 'admin';
};

// Get all projects (that the user is a member of)
router.get('/projects', authenticateToken, async (req, res) => {
  try {
    const userId = req.user.id;
    const projects = await getProjects();
    
    // Filter projects to only include those the user is a member of
    const userProjects = projects.filter(project => 
      isProjectMember(project, userId) || req.user.role === 'admin'
    );
    
    return res.json(userProjects);
  } catch (error) {
    console.error('Error getting projects:', error);
    return res.status(500).json({ error: 'Error getting projects' });
  }
});

// Get a specific project
router.get('/projects/:id', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    const userId = req.user.id;
    const projects = await getProjects();
    
    const project = projects.find(p => p.id === id);
    
    if (!project) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project or an admin
    if (!isProjectMember(project, userId) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    return res.json(project);
  } catch (error) {
    console.error('Error getting project:', error);
    return res.status(500).json({ error: 'Error getting project' });
  }
});

// Create a new project
router.post('/projects', authenticateToken, async (req, res) => {
  try {
    const { name, description } = req.body;
    
    if (!name) {
      return res.status(400).json({ error: 'Project name is required' });
    }
    
    const projects = await getProjects();
    
    // Create new project
    const newProject = {
      id: uuidv4(),
      name,
      description: description || '',
      createdBy: req.user.id,
      createdByUsername: req.user.username,
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      members: [
        {
          userId: req.user.id,
          username: req.user.username,
          role: 'admin',
          joinedAt: new Date().toISOString()
        }
      ],
      molecules: []
    };
    
    projects.push(newProject);
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      newProject.id,
      'PROJECT_CREATED',
      { projectName: name }
    );
    
    return res.status(201).json(newProject);
  } catch (error) {
    console.error('Error creating project:', error);
    return res.status(500).json({ error: 'Error creating project' });
  }
});

// Update a project
router.put('/projects/:id', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    const { name, description } = req.body;
    
    if (!name) {
      return res.status(400).json({ error: 'Project name is required' });
    }
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === id);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a project admin or site admin
    if (!isProjectAdmin(projects[projectIndex], req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Only project admins can update project details' });
    }
    
    // Update project
    projects[projectIndex].name = name;
    projects[projectIndex].description = description || '';
    projects[projectIndex].updatedAt = new Date().toISOString();
    
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      id,
      'PROJECT_UPDATED',
      { projectName: name }
    );
    
    return res.json(projects[projectIndex]);
  } catch (error) {
    console.error('Error updating project:', error);
    return res.status(500).json({ error: 'Error updating project' });
  }
});

// Delete a project
router.delete('/projects/:id', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === id);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a project admin or site admin
    if (!isProjectAdmin(projects[projectIndex], req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Only project admins can delete projects' });
    }
    
    const projectName = projects[projectIndex].name;
    
    // Remove project
    projects.splice(projectIndex, 1);
    
    await saveProjects(projects);
    
    // Delete project comments
    const comments = await getComments();
    const filteredComments = comments.filter(comment => comment.projectId !== id);
    await saveComments(filteredComments);
    
    // Log activity (even though the project is deleted, we keep the activity record)
    await logActivity(
      req.user.id,
      req.user.username,
      id,
      'PROJECT_DELETED',
      { projectName }
    );
    
    return res.json({ message: 'Project deleted successfully' });
  } catch (error) {
    console.error('Error deleting project:', error);
    return res.status(500).json({ error: 'Error deleting project' });
  }
});

// Add a molecule to a project
router.post('/projects/:id/molecules', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    const { molecule } = req.body;
    
    if (!molecule || !molecule.smiles) {
      return res.status(400).json({ error: 'Valid molecule data with SMILES is required' });
    }
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === id);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project
    if (!isProjectMember(projects[projectIndex], req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    // Add molecule with metadata
    const moleculeWithMeta = {
      ...molecule,
      id: molecule.id || uuidv4(),
      addedBy: req.user.id,
      addedByUsername: req.user.username,
      addedAt: new Date().toISOString()
    };
    
    projects[projectIndex].molecules.push(moleculeWithMeta);
    projects[projectIndex].updatedAt = new Date().toISOString();
    
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      id,
      'MOLECULE_ADDED',
      { 
        projectName: projects[projectIndex].name,
        moleculeName: molecule.name || 'Unnamed molecule',
        smiles: molecule.smiles
      }
    );
    
    return res.status(201).json(moleculeWithMeta);
  } catch (error) {
    console.error('Error adding molecule to project:', error);
    return res.status(500).json({ error: 'Error adding molecule to project' });
  }
});

// Remove a molecule from a project
router.delete('/projects/:projectId/molecules/:moleculeId', authenticateToken, async (req, res) => {
  try {
    const { projectId, moleculeId } = req.params;
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === projectId);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project
    if (!isProjectMember(projects[projectIndex], req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    const moleculeIndex = projects[projectIndex].molecules.findIndex(m => m.id === moleculeId);
    
    if (moleculeIndex === -1) {
      return res.status(404).json({ error: 'Molecule not found in project' });
    }
    
    const moleculeName = projects[projectIndex].molecules[moleculeIndex].name || 'Unnamed molecule';
    
    // Remove molecule
    projects[projectIndex].molecules.splice(moleculeIndex, 1);
    projects[projectIndex].updatedAt = new Date().toISOString();
    
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      projectId,
      'MOLECULE_REMOVED',
      { 
        projectName: projects[projectIndex].name,
        moleculeName
      }
    );
    
    return res.json({ message: 'Molecule removed from project successfully' });
  } catch (error) {
    console.error('Error removing molecule from project:', error);
    return res.status(500).json({ error: 'Error removing molecule from project' });
  }
});

// Add a member to a project
router.post('/projects/:id/members', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    const { userId, username, role = 'member' } = req.body;
    
    if (!userId || !username) {
      return res.status(400).json({ error: 'User ID and username are required' });
    }
    
    if (!['admin', 'member'].includes(role)) {
      return res.status(400).json({ error: 'Role must be either "admin" or "member"' });
    }
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === id);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a project admin or site admin
    if (!isProjectAdmin(projects[projectIndex], req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Only project admins can add members' });
    }
    
    // Check if user is already a member
    if (projects[projectIndex].members.some(member => member.userId === userId)) {
      return res.status(400).json({ error: 'User is already a member of this project' });
    }
    
    // Add member
    const newMember = {
      userId,
      username,
      role,
      joinedAt: new Date().toISOString()
    };
    
    projects[projectIndex].members.push(newMember);
    projects[projectIndex].updatedAt = new Date().toISOString();
    
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      id,
      'MEMBER_ADDED',
      { 
        projectName: projects[projectIndex].name,
        memberUsername: username,
        memberRole: role
      }
    );
    
    return res.status(201).json(newMember);
  } catch (error) {
    console.error('Error adding member to project:', error);
    return res.status(500).json({ error: 'Error adding member to project' });
  }
});

// Remove a member from a project
router.delete('/projects/:projectId/members/:userId', authenticateToken, async (req, res) => {
  try {
    const { projectId, userId } = req.params;
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === projectId);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a project admin or site admin or removing themselves
    const isAdmin = isProjectAdmin(projects[projectIndex], req.user.id) || req.user.role === 'admin';
    const isSelf = req.user.id === userId;
    
    if (!isAdmin && !isSelf) {
      return res.status(403).json({ error: 'Only project admins can remove members' });
    }
    
    const memberIndex = projects[projectIndex].members.findIndex(m => m.userId === userId);
    
    if (memberIndex === -1) {
      return res.status(404).json({ error: 'Member not found in project' });
    }
    
    // Prevent removing the last admin
    const isRemovingAdmin = projects[projectIndex].members[memberIndex].role === 'admin';
    const adminCount = projects[projectIndex].members.filter(m => m.role === 'admin').length;
    
    if (isRemovingAdmin && adminCount === 1) {
      return res.status(400).json({ error: 'Cannot remove the last admin from the project' });
    }
    
    const memberUsername = projects[projectIndex].members[memberIndex].username;
    
    // Remove member
    projects[projectIndex].members.splice(memberIndex, 1);
    projects[projectIndex].updatedAt = new Date().toISOString();
    
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      projectId,
      'MEMBER_REMOVED',
      { 
        projectName: projects[projectIndex].name,
        memberUsername
      }
    );
    
    return res.json({ message: 'Member removed from project successfully' });
  } catch (error) {
    console.error('Error removing member from project:', error);
    return res.status(500).json({ error: 'Error removing member from project' });
  }
});

// Update a member's role
router.put('/projects/:projectId/members/:userId', authenticateToken, async (req, res) => {
  try {
    const { projectId, userId } = req.params;
    const { role } = req.body;
    
    if (!role || !['admin', 'member'].includes(role)) {
      return res.status(400).json({ error: 'Role must be either "admin" or "member"' });
    }
    
    const projects = await getProjects();
    
    const projectIndex = projects.findIndex(p => p.id === projectId);
    
    if (projectIndex === -1) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a project admin or site admin
    if (!isProjectAdmin(projects[projectIndex], req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'Only project admins can update member roles' });
    }
    
    const memberIndex = projects[projectIndex].members.findIndex(m => m.userId === userId);
    
    if (memberIndex === -1) {
      return res.status(404).json({ error: 'Member not found in project' });
    }
    
    // Prevent downgrading the last admin
    const isDowngradingAdmin = projects[projectIndex].members[memberIndex].role === 'admin' && role === 'member';
    const adminCount = projects[projectIndex].members.filter(m => m.role === 'admin').length;
    
    if (isDowngradingAdmin && adminCount === 1) {
      return res.status(400).json({ error: 'Cannot downgrade the last admin of the project' });
    }
    
    const memberUsername = projects[projectIndex].members[memberIndex].username;
    const oldRole = projects[projectIndex].members[memberIndex].role;
    
    // Update role
    projects[projectIndex].members[memberIndex].role = role;
    projects[projectIndex].updatedAt = new Date().toISOString();
    
    await saveProjects(projects);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      projectId,
      'MEMBER_ROLE_UPDATED',
      { 
        projectName: projects[projectIndex].name,
        memberUsername,
        oldRole,
        newRole: role
      }
    );
    
    return res.json(projects[projectIndex].members[memberIndex]);
  } catch (error) {
    console.error('Error updating member role:', error);
    return res.status(500).json({ error: 'Error updating member role' });
  }
});

// Get project comments
router.get('/projects/:id/comments', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    
    const projects = await getProjects();
    
    const project = projects.find(p => p.id === id);
    
    if (!project) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project or an admin
    if (!isProjectMember(project, req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    const comments = await getComments();
    
    // Get comments for this project
    const projectComments = comments
      .filter(comment => comment.projectId === id)
      .sort((a, b) => new Date(a.timestamp) - new Date(b.timestamp));
    
    return res.json(projectComments);
  } catch (error) {
    console.error('Error getting project comments:', error);
    return res.status(500).json({ error: 'Error getting project comments' });
  }
});

// Add a comment to a project
router.post('/projects/:id/comments', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    const { text, moleculeId } = req.body;
    
    if (!text) {
      return res.status(400).json({ error: 'Comment text is required' });
    }
    
    const projects = await getProjects();
    
    const project = projects.find(p => p.id === id);
    
    if (!project) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project or an admin
    if (!isProjectMember(project, req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    // Check if moleculeId exists in the project if provided
    if (moleculeId && !project.molecules.some(m => m.id === moleculeId)) {
      return res.status(404).json({ error: 'Molecule not found in project' });
    }
    
    const comments = await getComments();
    
    // Create new comment
    const newComment = {
      id: uuidv4(),
      projectId: id,
      moleculeId: moleculeId || null,
      text,
      userId: req.user.id,
      username: req.user.username,
      timestamp: new Date().toISOString()
    };
    
    comments.push(newComment);
    
    await saveComments(comments);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      id,
      'COMMENT_ADDED',
      { 
        projectName: project.name,
        moleculeId,
        commentId: newComment.id
      }
    );
    
    return res.status(201).json(newComment);
  } catch (error) {
    console.error('Error adding comment:', error);
    return res.status(500).json({ error: 'Error adding comment' });
  }
});

// Delete a comment
router.delete('/projects/:projectId/comments/:commentId', authenticateToken, async (req, res) => {
  try {
    const { projectId, commentId } = req.params;
    
    const projects = await getProjects();
    
    const project = projects.find(p => p.id === projectId);
    
    if (!project) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project or an admin
    if (!isProjectMember(project, req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    const comments = await getComments();
    
    const commentIndex = comments.findIndex(c => c.id === commentId && c.projectId === projectId);
    
    if (commentIndex === -1) {
      return res.status(404).json({ error: 'Comment not found' });
    }
    
    // Check if user is the comment author, a project admin, or a site admin
    const isAuthor = comments[commentIndex].userId === req.user.id;
    const isProjectAdmin = project.members.some(m => m.userId === req.user.id && m.role === 'admin');
    
    if (!isAuthor && !isProjectAdmin && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You can only delete your own comments' });
    }
    
    // Remove comment
    comments.splice(commentIndex, 1);
    
    await saveComments(comments);
    
    // Log activity
    await logActivity(
      req.user.id,
      req.user.username,
      projectId,
      'COMMENT_DELETED',
      { 
        projectName: project.name,
        commentId
      }
    );
    
    return res.json({ message: 'Comment deleted successfully' });
  } catch (error) {
    console.error('Error deleting comment:', error);
    return res.status(500).json({ error: 'Error deleting comment' });
  }
});

// Get project activity
router.get('/projects/:id/activity', authenticateToken, async (req, res) => {
  try {
    const { id } = req.params;
    
    const projects = await getProjects();
    
    const project = projects.find(p => p.id === id);
    
    if (!project) {
      return res.status(404).json({ error: 'Project not found' });
    }
    
    // Check if user is a member of the project or an admin
    if (!isProjectMember(project, req.user.id) && req.user.role !== 'admin') {
      return res.status(403).json({ error: 'You do not have access to this project' });
    }
    
    const activities = await getActivities();
    
    // Get activities for this project, newest first
    const projectActivities = activities
      .filter(activity => activity.projectId === id)
      .sort((a, b) => new Date(b.timestamp) - new Date(a.timestamp));
    
    return res.json(projectActivities);
  } catch (error) {
    console.error('Error getting project activity:', error);
    return res.status(500).json({ error: 'Error getting project activity' });
  }
});

module.exports = router; 