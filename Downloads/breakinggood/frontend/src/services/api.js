import axios from "axios";

// Create axios instance with base URL and default headers
const api = axios.create({
  baseURL: process.env.REACT_APP_API_URL || "/api",
  headers: {
    "Content-Type": "application/json",
  },
});

// Add request interceptor for authentication
api.interceptors.request.use(
  (config) => {
    const token = localStorage.getItem("auth_token");
    if (token) {
      config.headers.Authorization = `Bearer ${token}`;
    }
    return config;
  },
  (error) => Promise.reject(error)
);

// Add response interceptor for error handling
api.interceptors.response.use(
  (response) => response,
  (error) => {
    // Handle authentication errors
    if (error.response && error.response.status === 401) {
      localStorage.removeItem("auth_token");
      window.location.href = "/login";
      return Promise.reject(
        new Error("Authentication failed. Please log in again.")
      );
    }
    return Promise.reject(error);
  }
);

// Authentication API
export const authAPI = {
  login: (credentials) => api.post("/auth/login", credentials),
  register: (userData) => api.post("/auth/register", userData),
  logout: () => api.post("/auth/logout"),
  getCurrentUser: () => api.get("/auth/me"),
};

// Drug Design API
export const drugDesignAPI = {
  // Generate molecule based on requirements
  generateMolecule: (requirements) =>
    api.post("/drug-design/generate", requirements),

  // Refine an existing molecule
  refineMolecule: (data) => api.post("/drug-design/refine", data),

  // Analyze regulatory pathway
  regulatoryAnalysis: (data) =>
    api.post("/drug-design/regulatory-analysis", data),

  // Save molecule to user's collection
  saveMolecule: (molecule) => api.post("/drug-design/molecules", molecule),

  // Get user's saved molecules
  getSavedMolecules: () => api.get("/drug-design/molecules"),

  // Get specific molecule by ID
  getMolecule: (id) => api.get(`/drug-design/molecules/${id}`),

  // Update molecule
  updateMolecule: (id, data) => api.put(`/drug-design/molecules/${id}`, data),

  // Delete molecule
  deleteMolecule: (id) => api.delete(`/drug-design/molecules/${id}`),
};

// Simulation API
export const simulationAPI = {
  // Run RDKit operations
  runRDKit: (data) => api.post("/simulation/rdkit", data),

  // Generate 3D structure
  generate3DStructure: (smiles) =>
    api.post("/simulation/3d-structure", { smiles }),

  // Compare two molecules
  compareMolecules: (smiles1, smiles2) =>
    api.post("/simulation/compare", { smiles1, smiles2 }),

  // Calculate binding affinity
  calculateBindingAffinity: (data) =>
    api.post("/simulation/binding-affinity", data),

  // Predict ADMET properties
  predictADMET: (smiles) => api.post("/simulation/admet", { smiles }),

  // Run molecular dynamics simulation
  runMolecularDynamics: (data) =>
    api.post("/simulation/molecular-dynamics", data),

  // Run docking simulation
  runDocking: (data) => api.post("/simulation/docking", data),

  // Get simulation results
  getSimulationResults: (id) => api.get(`/simulation/results/${id}`),

  // Save simulation results
  saveSimulationResults: (data) => api.post("/simulation/results", data),
};

// Literature API
export const literatureAPI = {
  // Search literature using Crossref
  searchLiterature: (query, options) =>
    api.get("/literature/search", { params: { query, ...options } }),

  // Create research note
  createNote: (noteData) => api.post("/literature/notes", noteData),

  // Get research notes for an article
  getNotes: (articleId) => api.get(`/literature/notes/${articleId}`),

  // Update research note
  updateNote: (noteId, noteData) =>
    api.put(`/literature/notes/${noteId}`, noteData),

  // Delete research note
  deleteNote: (noteId) => api.delete(`/literature/notes/${noteId}`),

  // Ask Claude to analyze literature
  analyzeLiterature: (articles, query) =>
    api.post("/literature/analyze", { articles, query }),
};

// Chemical Data API
export const chemicalDataAPI = {
  // Search ChEMBL
  searchChEMBL: (query, options) =>
    api.get("/chemical-data/chembl", { params: { query, ...options } }),

  // Get compound details from ChEMBL
  getChEMBLCompound: (chemblId) => api.get(`/chemical-data/chembl/${chemblId}`),

  // Search by structure similarity
  searchBySimilarity: (smiles, threshold) =>
    api.post("/chemical-data/similarity", { smiles, threshold }),

  // Search by substructure
  searchBySubstructure: (smiles) =>
    api.post("/chemical-data/substructure", { smiles }),

  // Get bioactivity data
  getBioactivityData: (chemblId) =>
    api.get(`/chemical-data/bioactivity/${chemblId}`),

  // Get drug-target interactions
  getDrugTargetInteractions: (chemblId) =>
    api.get(`/chemical-data/drug-target/${chemblId}`),

  // Similarity search with detailed options
  similaritySearch: (query, targets, options = {}) =>
    api.post("/similarity/search", {
      query,
      targets,
      fingerprint: options.fingerprint || "morgan",
      metric: options.metric || "tanimoto",
      threshold: options.threshold || 0.7,
      maxResults: options.maxResults || 50,
    }),

  // Multi-reference similarity search
  multiReferenceSimilaritySearch: (queries, targets, options = {}) =>
    api.post("/similarity/multi-reference", {
      queries,
      targets,
      fingerprint: options.fingerprint || "morgan",
      metric: options.metric || "tanimoto",
      threshold: options.threshold || 0.7,
      maxResults: options.maxResults || 50,
      aggregation: options.aggregation || "max",
    }),

  // Diversity selection from a set of molecules
  diversitySelection: (molecules, options = {}) =>
    api.post("/similarity/diversity", {
      molecules,
      numPicks: options.numPicks || 10,
      fingerprint: options.fingerprint || "morgan",
      metric: options.metric || "tanimoto",
    }),
};

// Regulatory Analysis API
export const regulatoryAPI = {
  // Analyze safety and toxicity
  analyzeSafety: (smiles) => api.post("/regulatory/safety", { smiles }),

  // Analyze pharmaceutical properties
  analyzePharmaceuticalProperties: (smiles) =>
    api.post("/regulatory/pharmaceutical", { smiles }),

  // Get regulatory pathway (DEPRECATED - Use generateReport)
  getRegulatoryPathway: (data) => api.post("/regulatory/pathway", data),

  // Get similar approved drugs (Might need refinement in backend)
  getSimilarApprovedDrugs: (smiles) =>
    api.post("/regulatory/similar-approved", { smiles }),

  // Get development timeline (DEPRECATED - Use generateReport)
  getDevelopmentTimeline: (data) => api.post("/regulatory/timeline", data),

  // Get regulatory requirements (Potentially DEPRECATED if covered by report)
  getRegulatoryRequirements: (indication) =>
    api.get(`/regulatory/requirements/${indication}`),

  // Generate the full regulatory report
  generateRegulatoryReport: (params) => api.post("/regulatory/report", params),
};

// Collaborative features API
export const collaborationAPI = {
  // Get projects
  getProjects: () => api.get("/collaboration/projects"),

  // Get project by ID
  getProject: (id) => api.get(`/collaboration/projects/${id}`),

  // Create project
  createProject: (project) => api.post("/collaboration/projects", project),

  // Update project
  updateProject: (id, data) => api.put(`/collaboration/projects/${id}`, data),

  // Delete project
  deleteProject: (id) => api.delete(`/collaboration/projects/${id}`),

  // Add member to project
  addMember: (projectId, userId, role) =>
    api.post(`/collaboration/projects/${projectId}/members`, { userId, role }),

  // Remove member from project
  removeMember: (projectId, userId) =>
    api.delete(`/collaboration/projects/${projectId}/members/${userId}`),

  // Get project comments
  getComments: (projectId) =>
    api.get(`/collaboration/projects/${projectId}/comments`),

  // Add comment
  addComment: (projectId, comment) =>
    api.post(`/collaboration/projects/${projectId}/comments`, comment),

  // Get activity log
  getActivityLog: (projectId) =>
    api.get(`/collaboration/projects/${projectId}/activity`),
};

// Claude AI Assistant API
export const claudeAPI = {
  // Ask Claude a question
  askQuestion: (question, context = null) =>
    api.post("/ai/ask", { question, context }),

  // Get molecule generation thinking process
  getMoleculeThinking: (requestId) =>
    api.get(`/ai/molecule-thinking/${requestId}`),

  // Generate molecule
  generateMolecule: (requirements) =>
    api.post("/ai/generate-molecule", requirements),

  // Analyze molecule
  analyzeMolecule: (smiles) => api.post("/ai/analyze-molecule", { smiles }),

  // Compare molecules
  compareMolecules: (molecules) =>
    api.post("/ai/compare-molecules", { molecules }),

  // Analyze literature
  analyzeLiterature: (articles, focus) =>
    api.post("/ai/analyze-literature", { articles, focus }),

  // Get detailed thinking process
  getThinkingProcess: (requestId) =>
    api.get(`/ai/thinking-process/${requestId}`),

  // Continue chat conversation
  continueChat: (messages, context = null) =>
    api.post("/ai/chat", { messages, context }),
};

export default {
  auth: authAPI,
  drugDesign: drugDesignAPI,
  simulation: simulationAPI,
  literature: literatureAPI,
  chemicalData: chemicalDataAPI,
  regulatory: regulatoryAPI,
  collaboration: collaborationAPI,
  claude: claudeAPI,
};
