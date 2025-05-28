import type { APIResponse } from "../types/monday"

export const sampleResponses: APIResponse[] = [
  {
    id: "1",
    mode: "basic",
    query: "Explain machine learning algorithms",
    content: `Machine learning algorithms are computational methods that enable systems to automatically learn and improve from experience without being explicitly programmed for every scenario. These algorithms can be broadly categorized into three main types: supervised learning, unsupervised learning, and reinforcement learning.

Supervised learning algorithms learn from labeled training data to make predictions on new, unseen data. Popular examples include linear regression for predicting continuous values, decision trees for classification tasks, and neural networks for complex pattern recognition. These algorithms require a dataset where both input features and desired outputs are known.

Unsupervised learning algorithms find hidden patterns in data without labeled examples. Clustering algorithms like K-means group similar data points together, while dimensionality reduction techniques like Principal Component Analysis (PCA) simplify complex datasets by identifying the most important features.

Reinforcement learning algorithms learn through interaction with an environment, receiving rewards or penalties for their actions. This approach has been particularly successful in game-playing AI, robotics, and autonomous systems where the algorithm must learn optimal strategies through trial and error.

The choice of algorithm depends on factors such as the type of problem, available data, computational resources, and desired accuracy. Modern machine learning often combines multiple algorithms in ensemble methods to achieve better performance than any single algorithm could provide.`,
    timestamp: new Date(),
  },
  {
    id: "2",
    mode: "reasoning",
    query: "Think about the relationship between quantum computing and cryptography",
    content: `The relationship between quantum computing and cryptography represents one of the most significant paradigm shifts in information security. Let me walk through the logical connections and implications.`,
    reasoning_steps: [
      {
        id: "r1",
        step: "Current Cryptographic Foundation",
        reasoning:
          "Modern cryptography relies on mathematical problems that are computationally difficult for classical computers, such as factoring large integers (RSA) and solving discrete logarithm problems (ECC).",
        confidence: 0.95,
      },
      {
        id: "r2",
        step: "Quantum Computational Advantage",
        reasoning:
          "Quantum computers can solve these specific mathematical problems exponentially faster using algorithms like Shor's algorithm, which can factor large integers in polynomial time.",
        confidence: 0.92,
      },
      {
        id: "r3",
        step: "Cryptographic Vulnerability",
        reasoning:
          "This quantum advantage means that sufficiently powerful quantum computers could break RSA, ECC, and other widely-used public key cryptographic systems that secure internet communications.",
        confidence: 0.88,
      },
      {
        id: "r4",
        step: "Timeline Considerations",
        reasoning:
          "While current quantum computers are not yet powerful enough to break real-world cryptographic systems, experts estimate that cryptographically relevant quantum computers may emerge within 10-30 years.",
        confidence: 0.75,
      },
      {
        id: "r5",
        step: "Post-Quantum Solutions",
        reasoning:
          "The cryptographic community is developing quantum-resistant algorithms based on mathematical problems that even quantum computers cannot solve efficiently, such as lattice-based and hash-based cryptography.",
        confidence: 0.85,
      },
    ],
    confidence: 0.87,
    timestamp: new Date(),
  },
  {
    id: "3",
    mode: "deep-research",
    query: "Deep dive on neural network architectures",
    content: `Neural network architectures have evolved dramatically since the perceptron model of the 1950s, with modern deep learning representing a convergence of computational power, algorithmic innovation, and vast datasets. This comprehensive analysis examines the fundamental principles, architectural innovations, and practical applications that define contemporary neural networks.

The foundational concept of artificial neural networks draws inspiration from biological neural systems, where interconnected nodes (neurons) process and transmit information through weighted connections. However, modern architectures have diverged significantly from their biological inspiration, optimizing for computational efficiency and specific task performance rather than biological plausibility.

Feedforward neural networks represent the simplest architecture, where information flows in one direction from input to output layers. These networks excel at function approximation and classification tasks, with the universal approximation theorem proving that sufficiently wide networks can approximate any continuous function. The depth versus width trade-off in these networks reveals that deeper architectures often achieve better generalization with fewer parameters.

Convolutional Neural Networks (CNNs) revolutionized computer vision by incorporating spatial locality and translation invariance through convolutional layers. The hierarchical feature extraction in CNNs mirrors the visual cortex's processing, with early layers detecting edges and textures while deeper layers recognize complex objects and scenes. Architectural innovations like ResNet's skip connections solved the vanishing gradient problem, enabling networks with hundreds of layers.

Recurrent Neural Networks (RNNs) and their variants (LSTM, GRU) address sequential data processing by maintaining internal state across time steps. These architectures excel at natural language processing, time series analysis, and any domain where temporal dependencies matter. The attention mechanism, initially developed for RNNs, has become a fundamental building block across architectures.

Transformer architectures represent the current state-of-the-art for many tasks, replacing recurrence with self-attention mechanisms that can process sequences in parallel. The multi-head attention mechanism allows the model to focus on different aspects of the input simultaneously, while positional encodings provide sequence order information without recurrence.

Graph Neural Networks (GNNs) extend deep learning to non-Euclidean data structures, enabling applications in social networks, molecular modeling, and knowledge graphs. These architectures aggregate information from neighboring nodes through message passing, with variants like Graph Convolutional Networks and GraphSAGE offering different approaches to neighborhood aggregation.

Generative architectures like Variational Autoencoders (VAEs) and Generative Adversarial Networks (GANs) focus on learning data distributions rather than discriminative tasks. VAEs provide a principled probabilistic framework for generation, while GANs use adversarial training to achieve remarkable generation quality in domains like image synthesis.

The emergence of foundation models and large language models represents a paradigm shift toward general-purpose architectures trained on massive datasets. These models demonstrate emergent capabilities that weren't explicitly trained for, suggesting that scale and architecture design can lead to qualitatively different behaviors.

Modern architectural design principles emphasize modularity, with components like attention blocks, normalization layers, and activation functions becoming standardized building blocks. The trend toward automated architecture search (NAS) uses machine learning to discover optimal architectures, though human insight remains crucial for novel architectural innovations.

Training dynamics and optimization landscapes vary significantly across architectures, with some requiring careful initialization schemes, specific learning rate schedules, or regularization techniques. The interplay between architecture design and optimization algorithms continues to drive performance improvements across domains.`,
    sources: [
      {
        id: "s1",
        title: "Deep Learning (Goodfellow, Bengio, Courville)",
        url: "https://www.deeplearningbook.org/",
        relevance: 0.95,
        excerpt: "Comprehensive treatment of neural network architectures and their mathematical foundations.",
      },
      {
        id: "s2",
        title: "Attention Is All You Need (Vaswani et al.)",
        url: "https://arxiv.org/abs/1706.03762",
        relevance: 0.92,
        excerpt: "Seminal paper introducing the Transformer architecture that revolutionized NLP.",
      },
      {
        id: "s3",
        title: "ResNet: Deep Residual Learning (He et al.)",
        url: "https://arxiv.org/abs/1512.03385",
        relevance: 0.88,
        excerpt: "Introduction of skip connections enabling very deep neural networks.",
      },
      {
        id: "s4",
        title: "Graph Neural Networks: A Review (Wu et al.)",
        url: "https://arxiv.org/abs/1901.00596",
        relevance: 0.85,
        excerpt: "Comprehensive survey of graph neural network architectures and applications.",
      },
    ],
    timestamp: new Date(),
  },
]
