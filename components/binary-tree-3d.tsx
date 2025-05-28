"use client"

import { useRef, useState, useEffect, useMemo } from "react"
import { Canvas, useFrame } from "@react-three/fiber"
import { Text, OrbitControls, Environment, Float } from "@react-three/drei"
import * as THREE from "three"
import type { JSX } from "react/jsx-runtime"

// Updated interface for general knowledge visualization
export interface ConceptNode {
  id: string | number;
  value: string | number;
  label?: string;
  type?: "main" | "concept" | "detail" | "connection";
  children?: ConceptNode[];
  position?: [number, number, number];
  color?: string;
  size?: number;
  importance?: number; // 0-1 scale for visual emphasis
}

interface TransformedConceptNode {
  id: string | number;
  value: string | number;
  label: string;
  position: [number, number, number];
  children: TransformedConceptNode[];
  isHighlighted?: boolean;
  isVisiting?: boolean;
  color?: string;
  size: number;
  depth: number;
  parentPosition?: [number, number, number];
  type: "main" | "concept" | "detail" | "connection";
  importance: number;
}

interface ConceptNodeProps {
  node: TransformedConceptNode;
  onNodeClick?: (value: string | number) => void;
}

function ConceptNodeComponent({ node, onNodeClick }: ConceptNodeProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [hovered, setHovered] = useState(false)

  useFrame((state) => {
    if (meshRef.current) {
      if (node.isVisiting) {
        meshRef.current.scale.setScalar(1.3 + Math.sin(state.clock.elapsedTime * 6) * 0.15)
      } else {
        meshRef.current.scale.setScalar(hovered ? 1.15 : 1)
      }
      
      // Gentle floating animation based on importance
      if (node.type === "main") {
        meshRef.current.position.y += Math.sin(state.clock.elapsedTime * 2) * 0.05
      }
    }
  })

  const getNodeColor = () => {
    if (node.color) return node.color;
    switch (node.type) {
      case "main": return "#20808D";
      case "concept": return "#3B82F6";
      case "detail": return "#64748B";
      case "connection": return "#10B981";
      default: return "#64748B";
    }
  }

  const baseColor = getNodeColor()
  const nodeColor = node.isVisiting ? "#10B981" : node.isHighlighted ? "#20808D" : hovered ? "#3B82F6" : baseColor
  const emissiveColor = node.isVisiting ? "#065F46" : node.isHighlighted ? "#0F4C5C" : "#000000"

  const getGeometry = () => {
    switch (node.type) {
      case "main":
        return <sphereGeometry args={[node.size * 1.2, 32, 32]} />
      case "concept":
        return <boxGeometry args={[node.size, node.size, node.size]} />
      case "detail":
        return <sphereGeometry args={[node.size * 0.8, 16, 16]} />
      case "connection":
        return <octahedronGeometry args={[node.size * 0.9]} />
      default:
        return <sphereGeometry args={[node.size, 16, 16]} />
    }
  }

  return (
    <group position={node.position}>
      <Float speed={1.5} rotationIntensity={0.2} floatIntensity={0.3}>
        <mesh
          ref={meshRef}
          onClick={() => onNodeClick?.(node.value)}
          onPointerOver={() => setHovered(true)}
          onPointerOut={() => setHovered(false)}
        >
          {getGeometry()}
          <meshStandardMaterial
            color={nodeColor}
            emissive={emissiveColor}
            emissiveIntensity={node.isVisiting ? 0.3 : 0.1}
            roughness={0.2}
            metalness={0.8}
            opacity={0.8 + node.importance * 0.2}
            transparent
          />
        </mesh>
        <Text
          position={[0, 0, node.size + 0.1]}
          fontSize={Math.max(0.15, node.size * 0.4)}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
          maxWidth={3}
          textAlign="center"
        >
          {String(node.label)}
        </Text>
      </Float>
    </group>
  )
}

function ConceptConnections({ nodes }: { nodes: TransformedConceptNode[] }) {
  const connections: Array<{ start: [number, number, number]; end: [number, number, number]; key: string; type: string }> = []

  nodes.forEach(node => {
    if (node.parentPosition) {
      connections.push({
        start: node.parentPosition,
        end: node.position,
        key: `conn-${node.id}-parent`,
        type: node.type
      })
    }
  })

  return (
    <>
      {connections.map((conn) => {
        const start = new THREE.Vector3(...conn.start)
        const end = new THREE.Vector3(...conn.end)
        const direction = end.clone().sub(start)
        const length = direction.length()
        const midpoint = start.clone().add(direction.clone().multiplyScalar(0.5))
        const orientation = new THREE.Matrix4()
        const up = new THREE.Vector3(0, 1, 0)
        orientation.lookAt(start, end, up)
        orientation.multiply(new THREE.Matrix4().makeRotationX(Math.PI / 2))

        const getConnectionColor = () => {
          switch (conn.type) {
            case "main": return "#20808D";
            case "concept": return "#3B82F6";
            case "detail": return "#64748B";
            case "connection": return "#10B981";
            default: return "#20808D";
          }
        }

        return (
          <mesh
            key={conn.key}
            position={midpoint.toArray()}
            quaternion={new THREE.Quaternion().setFromRotationMatrix(orientation)}
          >
            <cylinderGeometry args={[0.02, 0.02, length, 8]} />
            <meshStandardMaterial
              color={getConnectionColor()}
              opacity={0.6}
              transparent
              emissive="#0F4C5C"
              emissiveIntensity={0.1}
            />
          </mesh>
        )
      })}
    </>
  )
}

function ParticleField() {
  const particlesRef = useRef<THREE.Points>(null)

  useFrame((state) => {
    if (particlesRef.current) {
      particlesRef.current.rotation.y = state.clock.elapsedTime * 0.05
    }
  })

  const particleCount = 150
  const positions = new Float32Array(particleCount * 3)

  for (let i = 0; i < particleCount; i++) {
    positions[i * 3] = (Math.random() - 0.5) * 25
    positions[i * 3 + 1] = (Math.random() - 0.5) * 25
    positions[i * 3 + 2] = (Math.random() - 0.5) * 25
  }

  return (
    <points ref={particlesRef}>
      <bufferGeometry>
        <bufferAttribute 
          attach="attributes-position" 
          args={[positions, 3]}
          array={positions} 
          count={particleCount} 
          itemSize={3}
        />
      </bufferGeometry>
      <pointsMaterial color="#20808D" size={0.02} opacity={0.4} transparent />
    </points>
  )
}

// Function to generate knowledge visualization from query data
function generateKnowledgeVisualization(data: any): ConceptNode {
  if (!data) {
    return {
      id: "default",
      value: "Knowledge",
      label: "Knowledge Visualization",
      type: "main",
      children: [
        {
          id: "concept1",
          value: "Concepts",
          label: "Core Concepts",
          type: "concept",
          children: []
        },
        {
          id: "detail1",
          value: "Details",
          label: "Key Details",
          type: "detail",
          children: []
        }
      ]
    }
  }

  // Extract key concepts from the data
  const query = data.query || "Learning Topic"
  const mode = data.mode || "basic"
  const content = data.content || ""
  
  // Generate concepts based on the query and content
  const mainConcepts = extractConcepts(query, content, mode)
  
  return {
    id: "main",
    value: query,
    label: query.length > 20 ? query.substring(0, 20) + "..." : query,
    type: "main",
    importance: 1,
    children: mainConcepts
  }
}

function extractConcepts(query: string, content: string, mode: string): ConceptNode[] {
  const concepts: ConceptNode[] = []
  
  // Basic concept extraction based on query keywords
  const queryWords = query.toLowerCase().split(/\s+/).filter(word => word.length > 3)
  const importantWords = ['algorithm', 'theory', 'principle', 'method', 'concept', 'system', 'process', 'structure']
  
  queryWords.forEach((word, index) => {
    if (importantWords.some(imp => word.includes(imp) || imp.includes(word))) {
      concepts.push({
        id: `concept_${index}`,
        value: word,
        label: word.charAt(0).toUpperCase() + word.slice(1),
        type: "concept",
        importance: 0.8,
        children: []
      })
    }
  })
  
  // Add mode-specific concepts
  if (mode === "reasoning") {
    concepts.push({
      id: "reasoning",
      value: "Analysis",
      label: "Logical Analysis",
      type: "connection",
      importance: 0.9,
      children: []
    })
  } else if (mode === "deep-research") {
    concepts.push({
      id: "research",
      value: "Research",
      label: "Deep Research",
      type: "connection",
      importance: 0.9,
      children: []
    })
  }
  
  // Add some detail nodes
  concepts.push({
    id: "applications",
    value: "Applications",
    label: "Applications",
    type: "detail",
    importance: 0.6,
    children: []
  })
  
  concepts.push({
    id: "examples",
    value: "Examples",
    label: "Examples",
    type: "detail",
    importance: 0.5,
    children: []
  })
  
  return concepts.slice(0, 6) // Limit to 6 concepts for clean visualization
}

// Layout function for knowledge visualization
function layoutKnowledge(node: ConceptNode, depth = 0, xOffset = 0, yOffset = 0, parentPos?: [number, number, number]): TransformedConceptNode {
  const ySpacing = -2
  const xSpacing = 3
  const radiusSpacing = 4

  // Calculate position based on node type and depth
  let position: [number, number, number]
  
  if (node.position) {
    position = node.position
  } else if (depth === 0) {
    // Main node at center
    position = [0, 0, 0]
  } else if (depth === 1) {
    // Arrange first level in a circle around the main node
    const angle = (xOffset / (node.children?.length || 1)) * 2 * Math.PI
    const radius = radiusSpacing
    position = [
      Math.cos(angle) * radius,
      Math.sin(angle) * radius * 0.5,
      Math.sin(angle) * radius * 0.3
    ]
      } else {
    // Subsequent levels spread out more
    position = [xOffset, yOffset + depth * ySpacing, depth * 0.5]
  }

  const size = node.size || (
    node.type === "main" ? 0.5 :
    node.type === "concept" ? 0.35 :
    node.type === "connection" ? 0.3 :
    0.25
  )

  const transformedNode: TransformedConceptNode = {
    ...node,
    label: node.label || String(node.value),
    position,
    children: [],
    depth,
    parentPosition: parentPos,
    size,
    type: node.type || "concept",
    importance: node.importance || 0.5
  }

  if (node.children && node.children.length > 0) {
    const numChildren = node.children.length
    
    transformedNode.children = node.children.map((child, index) => {
      let childX = xOffset
      let childY = yOffset
      
      if (depth === 0) {
        // First level children arranged in circle
        childX = index
      } else {
        // Subsequent levels spread horizontally
        childX = xOffset + (index - (numChildren - 1) / 2) * xSpacing
      }
      
      return layoutKnowledge(child, depth + 1, childX, childY, position)
    })
  }
  
  return transformedNode
}

function flattenKnowledge(node: TransformedConceptNode): TransformedConceptNode[] {
  let nodes: TransformedConceptNode[] = [node]
  node.children.forEach(child => {
    nodes = nodes.concat(flattenKnowledge(child))
  })
  return nodes
}

interface BinaryTree3DProps {
  visualizationData?: any
  title?: string
}

export function BinaryTree3D({ visualizationData, title = "Knowledge Visualization" }: BinaryTree3DProps) {
  const [currentOperation, setCurrentOperation] = useState<string>("Initializing...")
  const [visitingPath, setVisitingPath] = useState<(string | number)[]>([])

  const knowledgeTree = useMemo(() => {
    const conceptData = generateKnowledgeVisualization(visualizationData)
    return layoutKnowledge(conceptData)
  }, [visualizationData])

  const allNodes: TransformedConceptNode[] = useMemo(() => {
    return flattenKnowledge(knowledgeTree)
  }, [knowledgeTree])

  useEffect(() => {
    if (allNodes.length > 0) {
      if (visualizationData?.query) {
        setCurrentOperation(`Visualizing: ${visualizationData.query}`)
      } else {
        setCurrentOperation("Ready for exploration")
      }
    } else {
      setCurrentOperation("Initializing visualization...")
    }
  }, [allNodes, visualizationData])

  const handleNodeClick = (value: string | number) => {
    console.log("Node clicked:", value)
    // Could trigger additional actions like expanding details or searching
  }

  // Apply highlighting based on visitingPath
  const nodesToRender = allNodes.map(node => ({
    ...node,
    isVisiting: visitingPath.includes(node.id),
    isHighlighted: visitingPath.includes(node.id),
  }))

  return (
    <div className="w-full h-full relative overflow-hidden rounded-2xl border border-[#20808D]/30">
      <div className="absolute top-6 left-6 z-10 bg-gradient-to-br from-[#20808D]/20 to-[#20808D]/10 backdrop-blur-xl rounded-2xl p-4 border border-[#20808D]/30 shadow-2xl max-w-xs">
        <div className="text-[#20808D] font-bold text-sm mb-2 flex items-center gap-2">
          <div className="w-2 h-2 bg-[#20808D] rounded-full animate-pulse"></div>
          {title}
        </div>
        <div className="text-[#FBFAF4] text-xs mb-1">Status: {currentOperation}</div>
        <div className="text-[#20808D]/80 text-xs">
          Mode: {visualizationData?.mode || "basic"} • Nodes: {allNodes.length}
        </div>
        {visitingPath.length > 0 && (
          <div className="text-[#20808D] text-xs mt-1">Path: {visitingPath.join(" → ")}</div>
        )}
      </div>

      <Canvas camera={{ position: [0, 3, 12], fov: 45 }}>
        <Environment preset="night" />
        <ambientLight intensity={0.4} />
        <pointLight position={[10, 10, 10]} intensity={1.5} color="#20808D" />
        <pointLight position={[-10, -10, -10]} intensity={0.8} color="#FBFAF4" />
        <spotLight position={[0, 15, 0]} intensity={1.2} color="#20808D" angle={0.3} penumbra={1} />

        <ParticleField />
        <ConceptConnections nodes={nodesToRender} />
        {nodesToRender.map(node => (
          <ConceptNodeComponent 
            key={node.id} 
            node={node} 
            onNodeClick={handleNodeClick}
          />
        ))}

        <OrbitControls
          enablePan={true}
          enableZoom={true}
          maxDistance={50}
          minDistance={3}
          autoRotate
          autoRotateSpeed={0.5}
          maxPolarAngle={Math.PI / 1.6}
          minPolarAngle={Math.PI / 6}
        />
      </Canvas>
    </div>
  )
}
