"use client"

import { useState, useEffect, useRef, useMemo } from "react"
import { Canvas, useFrame, useLoader } from "@react-three/fiber"
import { Text, OrbitControls, Environment, Float, Box, Sphere, Cylinder, Cone, Center } from "@react-three/drei"
import * as THREE from "three"
import { BinaryTree3D } from "./binary-tree-3d"

interface AdaptiveVisualizationPanelProps {
  mode: "basic" | "reasoning" | "deep-research"
  query?: string
  visualizationData?: any
  isThinking?: boolean
  progressUpdates?: string[]
  sources?: any[]
  reasoning?: any[]
  citations?: any[]
}

// Enhanced 3D object generator for basic mode queries
function generate3DVisualization(query: string) {
  const queryLower = query.toLowerCase()
  
  // Coffee-related queries
  if (queryLower.includes("coffee") || queryLower.includes("espresso") || queryLower.includes("latte")) {
    return {
      type: "coffee",
      objects: [
        { type: "coffee-cup", brand: "Starbucks", position: [-2, 0, 0] as [number, number, number] },
        { type: "coffee-cup", brand: "Dunkin", position: [0, 0, 0] as [number, number, number] },
        { type: "coffee-cup", brand: "Local Cafe", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Tree/nature queries
  if (queryLower.includes("tree") || queryLower.includes("forest") || queryLower.includes("nature")) {
    return {
      type: "nature",
      objects: [
        { type: "tree", variant: "oak", position: [-3, 0, 0] as [number, number, number] },
        { type: "tree", variant: "pine", position: [0, 0, 0] as [number, number, number] },
        { type: "tree", variant: "maple", position: [3, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Technology queries
  if (queryLower.includes("computer") || queryLower.includes("technology") || queryLower.includes("coding")) {
    return {
      type: "technology",
      objects: [
        { type: "laptop", position: [-2, 0, 0] as [number, number, number] },
        { type: "phone", position: [0, 0, 0] as [number, number, number] },
        { type: "server", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Space/astronomy queries
  if (queryLower.includes("space") || queryLower.includes("planet") || queryLower.includes("star")) {
    return {
      type: "space",
      objects: [
        { type: "planet", name: "Earth", position: [-3, 0, 0] as [number, number, number] },
        { type: "planet", name: "Mars", position: [0, 0, 0] as [number, number, number] },
        { type: "star", name: "Sun", position: [3, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Food queries
  if (queryLower.includes("food") || queryLower.includes("pizza") || queryLower.includes("burger") || queryLower.includes("fruit")) {
    return {
      type: "food",
      objects: [
        { type: "food", variant: "pizza", position: [-2, 0, 0] as [number, number, number] },
        { type: "food", variant: "burger", position: [0, 0, 0] as [number, number, number] },
        { type: "food", variant: "apple", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Animals/pets queries
  if (queryLower.includes("animal") || queryLower.includes("dog") || queryLower.includes("cat") || queryLower.includes("pet")) {
    return {
      type: "animals",
      objects: [
        { type: "animal", variant: "dog", position: [-2, 0, 0] as [number, number, number] },
        { type: "animal", variant: "cat", position: [0, 0, 0] as [number, number, number] },
        { type: "animal", variant: "bird", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Sports queries
  if (queryLower.includes("sport") || queryLower.includes("ball") || queryLower.includes("soccer") || queryLower.includes("basketball")) {
    return {
      type: "sports",
      objects: [
        { type: "ball", variant: "soccer", position: [-2, 0, 0] as [number, number, number] },
        { type: "ball", variant: "basketball", position: [0, 0, 0] as [number, number, number] },
        { type: "ball", variant: "tennis", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Music queries
  if (queryLower.includes("music") || queryLower.includes("instrument") || queryLower.includes("guitar") || queryLower.includes("piano")) {
    return {
      type: "music",
      objects: [
        { type: "instrument", variant: "guitar", position: [-2, 0, 0] as [number, number, number] },
        { type: "instrument", variant: "piano", position: [0, 0, 0] as [number, number, number] },
        { type: "instrument", variant: "drums", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Mathematics/geometry queries
  if (queryLower.includes("math") || queryLower.includes("geometry") || queryLower.includes("shape") || queryLower.includes("equation")) {
    return {
      type: "mathematics",
      objects: [
        { type: "shape", variant: "cube", position: [-2, 0, 0] as [number, number, number] },
        { type: "shape", variant: "sphere", position: [0, 0, 0] as [number, number, number] },
        { type: "shape", variant: "pyramid", position: [2, 0, 0] as [number, number, number] }
      ]
    }
  }
  
  // Default to abstract knowledge visualization
  return {
    type: "abstract",
    objects: []
  }
}

// Enhanced Coffee cup component with more detail
function CoffeeCup({ brand, position }: { brand: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  const steamRef = useRef<THREE.Points>(null)
  const [hovered, setHovered] = useState(false)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.5
      meshRef.current.scale.setScalar(hovered ? 1.1 : 1)
    }
    
    // Animate steam particles
    if (steamRef.current) {
      const positions = steamRef.current.geometry.attributes.position.array as Float32Array
      for (let i = 0; i < positions.length; i += 3) {
        positions[i + 1] += 0.01 // Move steam upward
        if (positions[i + 1] > 2) {
          positions[i + 1] = 0.8 // Reset to cup top
          positions[i] = (Math.random() - 0.5) * 0.4
          positions[i + 2] = (Math.random() - 0.5) * 0.4
        }
      }
      steamRef.current.geometry.attributes.position.needsUpdate = true
    }
  })
  
  const getBrandColor = () => {
    switch (brand) {
      case "Starbucks": return "#00704A"
      case "Dunkin": return "#FF6E1B"
      default: return "#8B4513"
    }
  }
  
  // Create steam particles
  const steamParticles = useMemo(() => {
    const particles = new Float32Array(150) // 50 particles * 3 coordinates
    for (let i = 0; i < particles.length; i += 3) {
      particles[i] = (Math.random() - 0.5) * 0.4
      particles[i + 1] = 0.8 + Math.random() * 1.2
      particles[i + 2] = (Math.random() - 0.5) * 0.4
    }
    return particles
  }, [])
  
  return (
    <Float speed={2} rotationIntensity={0.2} floatIntensity={0.3}>
      <group 
        ref={meshRef} 
        position={position}
        onPointerOver={() => setHovered(true)}
        onPointerOut={() => setHovered(false)}
      >
        {/* Cup body with gradient */}
        <mesh position={[0, 0, 0]} castShadow receiveShadow>
          <cylinderGeometry args={[0.6, 0.5, 1.5, 64, 1, false]} />
          <meshPhysicalMaterial 
            color={getBrandColor()} 
            metalness={0.2} 
            roughness={0.3}
            clearcoat={0.8}
            clearcoatRoughness={0.2}
            reflectivity={0.5}
          />
        </mesh>
        
        {/* Cup rim */}
        <mesh position={[0, 0.75, 0]}>
          <torusGeometry args={[0.6, 0.05, 16, 32]} />
          <meshStandardMaterial color="#CCCCCC" metalness={0.8} roughness={0.2} />
        </mesh>
        
        {/* Handle with better geometry */}
        <mesh position={[0.7, 0, 0]} rotation={[0, 0, Math.PI / 2]} castShadow>
          <torusGeometry args={[0.3, 0.08, 16, 32, Math.PI]} />
          <meshPhysicalMaterial 
            color={getBrandColor()} 
            metalness={0.2} 
            roughness={0.3}
            clearcoat={0.8}
          />
        </mesh>
        
        {/* Coffee inside with foam */}
        <mesh position={[0, 0.6, 0]}>
          <cylinderGeometry args={[0.55, 0.55, 0.1, 32]} />
          <meshStandardMaterial color="#3E2723" metalness={0.1} roughness={0.8} />
        </mesh>
        
        {/* Foam layer */}
        <mesh position={[0, 0.65, 0]}>
          <cylinderGeometry args={[0.54, 0.54, 0.05, 32]} />
          <meshStandardMaterial color="#F5DEB3" metalness={0} roughness={0.9} />
        </mesh>
        
        {/* Steam particles */}
        <points ref={steamRef}>
          <bufferGeometry>
            <bufferAttribute
              attach="attributes-position"
              count={50}
              array={steamParticles}
              itemSize={3}
              args={[steamParticles, 3]}
            />
          </bufferGeometry>
          <pointsMaterial
            size={0.05}
            color="#FFFFFF"
            transparent
            opacity={0.4}
            sizeAttenuation
            blending={THREE.AdditiveBlending}
          />
        </points>
        
        {/* Brand logo/text with background */}
        <group position={[0, 0, 0.61]}>
          <mesh>
            <planeGeometry args={[0.8, 0.3]} />
            <meshBasicMaterial color="#FFFFFF" opacity={0.9} transparent />
          </mesh>
          <Text
            position={[0, 0, 0.01]}
            fontSize={0.2}
            color={getBrandColor()}
            anchorX="center"
            anchorY="middle"
            font="/fonts/helvetica.woff"
          >
            {brand}
          </Text>
        </group>
        
        {/* Bottom text */}
        <Text
          position={[0, -1.2, 0]}
          fontSize={0.25}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {brand} Coffee
        </Text>
      </group>
    </Float>
  )
}

// Enhanced Tree component with leaves and branches
function Tree({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  const leavesRef = useRef<THREE.InstancedMesh>(null)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = Math.sin(state.clock.elapsedTime * 0.2) * 0.1
      
      // Animate leaves
      if (leavesRef.current && variant !== "pine") {
        const time = state.clock.elapsedTime
        for (let i = 0; i < 100; i++) {
          const matrix = new THREE.Matrix4()
          const x = (Math.random() - 0.5) * 3
          const y = 1 + Math.random() * 2
          const z = (Math.random() - 0.5) * 3
          const scale = 0.1 + Math.random() * 0.2
          
          matrix.makeTranslation(
            x + Math.sin(time + i) * 0.1,
            y + Math.sin(time * 0.5 + i) * 0.05,
            z + Math.cos(time + i) * 0.1
          )
          matrix.scale(new THREE.Vector3(scale, scale, scale))
          
          leavesRef.current.setMatrixAt(i, matrix)
        }
        leavesRef.current.instanceMatrix.needsUpdate = true
      }
    }
  })
  
  const getTreeColor = () => {
    switch (variant) {
      case "oak": return "#228B22"
      case "pine": return "#00563F"
      case "maple": return "#FF8C00"
      default: return "#228B22"
    }
  }
  
  return (
    <group ref={meshRef} position={position}>
      {/* Trunk with bark texture */}
      <mesh position={[0, -1, 0]} castShadow receiveShadow>
        <cylinderGeometry args={[0.3, 0.4, 2, 12, 4]} />
        <meshStandardMaterial 
          color="#8B4513" 
          roughness={0.9}
          normalScale={new THREE.Vector2(2, 2)}
        />
      </mesh>
      
      {/* Branches */}
      {[...Array(5)].map((_, i) => (
        <mesh 
          key={i}
          position={[
            Math.sin(i * 1.2) * 0.5,
            0.5 + i * 0.3,
            Math.cos(i * 1.2) * 0.5
          ]}
          rotation={[0, i * 0.5, Math.PI / 6]}
          castShadow
        >
          <cylinderGeometry args={[0.05, 0.02, 0.8, 8]} />
          <meshStandardMaterial color="#654321" roughness={0.8} />
        </mesh>
      ))}
      
      {/* Leaves or pine needles */}
      {variant === "pine" ? (
        <>
          {/* Pine cone shape */}
          <mesh position={[0, 1, 0]} castShadow receiveShadow>
            <coneGeometry args={[1.5, 3, 12]} />
            <meshStandardMaterial 
              color={getTreeColor()} 
              roughness={0.7}
              metalness={0.1}
            />
          </mesh>
          <mesh position={[0, 2.2, 0]} castShadow receiveShadow>
            <coneGeometry args={[1.2, 2, 12]} />
            <meshStandardMaterial 
              color={getTreeColor()} 
              roughness={0.7}
              metalness={0.1}
            />
          </mesh>
        </>
      ) : (
        <>
          {/* Leafy canopy */}
          <mesh position={[0, 1.5, 0]} castShadow receiveShadow>
            <sphereGeometry args={[1.8, 16, 16]} />
            <meshStandardMaterial 
              color={getTreeColor()} 
              roughness={0.8}
              metalness={0.1}
              opacity={0.9}
              transparent
            />
          </mesh>
          
          {/* Individual leaves */}
          <instancedMesh ref={leavesRef} args={[undefined, undefined, 100]} castShadow>
            <planeGeometry args={[0.2, 0.3]} />
            <meshStandardMaterial 
              color={getTreeColor()}
              side={THREE.DoubleSide}
              opacity={0.8}
              transparent
            />
          </instancedMesh>
        </>
      )}
      
      {/* Ground/shadow plane */}
      <mesh position={[0, -2, 0]} rotation={[-Math.PI / 2, 0, 0]} receiveShadow>
        <circleGeometry args={[2, 32]} />
        <meshStandardMaterial color="#3B5323" opacity={0.5} transparent />
      </mesh>
      
      <Text
        position={[0, -2.5, 0]}
        fontSize={0.3}
        color="#FBFAF4"
        anchorX="center"
        anchorY="middle"
      >
        {variant.charAt(0).toUpperCase() + variant.slice(1)} Tree
      </Text>
    </group>
  )
}

// Planet component
function Planet({ name, position }: { name: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Mesh>(null)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.3
      meshRef.current.rotation.x = state.clock.elapsedTime * 0.1
    }
  })
  
  const getPlanetColor = () => {
    switch (name) {
      case "Earth": return "#4169E1"
      case "Mars": return "#CD5C5C"
      case "Sun": return "#FDB813"
      default: return "#808080"
    }
  }
  
  const size = name === "Sun" ? 2 : 1
  
  return (
    <Float speed={1.5} rotationIntensity={0.3} floatIntensity={0.5}>
      <group position={position}>
        <mesh ref={meshRef}>
          <sphereGeometry args={[size, 32, 32]} />
          <meshStandardMaterial 
            color={getPlanetColor()} 
            metalness={0.3} 
            roughness={0.5}
            emissive={name === "Sun" ? "#FDB813" : "#000000"}
            emissiveIntensity={name === "Sun" ? 0.5 : 0}
          />
        </mesh>
        
        {name === "Earth" && (
          <mesh position={[0, 0, 0.01]}>
            <sphereGeometry args={[1.01, 32, 32]} />
            <meshStandardMaterial 
              color="#90EE90" 
              transparent 
              opacity={0.3}
              metalness={0.1} 
              roughness={0.7}
            />
          </mesh>
        )}
        
        <Text
          position={[0, -size - 0.5, 0]}
          fontSize={0.4}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {name}
        </Text>
      </group>
    </Float>
  )
}

// Food component
function Food({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  const [hovered, setHovered] = useState(false)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.3
      meshRef.current.scale.setScalar(hovered ? 1.1 : 1)
    }
  })
  
  const getFoodColor = () => {
    switch (variant) {
      case "pizza": return "#FF6347"
      case "burger": return "#8B4513"
      case "apple": return "#DC143C"
      default: return "#FFA500"
    }
  }
  
  return (
    <Float speed={1.5} rotationIntensity={0.2} floatIntensity={0.3}>
      <group 
        ref={meshRef} 
        position={position}
        onPointerOver={() => setHovered(true)}
        onPointerOut={() => setHovered(false)}
      >
        {variant === "pizza" && (
          <mesh>
            <cylinderGeometry args={[1, 1, 0.2, 8]} />
            <meshStandardMaterial color={getFoodColor()} roughness={0.6} />
          </mesh>
        )}
        
        {variant === "burger" && (
          <>
            <mesh position={[0, -0.3, 0]}>
              <cylinderGeometry args={[0.8, 0.8, 0.2, 16]} />
              <meshStandardMaterial color="#8B4513" roughness={0.7} />
            </mesh>
            <mesh position={[0, 0, 0]}>
              <cylinderGeometry args={[0.7, 0.7, 0.3, 16]} />
              <meshStandardMaterial color="#8B0000" roughness={0.5} />
            </mesh>
            <mesh position={[0, 0.3, 0]}>
              <cylinderGeometry args={[0.8, 0.8, 0.2, 16]} />
              <meshStandardMaterial color="#F4A460" roughness={0.7} />
            </mesh>
          </>
        )}
        
        {variant === "apple" && (
          <mesh>
            <sphereGeometry args={[0.8, 16, 16]} />
            <meshStandardMaterial color={getFoodColor()} roughness={0.4} metalness={0.1} />
          </mesh>
        )}
        
        <Text
          position={[0, -1.5, 0]}
          fontSize={0.3}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {variant.charAt(0).toUpperCase() + variant.slice(1)}
        </Text>
      </group>
    </Float>
  )
}

// Animal component  
function Animal({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.position.y = Math.sin(state.clock.elapsedTime * 2) * 0.1
    }
  })
  
  const getAnimalColor = () => {
    switch (variant) {
      case "dog": return "#8B4513"
      case "cat": return "#FF8C00"
      case "bird": return "#1E90FF"
      default: return "#808080"
    }
  }
  
  return (
    <group ref={meshRef} position={position}>
      {/* Simple animal representations */}
      <mesh>
        <boxGeometry args={variant === "bird" ? [0.5, 0.5, 1] : [1, 0.8, 1.5]} />
        <meshStandardMaterial color={getAnimalColor()} roughness={0.7} />
      </mesh>
      
      {/* Head */}
      <mesh position={variant === "bird" ? [0, 0.3, 0.3] : [0, 0.5, 0.6]}>
        <sphereGeometry args={variant === "bird" ? [0.3, 8, 8] : [0.4, 8, 8]} />
        <meshStandardMaterial color={getAnimalColor()} roughness={0.7} />
      </mesh>
      
      <Text
        position={[0, -1.5, 0]}
        fontSize={0.3}
        color="#FBFAF4"
        anchorX="center"
        anchorY="middle"
      >
        {variant.charAt(0).toUpperCase() + variant.slice(1)}
      </Text>
    </group>
  )
}

// Sports ball component
function Ball({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Mesh>(null)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.x = state.clock.elapsedTime * 0.5
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.3
    }
  })
  
  const getBallColor = () => {
    switch (variant) {
      case "soccer": return "#FFFFFF"
      case "basketball": return "#FF6700"
      case "tennis": return "#CCFF00"
      default: return "#FF0000"
    }
  }
  
  return (
    <Float speed={2} rotationIntensity={0.5} floatIntensity={0.5}>
      <group position={position}>
        <mesh ref={meshRef}>
          <sphereGeometry args={[0.8, 32, 32]} />
          <meshStandardMaterial 
            color={getBallColor()} 
            roughness={variant === "soccer" ? 0.2 : 0.4}
            metalness={0.1}
          />
        </mesh>
        
        <Text
          position={[0, -1.5, 0]}
          fontSize={0.3}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {variant.charAt(0).toUpperCase() + variant.slice(1)}
        </Text>
      </group>
    </Float>
  )
}

// Musical instrument component
function Instrument({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  
  useFrame((state) => {
    if (meshRef.current && variant === "guitar") {
      meshRef.current.rotation.z = Math.sin(state.clock.elapsedTime) * 0.1
    }
  })
  
  return (
    <group ref={meshRef} position={position}>
      {variant === "guitar" && (
        <>
          <mesh>
            <boxGeometry args={[0.8, 2, 0.2]} />
            <meshStandardMaterial color="#8B4513" roughness={0.6} />
          </mesh>
          <mesh position={[0, -0.8, 0]}>
            <cylinderGeometry args={[0.5, 0.6, 0.2, 16]} />
            <meshStandardMaterial color="#654321" roughness={0.5} />
          </mesh>
        </>
      )}
      
      {variant === "piano" && (
        <mesh>
          <boxGeometry args={[2, 0.5, 1]} />
          <meshStandardMaterial color="#000000" roughness={0.2} metalness={0.8} />
        </mesh>
      )}
      
      {variant === "drums" && (
        <>
          <mesh position={[0, 0, 0]}>
            <cylinderGeometry args={[0.6, 0.6, 0.8, 16]} />
            <meshStandardMaterial color="#FF0000" roughness={0.5} />
          </mesh>
          <mesh position={[-0.8, 0.2, 0]}>
            <cylinderGeometry args={[0.4, 0.4, 0.6, 16]} />
            <meshStandardMaterial color="#0000FF" roughness={0.5} />
          </mesh>
        </>
      )}
      
      <Text
        position={[0, -1.5, 0]}
        fontSize={0.3}
        color="#FBFAF4"
        anchorX="center"
        anchorY="middle"
      >
        {variant.charAt(0).toUpperCase() + variant.slice(1)}
      </Text>
    </group>
  )
}

// Mathematical shape component
function Shape({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Mesh>(null)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.x = state.clock.elapsedTime * 0.5
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.3
    }
  })
  
  return (
    <Float speed={1} rotationIntensity={0.5} floatIntensity={0.3}>
      <group position={position}>
        <mesh ref={meshRef}>
          {variant === "cube" && <boxGeometry args={[1, 1, 1]} />}
          {variant === "sphere" && <sphereGeometry args={[0.8, 32, 32]} />}
          {variant === "pyramid" && <coneGeometry args={[0.8, 1.2, 4]} />}
          <meshStandardMaterial 
            color="#20808D" 
            metalness={0.5} 
            roughness={0.2}
            emissive="#20808D"
            emissiveIntensity={0.1}
          />
        </mesh>
        
        <Text
          position={[0, -1.5, 0]}
          fontSize={0.3}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {variant.charAt(0).toUpperCase() + variant.slice(1)}
        </Text>
      </group>
    </Float>
  )
}

// Enhanced Basic mode 3D visualization component with better lighting
function Basic3DVisualization({ query }: { query: string }) {
  const visualization = generate3DVisualization(query)
  
  return (
    <Canvas 
      camera={{ position: [0, 3, 10], fov: 50 }}
      shadows
      gl={{ antialias: true, alpha: true }}
    >
      <Environment preset="sunset" />
      <fog attach="fog" args={['#091717', 10, 30]} />
      
      {/* Enhanced lighting setup */}
      <ambientLight intensity={0.3} />
      <directionalLight
        position={[5, 10, 5]}
        intensity={1}
        castShadow
        shadow-mapSize={[2048, 2048]}
        shadow-camera-far={50}
        shadow-camera-left={-10}
        shadow-camera-right={10}
        shadow-camera-top={10}
        shadow-camera-bottom={-10}
      />
      <pointLight position={[10, 10, 10]} intensity={0.5} color="#20808D" />
      <pointLight position={[-10, -10, -10]} intensity={0.3} color="#FBFAF4" />
      <spotLight 
        position={[0, 15, 0]} 
        intensity={0.8} 
        color="#20808D" 
        angle={0.3} 
        penumbra={1}
        castShadow
      />
      
      {/* Ground plane */}
      <mesh rotation={[-Math.PI / 2, 0, 0]} position={[0, -2.5, 0]} receiveShadow>
        <planeGeometry args={[50, 50]} />
        <meshStandardMaterial color="#0A1A1A" roughness={0.9} />
      </mesh>
      
      {visualization.type === "coffee" && visualization.objects.map((obj: any, idx) => (
        <CoffeeCup key={idx} brand={obj.brand} position={obj.position} />
      ))}
      
      {visualization.type === "nature" && visualization.objects.map((obj: any, idx) => (
        <Tree key={idx} variant={obj.variant} position={obj.position} />
      ))}
      
      {visualization.type === "space" && visualization.objects.map((obj: any, idx) => (
        <Planet key={idx} name={obj.name} position={obj.position} />
      ))}
      
      {visualization.type === "food" && visualization.objects.map((obj: any, idx) => (
        <Food key={idx} variant={obj.variant} position={obj.position} />
      ))}
      
      {visualization.type === "animals" && visualization.objects.map((obj: any, idx) => (
        <Animal key={idx} variant={obj.variant} position={obj.position} />
      ))}
      
      {visualization.type === "sports" && visualization.objects.map((obj: any, idx) => (
        <Ball key={idx} variant={obj.variant} position={obj.position} />
      ))}
      
      {visualization.type === "music" && visualization.objects.map((obj: any, idx) => (
        <Instrument key={idx} variant={obj.variant} position={obj.position} />
      ))}
      
      {visualization.type === "mathematics" && visualization.objects.map((obj: any, idx) => (
        <Shape key={idx} variant={obj.variant} position={obj.position} />
      ))}
      
      <OrbitControls
        enablePan={true}
        enableZoom={true}
        maxDistance={20}
        minDistance={5}
        autoRotate
        autoRotateSpeed={0.5}
      />
    </Canvas>
  )
}

// Enhanced Progressive thinking/research display with progress bar
function ProgressiveThinkingDisplay({ 
  mode, 
  isThinking, 
  progressUpdates,
  sources,
  reasoning,
  citations,
  isComplete = false
}: {
  mode: "reasoning" | "deep-research"
  isThinking?: boolean
  progressUpdates?: string[]
  sources?: any[]
  reasoning?: any[]
  citations?: any[]
  isComplete?: boolean
}) {
  const [displayedText, setDisplayedText] = useState("")
  const [currentLine, setCurrentLine] = useState(0)
  const [progress, setProgress] = useState(0)
  const progressRef = useRef<HTMLDivElement>(null)
  
  const getHeaderText = () => {
    if (mode === "reasoning") {
      return isThinking ? "Monday is thinking..." : "Reasoning Process"
    } else {
      return isThinking ? "Monday is researching..." : "Research Process"
    }
  }
  
  const getHeaderIcon = () => {
    return mode === "reasoning" ? "🧠" : "🔍"
  }
  
  // Update progress based on updates
  useEffect(() => {
    if (isThinking && progressUpdates && progressUpdates.length > 0) {
      // Calculate progress based on number of updates
      const estimatedTotal = mode === "reasoning" ? 5 : 8 // More steps for research
      const currentProgress = Math.min((progressUpdates.length / estimatedTotal) * 100, 90)
      setProgress(currentProgress)
    }
    
    if (isComplete) {
      setProgress(100)
    }
  }, [progressUpdates, isThinking, isComplete, mode])
  
  useEffect(() => {
    if (progressUpdates && progressUpdates.length > 0) {
      const latestUpdate = progressUpdates[progressUpdates.length - 1]
      setDisplayedText(prev => prev + "\n" + latestUpdate)
    }
  }, [progressUpdates])
  
  return (
    <div className="w-full h-full flex flex-col bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] rounded-2xl border border-[#20808D]/30 overflow-hidden">
      {/* Progress bar */}
      <div className="relative h-1 bg-[#0A1A1A] overflow-hidden">
        <div 
          ref={progressRef}
          className="absolute top-0 left-0 h-full bg-gradient-to-r from-blue-500 to-blue-400 transition-all duration-500 ease-out"
          style={{ width: `${progress}%` }}
        >
          {/* Glowing effect */}
          <div className="absolute right-0 top-0 w-8 h-full bg-gradient-to-r from-transparent to-white/30 animate-pulse"></div>
        </div>
        
        {/* Progress complete indicator */}
        {progress === 100 && (
          <div className="absolute right-2 top-1/2 -translate-y-1/2 text-xs font-bold text-blue-400 animate-bounce">
            DONE!
          </div>
        )}
      </div>
      
      <div className="flex-1 p-6 overflow-hidden flex flex-col">
        {/* Header */}
        <div className="flex items-center gap-3 mb-4">
          <div className="text-3xl animate-pulse">{getHeaderIcon()}</div>
          <h3 className="text-xl font-bold text-[#20808D]">{getHeaderText()}</h3>
          {isThinking && progress < 100 && (
            <div className="flex gap-1">
              <div className="w-2 h-2 bg-[#20808D] rounded-full animate-bounce" style={{ animationDelay: "0ms" }}></div>
              <div className="w-2 h-2 bg-[#20808D] rounded-full animate-bounce" style={{ animationDelay: "150ms" }}></div>
              <div className="w-2 h-2 bg-[#20808D] rounded-full animate-bounce" style={{ animationDelay: "300ms" }}></div>
            </div>
          )}
        </div>
        
        {/* Progress content - with max height to prevent stretching */}
        <div className="flex-1 overflow-y-auto space-y-4 pr-2 custom-scrollbar max-h-[calc(100vh-300px)]">
          {/* Progress updates */}
          {progressUpdates && progressUpdates.length > 0 && (
            <div className="space-y-2">
              <h4 className="text-sm font-semibold text-[#20808D]/80">Progress</h4>
              <div className="space-y-1 max-h-40 overflow-y-auto">
                {progressUpdates.map((update, idx) => (
                  <div 
                    key={idx} 
                    className="text-sm text-[#FBFAF4]/80 pl-4 border-l-2 border-[#20808D]/30 animate-fadeIn"
                    style={{ animationDelay: `${idx * 100}ms` }}
                  >
                    {update}
                  </div>
                ))}
              </div>
            </div>
          )}
          
          {/* Sources (for research mode) - Enhanced display */}
          {mode === "deep-research" && sources && sources.length > 0 && (
            <div className="space-y-2">
              <h4 className="text-sm font-semibold text-[#20808D]/80">Sources Being Accessed ({sources.length})</h4>
              <div className="space-y-1 max-h-60 overflow-y-auto">
                {sources.map((source, idx) => (
                  <div 
                    key={idx} 
                    className="text-sm text-[#FBFAF4]/70 pl-4 border-l-2 border-blue-500/30 hover:border-blue-500/50 transition-colors"
                  >
                    <div className="flex items-center gap-2">
                      <div className="w-2 h-2 bg-blue-500 rounded-full animate-pulse"></div>
                      <div className="font-medium truncate">{source.domain || source.title || `Source ${idx + 1}`}</div>
                    </div>
                    {source.url && (
                      <div className="text-xs text-[#20808D]/60 mt-1 font-mono truncate">{source.url}</div>
                    )}
                    {source.description && (
                      <div className="text-xs text-[#FBFAF4]/50 mt-1 line-clamp-2">{source.description}</div>
                    )}
                  </div>
                ))}
              </div>
            </div>
          )}
          
          {/* Reasoning steps (for reasoning mode) */}
          {mode === "reasoning" && reasoning && reasoning.length > 0 && (
            <div className="space-y-2">
              <h4 className="text-sm font-semibold text-[#20808D]/80">Reasoning Steps</h4>
              <div className="space-y-1 max-h-60 overflow-y-auto">
                {reasoning.map((step, idx) => (
                  <div 
                    key={idx} 
                    className="text-sm text-[#FBFAF4]/70 pl-4 border-l-2 border-purple-500/30 hover:border-purple-500/50 transition-colors"
                  >
                    <div className="flex items-center gap-2 mb-1">
                      <div className="w-6 h-6 bg-purple-500/20 rounded-full flex items-center justify-center text-xs font-bold text-purple-400">
                        {step.step || idx + 1}
                      </div>
                      <div className="font-medium text-purple-400">Step {step.step || idx + 1}</div>
                    </div>
                    <div className="ml-8 text-xs">{step.content || step}</div>
                  </div>
                ))}
              </div>
            </div>
          )}
          
          {/* Citations */}
          {citations && citations.length > 0 && (
            <div className="space-y-2">
              <h4 className="text-sm font-semibold text-[#20808D]/80">References</h4>
              <div className="grid gap-2 max-h-40 overflow-y-auto">
                {citations.map((citation, idx) => (
                  <div 
                    key={idx} 
                    className="text-xs bg-[#20808D]/10 rounded p-2 border border-[#20808D]/20 hover:bg-[#20808D]/20 transition-colors"
                  >
                    <span className="font-bold text-[#20808D]">[{idx + 1}]</span> {citation.title || citation}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
        
        {/* Status bar */}
        <div className="mt-4 pt-4 border-t border-[#20808D]/20">
          <div className="flex justify-between items-center text-xs text-[#20808D]/60">
            <span>Mode: {mode === "reasoning" ? "Reasoning Pro" : "Deep Research"}</span>
            <span className="flex items-center gap-2">
              {progress < 100 ? (
                <>
                  <div className="w-4 h-4 border-2 border-[#20808D] border-t-transparent rounded-full animate-spin"></div>
                  Processing... {Math.round(progress)}%
                </>
              ) : (
                <>
                  <div className="w-4 h-4 bg-green-500 rounded-full flex items-center justify-center">
                    <span className="text-white text-[8px]">✓</span>
                  </div>
                  Complete
                </>
              )}
            </span>
          </div>
        </div>
      </div>
    </div>
  )
}

// Default state component for when no query is active
function DefaultVisualizationState() {
  return (
    <div className="w-full h-full relative overflow-hidden rounded-2xl border border-[#20808D]/30 bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717]">
      <div className="w-full h-full flex items-center justify-center p-8">
        <div className="text-center max-w-md">
          <div className="mb-6">
            <div className="w-20 h-20 mx-auto bg-gradient-to-br from-[#20808D]/20 to-[#20808D]/10 rounded-full flex items-center justify-center">
              <svg className="w-10 h-10 text-[#20808D]" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
              </svg>
            </div>
          </div>
          <h3 className="text-xl font-bold text-[#20808D] mb-3">Query Visualization</h3>
          <p className="text-[#FBFAF4]/60 text-sm leading-relaxed mb-4">
            Ask Monday a question to see beautiful 3D visualizations related to your query
          </p>
          <div className="flex items-center justify-center gap-2 text-xs text-[#20808D]/50">
            <div className="w-2 h-2 bg-[#20808D]/30 rounded-full"></div>
            <span>Ready for exploration</span>
          </div>
        </div>
      </div>
    </div>
  )
}

// Add custom scrollbar styles
const scrollbarStyles = `
  .custom-scrollbar::-webkit-scrollbar {
    width: 6px;
  }
  
  .custom-scrollbar::-webkit-scrollbar-track {
    background: rgba(32, 128, 141, 0.1);
    border-radius: 3px;
  }
  
  .custom-scrollbar::-webkit-scrollbar-thumb {
    background: rgba(32, 128, 141, 0.3);
    border-radius: 3px;
  }
  
  .custom-scrollbar::-webkit-scrollbar-thumb:hover {
    background: rgba(32, 128, 141, 0.5);
  }
  
  @keyframes fadeIn {
    from {
      opacity: 0;
      transform: translateY(10px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
  
  .animate-fadeIn {
    animation: fadeIn 0.5s ease-out forwards;
  }
`

// Inject styles
if (typeof document !== 'undefined') {
  const styleSheet = document.createElement('style')
  styleSheet.textContent = scrollbarStyles
  document.head.appendChild(styleSheet)
}

export function AdaptiveVisualizationPanel({
  mode,
  query = "",
  visualizationData,
  isThinking = false,
  progressUpdates = [],
  sources = [],
  reasoning = [],
  citations = []
}: AdaptiveVisualizationPanelProps) {
  // Determine if the process is complete
  const isComplete = !isThinking && (progressUpdates.length > 0 || sources.length > 0 || reasoning.length > 0)
  
  // If no query yet, show default state
  if (!query && !isThinking && progressUpdates.length === 0) {
    return <DefaultVisualizationState />
  }
  
  // For basic mode, show 3D visualization
  if (mode === "basic") {
    // Check if we have a specific 3D visualization for this query
    const visualization = generate3DVisualization(query)
    
    if (visualization.type !== "abstract" && visualization.objects.length > 0) {
      return (
        <div className="w-full h-full relative overflow-hidden rounded-2xl border border-[#20808D]/30">
          <div className="absolute top-6 left-6 z-10 bg-gradient-to-br from-[#20808D]/20 to-[#20808D]/10 backdrop-blur-xl rounded-2xl p-4 border border-[#20808D]/30 shadow-2xl max-w-xs">
            <div className="text-[#20808D] font-bold text-sm mb-2 flex items-center gap-2">
              <div className="w-2 h-2 bg-[#20808D] rounded-full animate-pulse"></div>
              Query Visualization
            </div>
            <div className="text-[#FBFAF4] text-xs mb-1">Query: {query}</div>
            <div className="text-[#20808D]/80 text-xs">Mode: Basic</div>
          </div>
          <Basic3DVisualization query={query} />
        </div>
      )
    } else {
      // Fall back to knowledge tree visualization
      return <BinaryTree3D visualizationData={visualizationData} title="Query Visualization" />
    }
  }
  
  // For reasoning and deep-research modes, show progressive text display
  return (
    <ProgressiveThinkingDisplay
      mode={mode}
      isThinking={isThinking}
      progressUpdates={progressUpdates}
      sources={sources}
      reasoning={reasoning}
      citations={citations}
      isComplete={isComplete}
    />
  )
} 