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
  thinking?: string
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

// Enhanced Coffee cup component with photorealistic details
function CoffeeCup({ brand, position }: { brand: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  const steamRef = useRef<THREE.Points>(null)
  const [hovered, setHovered] = useState(false)
  
  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = Math.sin(state.clock.elapsedTime * 0.3) * 0.2
      meshRef.current.scale.setScalar(hovered ? 1.1 : 1)
      meshRef.current.position.y = Math.sin(state.clock.elapsedTime * 2) * 0.02
    }
    
    // Animate steam particles with realistic physics
    if (steamRef.current) {
      const positions = steamRef.current.geometry.attributes.position.array as Float32Array
      const velocities = steamRef.current.geometry.attributes.velocity?.array as Float32Array || new Float32Array(positions.length)
      
      for (let i = 0; i < positions.length; i += 3) {
        // Add turbulence and wind effects
        velocities[i] += (Math.random() - 0.5) * 0.001 // X turbulence
        velocities[i + 1] += 0.015 + Math.random() * 0.005 // Y velocity (upward)
        velocities[i + 2] += (Math.random() - 0.5) * 0.001 // Z turbulence
        
        // Apply velocities
        positions[i] += velocities[i]
        positions[i + 1] += velocities[i + 1]
        positions[i + 2] += velocities[i + 2]
        
        // Reset particles that go too high
        if (positions[i + 1] > 3) {
          positions[i] = (Math.random() - 0.5) * 0.3
          positions[i + 1] = 0.8
          positions[i + 2] = (Math.random() - 0.5) * 0.3
          velocities[i] = 0
          velocities[i + 1] = 0
          velocities[i + 2] = 0
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
  
  // Create enhanced steam particles
  const steamParticles = useMemo(() => {
    const particles = new Float32Array(60) // Doubled particles for better effect
    for (let i = 0; i < particles.length; i += 3) {
      particles[i] = (Math.random() - 0.5) * 0.3
      particles[i + 1] = 0.8 + Math.random() * 0.5
      particles[i + 2] = (Math.random() - 0.5) * 0.3
    }
    return particles
  }, [])
  
  return (
    <Float speed={1.5} rotationIntensity={0.1} floatIntensity={0.2}>
      <group 
        ref={meshRef} 
        position={position}
        onPointerOver={() => setHovered(true)}
        onPointerOut={() => setHovered(false)}
      >
        {/* Enhanced cup body with ceramic material */}
        <mesh position={[0, 0, 0]} castShadow receiveShadow>
          <cylinderGeometry args={[0.6, 0.5, 1.5, 64, 4, false]} />
          <meshPhysicalMaterial 
            color={getBrandColor()} 
            metalness={0.1} 
            roughness={0.2}
            clearcoat={1}
            clearcoatRoughness={0.1}
            reflectivity={0.8}
            envMapIntensity={1.5}
            normalScale={new THREE.Vector2(0.5, 0.5)}
          />
        </mesh>
        
        {/* Ceramic rim with gold accent */}
        <mesh position={[0, 0.75, 0]} castShadow>
          <torusGeometry args={[0.6, 0.06, 16, 64]} />
          <meshPhysicalMaterial 
            color="#D4AF37" 
            metalness={0.9} 
            roughness={0.1}
            clearcoat={1}
            envMapIntensity={2}
          />
        </mesh>
        
        {/* Enhanced handle with detailed geometry */}
        <mesh position={[0.7, 0, 0]} rotation={[0, 0, Math.PI / 2]} castShadow>
          <torusGeometry args={[0.3, 0.08, 16, 64, Math.PI]} />
          <meshPhysicalMaterial 
            color={getBrandColor()} 
            metalness={0.1} 
            roughness={0.2}
            clearcoat={1}
            clearcoatRoughness={0.1}
          />
        </mesh>
        
        {/* Realistic coffee surface with foam */}
        <mesh position={[0, 0.6, 0]} receiveShadow>
          <cylinderGeometry args={[0.55, 0.55, 0.1, 64]} />
          <meshPhysicalMaterial 
            color="#2B1810" 
            metalness={0} 
            roughness={0.9}
            clearcoat={0.3}
            normalScale={new THREE.Vector2(2, 2)}
          />
        </mesh>
        
        {/* Foam layer with micro-foam texture */}
        <mesh position={[0, 0.65, 0]} receiveShadow>
          <cylinderGeometry args={[0.54, 0.54, 0.05, 64]} />
          <meshPhysicalMaterial 
            color="#F5DEB3" 
            metalness={0} 
            roughness={1}
            transmission={0.1}
            opacity={0.95}
            transparent
          />
        </mesh>
        
        {/* Enhanced steam with volume rendering effect */}
        <points ref={steamRef}>
          <bufferGeometry>
            <bufferAttribute
              attach="attributes-position"
              count={20}
              array={steamParticles}
              itemSize={3}
              args={[steamParticles, 3]}
            />
          </bufferGeometry>
          <pointsMaterial
            size={0.08}
            color="#FFFFFF"
            transparent
            opacity={0.6}
            sizeAttenuation
            blending={THREE.AdditiveBlending}
            alphaMap={new THREE.TextureLoader().load('data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iNjQiIGhlaWdodD0iNjQiIHZpZXdCb3g9IjAgMCA2NCA2NCIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPGNpcmNsZSBjeD0iMzIiIGN5PSIzMiIgcj0iMzIiIGZpbGw9InVybCgjZ3JhZGllbnQwX3JhZGlhbF8xXzIpIiBmaWxsLW9wYWNpdHk9IjAuOCIvPgo8ZGVmcz4KPHJhZGlhbEdyYWRpZW50IGlkPSJncmFkaWVudDBfcmFkaWFsXzFfMiIgY3g9IjAiIGN5PSIwIiByPSIxIiBncmFkaWVudFVuaXRzPSJ1c2VyU3BhY2VPblVzZSIgZ3JhZGllbnRUcmFuc2Zvcm09InRyYW5zbGF0ZSgzMiAzMikgcm90YXRlKDkwKSBzY2FsZSgzMikiPgo8c3RvcCBzdG9wLWNvbG9yPSJ3aGl0ZSIvPgo8c3RvcCBvZmZzZXQ9IjEiIHN0b3AtY29sb3I9IndoaXRlIiBzdG9wLW9wYWNpdHk9IjAiLz4KPC9yYWRpYWxHcmFkaWVudD4KPC9kZWZzPgo8L3N2Zz4K')}
          />
        </points>
        
        {/* Brand logo with 3D embossed effect */}
        <group position={[0, 0, 0.61]}>
          <mesh>
            <planeGeometry args={[0.8, 0.3]} />
            <meshPhysicalMaterial 
              color="#FFFFFF" 
              opacity={0.95} 
              transparent
              clearcoat={0.5}
              normalScale={new THREE.Vector2(1, 1)}
            />
          </mesh>
          <Text
            position={[0, 0, 0.01]}
            fontSize={0.15}
            color={getBrandColor()}
            anchorX="center"
            anchorY="middle"
            font="/fonts/Roboto-Bold.woff"
            outlineWidth={0.005}
            outlineColor="#000000"
          >
            {brand}
          </Text>
        </group>
        
        {/* Ceramic base with shadow */}
        <mesh position={[0, -0.8, 0]} rotation={[-Math.PI / 2, 0, 0]} receiveShadow>
          <circleGeometry args={[0.7, 64]} />
          <meshPhysicalMaterial 
            color="#FFFFFF" 
            metalness={0} 
            roughness={0.3}
            opacity={0.8}
            transparent
          />
        </mesh>
        
        {/* Ambient lighting for the cup */}
        <pointLight
          position={[0, 1, 1]}
          intensity={0.3}
          color="#FFFACD"
          distance={3}
          decay={2}
        />
        
        {/* Bottom label with embossed text */}
        <Text
          position={[0, -1.2, 0]}
          fontSize={0.2}
          color="#D4AF37"
          anchorX="center"
          anchorY="middle"
          font="/fonts/Roboto-Bold.woff"
          outlineWidth={0.01}
          outlineColor="#8B4513"
        >
          Premium {brand}
        </Text>
      </group>
    </Float>
  )
}

// Enhanced Tree component with photorealistic bark, dynamic leaves, and seasonal effects
function Tree({ variant, position }: { variant: string, position: [number, number, number] }) {
  const meshRef = useRef<THREE.Group>(null)
  const leavesRef = useRef<THREE.InstancedMesh>(null)
  const windRef = useRef({ time: 0, strength: 0.5 })
  
  useFrame((state) => {
    windRef.current.time = state.clock.elapsedTime
    windRef.current.strength = 0.3 + Math.sin(state.clock.elapsedTime * 0.5) * 0.2
    
    if (meshRef.current) {
      // Gentle swaying motion based on wind
      const windX = Math.sin(state.clock.elapsedTime * 0.8) * windRef.current.strength * 0.02
      const windZ = Math.cos(state.clock.elapsedTime * 0.6) * windRef.current.strength * 0.015
      meshRef.current.rotation.x = windX
      meshRef.current.rotation.z = windZ
      
      // Animate leaves with wind physics
      if (leavesRef.current && variant !== "pine") {
        const time = state.clock.elapsedTime
        const leafCount = 40 // Increased leaf count
        
        for (let i = 0; i < leafCount; i++) {
          const matrix = new THREE.Matrix4()
          
          // Base position on branches
          const angle = (i / leafCount) * Math.PI * 4
          const radius = 0.8 + Math.sin(i * 0.5) * 0.5
          const height = 1.5 + Math.sin(i * 0.3) * 1.2
          
          // Wind effect on individual leaves
          const windEffect = windRef.current.strength * Math.sin(time * 2 + i * 0.1)
          
          const x = Math.cos(angle) * radius + windEffect * 0.1
          const y = height + Math.sin(time * 1.5 + i * 0.2) * 0.05
          const z = Math.sin(angle) * radius + windEffect * 0.08
          
          // Varying leaf sizes based on position
          const scale = 0.15 + Math.sin(i * 0.7) * 0.05
          
          // Leaf rotation with wind
          const rotX = windEffect * 0.3
          const rotY = angle + windEffect * 0.2
          const rotZ = Math.sin(time + i * 0.1) * 0.1
          
          matrix.makeRotationFromEuler(new THREE.Euler(rotX, rotY, rotZ))
          matrix.setPosition(x, y, z)
          matrix.scale(new THREE.Vector3(scale, scale, scale))
          
          leavesRef.current.setMatrixAt(i, matrix)
        }
        leavesRef.current.instanceMatrix.needsUpdate = true
      }
    }
  })
  
  const getTreeColors = () => {
    switch (variant) {
      case "oak": return { 
        leaves: "#2D5016", 
        bark: "#8B4513", 
        accent: "#228B22",
        season: "spring"
      }
      case "pine": return { 
        leaves: "#0F4B0F", 
        bark: "#654321", 
        accent: "#00563F",
        season: "evergreen"
      }
      case "maple": return { 
        leaves: "#FF8C00", 
        bark: "#A0522D", 
        accent: "#FF6347",
        season: "autumn"
      }
      default: return { 
        leaves: "#228B22", 
        bark: "#8B4513", 
        accent: "#32CD32",
        season: "summer"
      }
    }
  }
  
  const colors = getTreeColors()
  
  return (
    <group ref={meshRef} position={position}>
      {/* Enhanced trunk with realistic bark texture */}
      <mesh position={[0, -1, 0]} castShadow receiveShadow>
        <cylinderGeometry args={[0.35, 0.45, 2.2, 24, 8]} />
        <meshPhysicalMaterial 
          color={colors.bark}
          roughness={0.95}
          metalness={0}
          normalScale={new THREE.Vector2(3, 3)}
          bumpScale={0.02}
          envMapIntensity={0.3}
        />
      </mesh>
      
      {/* Detailed root system */}
      {[...Array(6)].map((_, i) => (
        <mesh 
          key={`root-${i}`}
          position={[
            Math.cos(i * 1.047) * 0.6,
            -1.9,
            Math.sin(i * 1.047) * 0.6
          ]}
          rotation={[Math.PI / 6, i * 1.047, 0]}
          castShadow
        >
          <cylinderGeometry args={[0.08, 0.15, 0.8, 12]} />
          <meshPhysicalMaterial 
            color={colors.bark} 
            roughness={0.9}
            normalScale={new THREE.Vector2(2, 2)}
          />
        </mesh>
      ))}
      
      {/* Complex branch system */}
      {[...Array(8)].map((_, i) => {
        const angle = (i / 8) * Math.PI * 2
        const branchLength = 1.2 + Math.sin(i * 0.5) * 0.3
        return (
          <group key={`branch-${i}`}>
            {/* Main branch */}
            <mesh 
              position={[
                Math.cos(angle) * 0.3,
                0.3 + i * 0.25,
                Math.sin(angle) * 0.3
              ]}
              rotation={[0, angle, Math.PI / 4 + Math.sin(i) * 0.2]}
              castShadow
            >
              <cylinderGeometry args={[0.03, 0.08, branchLength, 12]} />
              <meshPhysicalMaterial 
                color={colors.bark} 
                roughness={0.8}
                metalness={0}
              />
            </mesh>
            
            {/* Sub-branches */}
            {[...Array(3)].map((_, j) => (
              <mesh 
                key={`subbranch-${i}-${j}`}
                position={[
                  Math.cos(angle) * (0.5 + j * 0.2),
                  0.8 + i * 0.25 + j * 0.1,
                  Math.sin(angle) * (0.5 + j * 0.2)
                ]}
                rotation={[
                  Math.sin(j) * 0.3, 
                  angle + j * 0.5, 
                  Math.PI / 3 + j * 0.2
                ]}
                castShadow
              >
                <cylinderGeometry args={[0.01, 0.04, 0.6, 8]} />
                <meshPhysicalMaterial 
                  color={colors.bark} 
                  roughness={0.75}
                />
              </mesh>
            ))}
          </group>
        )
      })}
      
      {/* Pine cone shape for evergreen or enhanced canopy for deciduous */}
      {variant === "pine" ? (
        <>
          {/* Layered pine needle structure */}
          <mesh position={[0, 1, 0]} castShadow receiveShadow>
            <coneGeometry args={[1.6, 3.2, 16]} />
            <meshPhysicalMaterial 
              color={colors.leaves}
              roughness={0.8}
              metalness={0.1}
              normalScale={new THREE.Vector2(1.5, 1.5)}
              envMapIntensity={0.5}
            />
          </mesh>
          <mesh position={[0, 2.4, 0]} castShadow receiveShadow>
            <coneGeometry args={[1.3, 2.6, 16]} />
            <meshPhysicalMaterial 
              color={colors.accent}
              roughness={0.8}
              metalness={0.1}
            />
          </mesh>
          <mesh position={[0, 3.5, 0]} castShadow receiveShadow>
            <coneGeometry args={[1, 2, 16]} />
            <meshPhysicalMaterial 
              color={colors.leaves}
              roughness={0.75}
            />
          </mesh>
        </>
      ) : (
        <>
          {/* Volumetric leafy canopy with seasonal colors */}
          <mesh position={[0, 1.8, 0]} castShadow receiveShadow>
            <sphereGeometry args={[2.2, 32, 32]} />
            <meshPhysicalMaterial 
              color={colors.leaves}
              roughness={0.9}
              metalness={0}
              transmission={0.1}
              opacity={0.85}
              transparent
              normalScale={new THREE.Vector2(0.8, 0.8)}
            />
          </mesh>
          
          {/* Secondary leaf layer with accent color */}
          <mesh position={[0.2, 2.1, -0.1]} castShadow receiveShadow>
            <sphereGeometry args={[1.8, 24, 24]} />
            <meshPhysicalMaterial 
              color={colors.accent}
              roughness={0.85}
              opacity={0.6}
              transparent
            />
          </mesh>
          
          {/* Individual animated leaves */}
          <instancedMesh ref={leavesRef} args={[undefined, undefined, 40]} castShadow>
            <planeGeometry args={[0.25, 0.35]} />
            <meshPhysicalMaterial 
              color={colors.leaves}
              side={THREE.DoubleSide}
              opacity={0.9}
              transparent
              roughness={0.7}
              normalScale={new THREE.Vector2(0.5, 0.5)}
            />
          </instancedMesh>
          
          {/* Falling leaves effect for autumn trees */}
          {colors.season === "autumn" && (
            <points>
              <bufferGeometry>
                <bufferAttribute
                  attach="attributes-position"
                  count={15}
                  array={new Float32Array(Array.from({length: 45}, () => (Math.random() - 0.5) * 6))}
                  itemSize={3}
                  args={[new Float32Array(Array.from({length: 45}, () => (Math.random() - 0.5) * 6)), 3]}
                />
              </bufferGeometry>
              <pointsMaterial
                size={0.15}
                color={colors.accent}
                transparent
                opacity={0.7}
                sizeAttenuation
              />
            </points>
          )}
        </>
      )}
      
      {/* Environmental ground vegetation */}
      <mesh position={[0, -2.1, 0]} rotation={[-Math.PI / 2, 0, 0]} receiveShadow>
        <ringGeometry args={[0.5, 2.5, 32]} />
        <meshPhysicalMaterial 
          color="#2F4F2F" 
          opacity={0.7} 
          transparent
          roughness={0.9}
        />
      </mesh>
      
      {/* Ambient tree lighting */}
      <pointLight
        position={[0, 2, 1]}
        intensity={0.2}
        color={colors.season === "autumn" ? "#FFB347" : "#90EE90"}
        distance={6}
        decay={2}
      />
      
      {/* Enhanced species label */}
      <Text
        position={[0, -2.8, 0]}
        fontSize={0.25}
        color={colors.accent}
        anchorX="center"
        anchorY="middle"
        font="/fonts/Roboto-Bold.woff"
        outlineWidth={0.01}
        outlineColor="#000000"
      >
        {variant.charAt(0).toUpperCase() + variant.slice(1)} Tree
      </Text>
      
      {/* Seasonal indicator */}
      <Text
        position={[0, -3.2, 0]}
        fontSize={0.15}
        color="#FBFAF4"
        anchorX="center"
        anchorY="middle"
      >
        {colors.season.charAt(0).toUpperCase() + colors.season.slice(1)} Season
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

// Basic 3D visualization wrapper
function Basic3DVisualization({ query }: { query: string }) {
  const [error, setError] = useState(false)
  
  // Handle WebGL context lost
  useEffect(() => {
    const handleContextLost = (e: Event) => {
      e.preventDefault()
      console.error('WebGL context lost')
      setError(true)
    }
    
    const handleContextRestored = () => {
      console.log('WebGL context restored')
      setError(false)
    }
    
    window.addEventListener('webglcontextlost', handleContextLost)
    window.addEventListener('webglcontextrestored', handleContextRestored)
    
    return () => {
      window.removeEventListener('webglcontextlost', handleContextLost)
      window.removeEventListener('webglcontextrestored', handleContextRestored)
    }
  }, [])
  
  const visualization = generate3DVisualization(query)
  
  if (error) {
    return (
      <div className="w-full h-full flex items-center justify-center bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717]">
        <div className="text-center">
          <p className="text-[#20808D] mb-4">3D visualization temporarily unavailable</p>
          <button 
            onClick={() => setError(false)}
            className="px-4 py-2 bg-[#20808D]/20 hover:bg-[#20808D]/30 rounded-lg text-[#FBFAF4] text-sm"
          >
            Retry
          </button>
        </div>
      </div>
    )
  }
  
  return (
    <div className="w-full h-full relative">
      <Canvas
        camera={{ position: [0, 2, 8], fov: 50 }}
        shadows
        dpr={[1, 1.5]} // Limit pixel ratio to prevent performance issues
        gl={{ 
          antialias: true,
          alpha: true,
          powerPreference: 'high-performance',
          failIfMajorPerformanceCaveat: false
        }}
        onCreated={({ gl }) => {
          gl.shadowMap.enabled = true
          gl.shadowMap.type = THREE.PCFSoftShadowMap
        }}
        onError={(error) => {
          console.error('Canvas error:', error)
          setError(true)
        }}
      >
        {/* Simplified lighting */}
        <ambientLight intensity={0.5} />
        <directionalLight
          position={[5, 10, 5]}
          intensity={1}
          castShadow
          shadow-mapSize={[1024, 1024]}
          shadow-camera-near={0.1}
          shadow-camera-far={50}
          shadow-camera-left={-10}
          shadow-camera-right={10}
          shadow-camera-top={10}
          shadow-camera-bottom={-10}
        />
        
        {/* Environment with lower resolution */}
        <Environment preset="city" resolution={256} />
        
        {/* Fog for depth */}
        <fog attach="fog" args={['#091717', 10, 30]} />
        
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
    </div>
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
  isComplete = false,
  thinking
}: {
  mode: "reasoning" | "deep-research"
  isThinking?: boolean
  progressUpdates?: string[]
  sources?: any[]
  reasoning?: any[]
  citations?: any[]
  isComplete?: boolean
  thinking?: string
}) {
  const [displayedText, setDisplayedText] = useState("")
  const [currentLine, setCurrentLine] = useState(0)
  const [progress, setProgress] = useState(0)
  const progressRef = useRef<HTMLDivElement>(null)
  const scrollContainerRef = useRef<HTMLDivElement>(null)
  
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
  
  // Auto-scroll to bottom when new content is added
  useEffect(() => {
    if (scrollContainerRef.current) {
      scrollContainerRef.current.scrollTop = scrollContainerRef.current.scrollHeight
    }
  }, [progressUpdates, sources, reasoning, citations])
  
  // Update progress based on updates
  useEffect(() => {
    if (isThinking && (progressUpdates?.length || reasoning?.length || sources?.length)) {
      const totalUpdates = (progressUpdates?.length || 0) + (reasoning?.length || 0) + (sources?.length || 0)
      const estimatedTotal = mode === "reasoning" ? 8 : 12 // More steps for research
      const currentProgress = Math.min((totalUpdates / estimatedTotal) * 100, 90)
      setProgress(currentProgress)
    }
    
    if (isComplete) {
      setProgress(100)
    }
  }, [progressUpdates, reasoning, sources, isThinking, isComplete, mode])
  
  // Show content immediately when it arrives
  useEffect(() => {
    if (progressUpdates && progressUpdates.length > 0) {
      setDisplayedText(progressUpdates.join('\n'))
    }
  }, [progressUpdates])
  
  // Check if we have any streaming content to display
  const hasStreamingContent = (progressUpdates && progressUpdates.length > 0) || 
                              (reasoning && reasoning.length > 0) || 
                              (sources && sources.length > 0) ||
                              (citations && citations.length > 0)
  
  return (
    <div className="w-full h-full flex flex-col bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] rounded-2xl border border-[#20808D]/30 overflow-hidden">
      {/* Enhanced progress bar with glow effect */}
      <div className="relative h-2 bg-[#0A1A1A] overflow-hidden flex-shrink-0">
        <div 
          ref={progressRef}
          className="absolute top-0 left-0 h-full bg-gradient-to-r from-blue-500 via-purple-500 to-blue-400 transition-all duration-500 ease-out shadow-lg"
          style={{ width: `${progress}%` }}
        >
          {/* Animated scanning effect */}
          <div className="absolute right-0 top-0 w-12 h-full bg-gradient-to-r from-transparent via-white/40 to-transparent animate-pulse"></div>
        </div>
        
        {/* Progress complete indicator */}
        {progress === 100 && (
          <div className="absolute right-3 top-1/2 -translate-y-1/2 text-xs font-bold text-green-400 animate-bounce">
            ✓ COMPLETE
          </div>
        )}
      </div>
      
      {/* Main content container with flex */}
      <div className="flex-1 flex flex-col min-h-0 p-6">
        {/* Enhanced header with status */}
        <div className="flex items-center gap-3 mb-6 flex-shrink-0">
          <div className="text-4xl animate-pulse">{getHeaderIcon()}</div>
          <div className="flex-1">
            <h3 className="text-xl font-bold text-[#20808D] mb-1">{getHeaderText()}</h3>
            <div className="text-sm text-[#20808D]/70">
              {hasStreamingContent ? `Processing... ${Math.round(progress)}%` : "Waiting for content..."}
            </div>
          </div>
          {isThinking && progress < 100 && (
            <div className="flex gap-1">
              <div className="w-3 h-3 bg-[#20808D] rounded-full animate-bounce" style={{ animationDelay: "0ms" }}></div>
              <div className="w-3 h-3 bg-[#20808D] rounded-full animate-bounce" style={{ animationDelay: "200ms" }}></div>
              <div className="w-3 h-3 bg-[#20808D] rounded-full animate-bounce" style={{ animationDelay: "400ms" }}></div>
            </div>
          )}
        </div>
        
        {/* Content display area */}
        {hasStreamingContent ? (
          <div 
            ref={scrollContainerRef}
            className="flex-1 overflow-y-auto space-y-4 pr-3 custom-scrollbar min-h-0"
            style={{ maxHeight: 'calc(100% - 140px)' }}
          >
            {/* Complete thinking process (when available) */}
            {thinking && thinking.trim() !== "" && (
              <div className="space-y-3">
                <h4 className="text-sm font-semibold text-[#20808D] flex items-center gap-2">
                  <span className="w-2 h-2 bg-indigo-500 rounded-full animate-pulse"></span>
                  {mode === "reasoning" ? "Complete Reasoning Process" : "Complete Research Process"}
                </h4>
                <div className="text-sm text-[#FBFAF4]/90 pl-4 border-l-2 border-indigo-500/40 bg-gradient-to-r from-indigo-500/10 to-transparent p-4 rounded-r-lg">
                  <pre className="whitespace-pre-wrap font-mono text-xs leading-relaxed">
                    {thinking}
                  </pre>
                </div>
              </div>
            )}
            
            {/* Progress updates - Enhanced display */}
            {progressUpdates && progressUpdates.length > 0 && (
              <div className="space-y-3">
                <h4 className="text-sm font-semibold text-[#20808D] flex items-center gap-2">
                  <span className="w-2 h-2 bg-blue-500 rounded-full animate-pulse"></span>
                  Live Progress
                </h4>
                {progressUpdates.slice(-8).map((update, idx) => (
                  <div 
                    key={idx} 
                    className="text-sm text-[#FBFAF4]/90 pl-4 border-l-2 border-blue-500/40 bg-gradient-to-r from-blue-500/10 to-transparent p-3 rounded-r-lg animate-fadeIn"
                    style={{ animationDelay: `${idx * 150}ms` }}
                  >
                    <div className="flex items-start gap-2">
                      <span className="text-blue-400 font-mono text-xs">#{idx + 1}</span>
                      <span className="flex-1">{update}</span>
                    </div>
                  </div>
                ))}
                {progressUpdates.length > 8 && (
                  <div className="text-xs text-[#20808D]/60 pl-4 italic">
                    + {progressUpdates.length - 8} more updates...
                  </div>
                )}
              </div>
            )}
            
            {/* Sources (for research mode) - Enhanced display */}
            {mode === "deep-research" && sources && sources.length > 0 && (
              <div className="space-y-3">
                <h4 className="text-sm font-semibold text-[#20808D] flex items-center gap-2">
                  <span className="w-2 h-2 bg-green-500 rounded-full animate-pulse"></span>
                  Sources Found ({sources.length})
                </h4>
                {sources.slice(-6).map((source, idx) => (
                  <div 
                    key={idx} 
                    className="text-sm text-[#FBFAF4]/80 pl-4 border-l-2 border-green-500/40 bg-gradient-to-r from-green-500/10 to-transparent p-3 rounded-r-lg hover:from-green-500/20 transition-colors"
                  >
                    <div className="flex items-center gap-3 mb-2">
                      <div className="w-3 h-3 bg-green-500 rounded-full animate-pulse"></div>
                      <div className="font-medium text-green-400 truncate flex-1">
                        {source.title || `Source ${idx + 1}`}
                      </div>
                    </div>
                    {source.url && (
                      <div className="text-xs text-[#20808D]/70 mb-1 font-mono truncate pl-6">
                        🔗 {source.url}
                      </div>
                    )}
                    {source.description && (
                      <div className="text-xs text-[#FBFAF4]/60 pl-6 line-clamp-2">
                        {source.description}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            )}
            
            {/* Reasoning steps (for reasoning mode) - Enhanced display */}
            {mode === "reasoning" && reasoning && reasoning.length > 0 && (
              <div className="space-y-3">
                <h4 className="text-sm font-semibold text-[#20808D] flex items-center gap-2">
                  <span className="w-2 h-2 bg-purple-500 rounded-full animate-pulse"></span>
                  Reasoning Steps
                </h4>
                {reasoning.slice(-6).map((step, idx) => (
                  <div 
                    key={idx} 
                    className="text-sm text-[#FBFAF4]/80 pl-4 border-l-2 border-purple-500/40 bg-gradient-to-r from-purple-500/10 to-transparent p-3 rounded-r-lg hover:from-purple-500/20 transition-colors"
                  >
                    <div className="flex items-center gap-3 mb-2">
                      <div className="w-6 h-6 bg-purple-500/30 rounded-full flex items-center justify-center text-xs font-bold text-purple-300 border border-purple-500/50">
                        {idx + 1}
                      </div>
                      <div className="font-medium text-purple-400">Analysis Step {idx + 1}</div>
                    </div>
                    <div className="ml-9 text-[#FBFAF4]/90 line-clamp-4">
                      {step.content || step}
                    </div>
                  </div>
                ))}
              </div>
            )}
            
            {/* Citations - Enhanced display */}
            {citations && citations.length > 0 && (
              <div className="space-y-3">
                <h4 className="text-sm font-semibold text-[#20808D] flex items-center gap-2">
                  <span className="w-2 h-2 bg-yellow-500 rounded-full animate-pulse"></span>
                  References ({citations.length})
                </h4>
                <div className="grid gap-2">
                  {citations.slice(-6).map((citation, idx) => (
                    <div 
                      key={idx} 
                      className="text-xs bg-gradient-to-r from-yellow-500/20 to-yellow-500/5 rounded-lg p-3 border border-yellow-500/30 hover:from-yellow-500/30 transition-colors"
                    >
                      <div className="flex items-start gap-2">
                        <span className="font-bold text-yellow-400 text-sm">[{idx + 1}]</span> 
                        <span className="flex-1 line-clamp-2 text-[#FBFAF4]/80">
                          {citation.title || citation}
                        </span>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        ) : (
          // No content message with enhanced styling
          <div className="flex-1 flex items-center justify-center">
            <div className="text-center">
              <div className="w-16 h-16 mx-auto mb-4 bg-gradient-to-br from-[#20808D]/30 to-[#20808D]/10 rounded-full flex items-center justify-center">
                <div className="text-2xl">⚡</div>
              </div>
              <h3 className="text-lg font-semibold text-[#20808D] mb-2">Processing Started</h3>
              <p className="text-sm text-[#20808D]/60 max-w-xs mx-auto">
                {mode === "reasoning" ? 
                  "Monday is analyzing your query with deep reasoning..." : 
                  "Monday is conducting comprehensive research..."}
              </p>
              <div className="mt-4 flex justify-center">
                <div className="w-8 h-8 border-2 border-[#20808D] border-t-transparent rounded-full animate-spin"></div>
              </div>
            </div>
          </div>
        )}
        
        {/* Enhanced status bar */}
        <div className="mt-4 pt-4 border-t border-[#20808D]/20 flex-shrink-0">
          <div className="flex justify-between items-center text-xs">
            <div className="flex items-center gap-3">
              <span className="text-[#20808D]/70">Mode:</span>
              <span className="px-2 py-1 bg-gradient-to-r from-[#20808D]/20 to-[#20808D]/10 rounded-full text-[#20808D] font-medium">
                {mode === "reasoning" ? "Reasoning Pro" : "Deep Research"}
              </span>
            </div>
            <div className="flex items-center gap-2 text-[#20808D]/60">
              {progress < 100 ? (
                <>
                  <div className="w-3 h-3 border border-[#20808D] border-t-transparent rounded-full animate-spin"></div>
                  <span>Processing {Math.round(progress)}%</span>
                </>
              ) : (
                <>
                  <div className="w-3 h-3 bg-green-500 rounded-full flex items-center justify-center">
                    <span className="text-white text-[6px]">✓</span>
                  </div>
                  <span className="text-green-400 font-medium">Complete</span>
                </>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

// Placeholder component for visualization
function VisualizationPlaceholder() {
  return (
    <div className="w-full h-full flex items-center justify-center bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] rounded-2xl border border-[#20808D]/30">
      <div className="text-center">
        <div className="mb-4">
          <div className="w-20 h-20 mx-auto bg-[#20808D]/20 rounded-full flex items-center justify-center">
            <svg className="w-10 h-10 text-[#20808D]" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" />
            </svg>
          </div>
        </div>
        <h3 className="text-lg font-semibold text-[#20808D] mb-2">Query Monday to visualize concepts</h3>
        <p className="text-sm text-[#20808D]/60 max-w-xs mx-auto">
          Ask about any topic and I'll create interactive 3D visualizations to help you understand better
        </p>
      </div>
    </div>
  )
}

export function AdaptiveVisualizationPanel({
  mode,
  query = "",
  visualizationData,
  isThinking = false,
  progressUpdates = [],
  sources = [],
  reasoning = [],
  citations = [],
  thinking
}: AdaptiveVisualizationPanelProps) {
  const [use3D, setUse3D] = useState(true)
  const [renderError, setRenderError] = useState(false)
  
  // Error boundary for 3D rendering
  useEffect(() => {
    const handleError = (e: ErrorEvent) => {
      if (e.message.includes('WebGL') || e.message.includes('THREE')) {
        console.error('3D rendering error:', e)
        setRenderError(true)
        setUse3D(false)
      }
    }
    
    window.addEventListener('error', handleError)
    return () => window.removeEventListener('error', handleError)
  }, [])
  
  // Determine if the process is complete
  const isComplete = !isThinking && (progressUpdates.length > 0 || sources.length > 0 || reasoning.length > 0)
  
  // For basic mode, show 3D visualization or placeholder
  if (mode === "basic") {
    // If no query yet, show placeholder
    if (!query || query.trim() === "") {
      return <VisualizationPlaceholder />
    }
    
    // If 3D rendering has errors, show simple visualization
    if (renderError) {
      return (
        <div className="w-full h-full bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] rounded-2xl border border-[#20808D]/30 overflow-hidden">
          <div className="h-full flex flex-col items-center justify-center p-8">
            <div className="text-center mb-8">
              <div className="w-24 h-24 mx-auto mb-4 bg-gradient-to-br from-[#20808D]/20 to-[#20808D]/10 rounded-full flex items-center justify-center">
                <div className="text-4xl">
                  {query.toLowerCase().includes('coffee') ? '☕' : 
                   query.toLowerCase().includes('tree') ? '🌳' :
                   query.toLowerCase().includes('space') ? '🌌' :
                   query.toLowerCase().includes('food') ? '🍕' :
                   query.toLowerCase().includes('music') ? '🎵' :
                   query.toLowerCase().includes('sport') ? '⚽' :
                   '💡'}
                </div>
              </div>
              <h3 className="text-xl font-bold text-[#20808D] mb-2">Visualization Active</h3>
              <p className="text-[#FBFAF4]/60 text-sm">Query: {query}</p>
              <p className="text-[#20808D]/60 text-xs mt-2">Mode: Basic Search</p>
            </div>
            
            <button
              onClick={() => {
                setRenderError(false)
                setUse3D(true)
              }}
              className="px-4 py-2 bg-[#20808D]/20 hover:bg-[#20808D]/30 rounded-lg text-[#FBFAF4] text-sm transition-colors"
            >
              Retry 3D Visualization
            </button>
          </div>
        </div>
      )
    }
    
    // Check if we have a specific 3D visualization for this query
    const visualization = generate3DVisualization(query)
    
    // Always try to show 3D visualization if use3D is true
    if (use3D) {
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
      thinking={thinking}
    />
  )
} 