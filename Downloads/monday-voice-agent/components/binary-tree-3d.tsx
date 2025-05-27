"use client"

import { useRef, useState, useEffect } from "react"
import { Canvas, useFrame } from "@react-three/fiber"
import { Text, OrbitControls, Environment, Float } from "@react-three/drei"
import * as THREE from "three"
import type { JSX } from "react/jsx-runtime"

interface TreeNode {
  value: number
  left?: TreeNode
  right?: TreeNode
  position: [number, number, number]
  isHighlighted?: boolean
  isVisiting?: boolean
}

interface TreeNodeProps {
  node: TreeNode
  onNodeClick?: (value: number) => void
}

function TreeNodeComponent({ node, onNodeClick }: TreeNodeProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [hovered, setHovered] = useState(false)

  useFrame((state) => {
    if (meshRef.current) {
      if (node.isVisiting) {
        meshRef.current.scale.setScalar(1.3 + Math.sin(state.clock.elapsedTime * 6) * 0.15)
      } else {
        meshRef.current.scale.setScalar(hovered ? 1.15 : 1)
      }
    }
  })

  const nodeColor = node.isVisiting ? "#10B981" : node.isHighlighted ? "#20808D" : hovered ? "#3B82F6" : "#64748B"
  const emissiveColor = node.isVisiting ? "#065F46" : node.isHighlighted ? "#0F4C5C" : "#000000"

  return (
    <group position={node.position}>
      <Float speed={1.5} rotationIntensity={0.2} floatIntensity={0.3}>
        <mesh
          ref={meshRef}
          onClick={() => onNodeClick?.(node.value)}
          onPointerOver={() => setHovered(true)}
          onPointerOut={() => setHovered(false)}
        >
          <sphereGeometry args={[0.35, 32, 32]} />
          <meshStandardMaterial
            color={nodeColor}
            emissive={emissiveColor}
            emissiveIntensity={node.isVisiting ? 0.3 : 0.1}
            roughness={0.2}
            metalness={0.8}
          />
        </mesh>
        <Text
          position={[0, 0, 0.36]}
          fontSize={0.25}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
          font="/fonts/Inter-Bold.ttf"
        >
          {node.value}
        </Text>
      </Float>
    </group>
  )
}

function TreeConnections({ tree }: { tree: TreeNode }) {
  const connections: Array<{ start: [number, number, number]; end: [number, number, number] }> = []

  function addConnections(node: TreeNode) {
    if (node.left) {
      connections.push({ start: node.position, end: node.left.position })
      addConnections(node.left)
    }
    if (node.right) {
      connections.push({ start: node.position, end: node.right.position })
      addConnections(node.right)
    }
  }

  addConnections(tree)

  return (
    <>
      {connections.map((conn, index) => {
        const start = new THREE.Vector3(...conn.start)
        const end = new THREE.Vector3(...conn.end)
        const direction = end.clone().sub(start)
        const length = direction.length()
        const midpoint = start.clone().add(direction.clone().multiplyScalar(0.5))

        return (
          <mesh key={index} position={midpoint.toArray()}>
            <cylinderGeometry args={[0.03, 0.03, length]} />
            <meshStandardMaterial
              color="#20808D"
              opacity={0.8}
              transparent
              emissive="#0F4C5C"
              emissiveIntensity={0.2}
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

  const particleCount = 100
  const positions = new Float32Array(particleCount * 3)

  for (let i = 0; i < particleCount; i++) {
    positions[i * 3] = (Math.random() - 0.5) * 20
    positions[i * 3 + 1] = (Math.random() - 0.5) * 20
    positions[i * 3 + 2] = (Math.random() - 0.5) * 20
  }

  return (
    <points ref={particlesRef}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" array={positions} count={particleCount} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial color="#20808D" size={0.02} opacity={0.6} transparent />
    </points>
  )
}

export function BinaryTree3D() {
  const [tree, setTree] = useState<TreeNode>({
    value: 50,
    position: [0, 2, 0],
    left: {
      value: 30,
      position: [-2.5, 0.5, 0],
      left: {
        value: 20,
        position: [-3.5, -1, 0],
      },
      right: {
        value: 40,
        position: [-1.5, -1, 0],
      },
    },
    right: {
      value: 70,
      position: [2.5, 0.5, 0],
      left: {
        value: 60,
        position: [1.5, -1, 0],
      },
      right: {
        value: 80,
        position: [3.5, -1, 0],
      },
    },
  })

  const [currentOperation, setCurrentOperation] = useState<string>("Initializing...")
  const [visitingPath, setVisitingPath] = useState<number[]>([])

  useEffect(() => {
    if (!tree) return

    const searchValue = 60
    const path: number[] = []

    function findPath(node: TreeNode | undefined, target: number): boolean {
      if (!node) return false
      path.push(node.value)
      if (node.value === target) return true
      if (target < node.value) {
        return findPath(node.left, target)
      } else {
        return findPath(node.right, target)
      }
    }

    const animateSearch = async () => {
      setCurrentOperation(`Searching for ${searchValue}`)
      findPath(tree, searchValue)

      for (let i = 0; i < path.length; i++) {
        await new Promise((resolve) => setTimeout(resolve, 1800))
        setVisitingPath(path.slice(0, i + 1))
      }

      await new Promise((resolve) => setTimeout(resolve, 2500))
      setCurrentOperation("Target Found!")

      await new Promise((resolve) => setTimeout(resolve, 1500))
      setVisitingPath([])
      setCurrentOperation("Ready for next search")
    }

    const interval = setInterval(animateSearch, 10000)
    return () => clearInterval(interval)
  }, [tree])

  function updateTreeWithVisiting(node: TreeNode): TreeNode {
    return {
      ...node,
      isVisiting: visitingPath.includes(node.value),
      isHighlighted: visitingPath.includes(node.value),
      left: node.left ? updateTreeWithVisiting(node.left) : undefined,
      right: node.right ? updateTreeWithVisiting(node.right) : undefined,
    }
  }

  const animatedTree = updateTreeWithVisiting(tree)

  function renderTree(node: TreeNode): JSX.Element[] {
    const elements = [<TreeNodeComponent key={node.value} node={node} />]
    if (node.left) {
      elements.push(...renderTree(node.left))
    }
    if (node.right) {
      elements.push(...renderTree(node.right))
    }
    return elements
  }

  return (
    <div className="w-full h-full relative overflow-hidden">
      <div className="absolute top-6 left-6 z-10 bg-gradient-to-br from-[#20808D]/20 to-[#20808D]/10 backdrop-blur-xl rounded-2xl p-4 border border-[#20808D]/30 shadow-2xl">
        <div className="text-[#20808D] font-bold text-sm mb-2 flex items-center gap-2">
          <div className="w-2 h-2 bg-[#20808D] rounded-full animate-pulse"></div>
          Binary Search Tree
        </div>
        <div className="text-[#FBFAF4] text-xs mb-1">Status: {currentOperation}</div>
        {visitingPath.length > 0 && <div className="text-[#20808D] text-xs">Path: {visitingPath.join(" → ")}</div>}
      </div>

      <Canvas camera={{ position: [0, 3, 10], fov: 45 }}>
        <Environment preset="night" />
        <ambientLight intensity={0.3} />
        <pointLight position={[10, 10, 10]} intensity={1.5} color="#20808D" />
        <pointLight position={[-10, -10, -10]} intensity={0.8} color="#FBFAF4" />
        <spotLight position={[0, 15, 0]} intensity={1} color="#20808D" angle={0.3} penumbra={1} />

        <ParticleField />
        <TreeConnections tree={animatedTree} />
        {renderTree(animatedTree)}

        <OrbitControls
          enablePan={false}
          enableZoom={true}
          maxDistance={20}
          minDistance={6}
          autoRotate
          autoRotateSpeed={0.3}
          maxPolarAngle={Math.PI / 1.8}
          minPolarAngle={Math.PI / 4}
        />
      </Canvas>
    </div>
  )
}
