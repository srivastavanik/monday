import { useRef } from 'react'
import { useFrame } from '@react-three/fiber'
import { Text, Sphere } from '@react-three/drei'

const MondayScene = () => {
  const asteriskRef = useRef<any>()

  // Animate Monday's presence indicator (asterisk)
  useFrame((state) => {
    if (asteriskRef.current) {
      asteriskRef.current.rotation.z = Math.sin(state.clock.elapsedTime) * 0.1
      asteriskRef.current.scale.setScalar(1 + Math.sin(state.clock.elapsedTime * 2) * 0.05)
    }
  })

  return (
    <group>
      {/* Monday's presence indicator */}
      <group ref={asteriskRef} position={[0, 1.6, -2]}>
        <Text
          fontSize={0.3}
          color="#20808D"
          anchorX="center"
          anchorY="middle"
        >
          ✱
        </Text>
        
        {/* Glow effect */}
        <Sphere args={[0.5]} position={[0, 0, 0]}>
          <meshBasicMaterial 
            color="#20808D" 
            transparent 
            opacity={0.1}
          />
        </Sphere>
      </group>

      {/* Environment lighting */}
      <ambientLight intensity={0.3} color="#20808D" />
      <directionalLight 
        position={[5, 5, 5]} 
        intensity={0.5} 
        color="#FBFAF4" 
      />

      {/* Floor grid for spatial reference */}
      <gridHelper 
        args={[20, 20, '#20808D', '#20808D']} 
        position={[0, 0, 0]}
        material-opacity={0.2}
        material-transparent={true}
      />
    </group>
  )
}

export default MondayScene 