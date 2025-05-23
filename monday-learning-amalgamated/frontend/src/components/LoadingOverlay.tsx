interface LoadingOverlayProps {
  message?: string
}

const LoadingOverlay = ({ message = "Loading..." }: LoadingOverlayProps) => {
  return (
    <div style={{
      position: 'fixed',
      top: 0,
      left: 0,
      width: '100%',
      height: '100%',
      backgroundColor: 'var(--offblack)',
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'center',
      alignItems: 'center',
      zIndex: 9999
    }}>
      <div className="loading-spinner" />
      <div style={{
        color: 'var(--paper-white)',
        fontSize: '18px',
        textAlign: 'center',
        marginTop: '20px'
      }}>
        {message}
      </div>
    </div>
  )
}

export default LoadingOverlay 