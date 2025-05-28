import React from 'react'

interface Props {
  children: React.ReactNode
}

interface State {
  hasError: boolean
  error?: Error
}

class ErrorBoundary extends React.Component<Props, State> {
  constructor(props: Props) {
    super(props)
    this.state = { hasError: false }
  }

  static getDerivedStateFromError(error: Error): State {
    return { hasError: true, error }
  }

  componentDidCatch(error: Error, errorInfo: React.ErrorInfo) {
    console.error('ErrorBoundary caught an error:', error, errorInfo)
  }

  render() {
    if (this.state.hasError) {
      return (
        <div style={{
          padding: '2rem',
          textAlign: 'center',
          color: 'var(--paper-white)',
          backgroundColor: 'var(--offblack)',
          height: '100vh',
          display: 'flex',
          flexDirection: 'column',
          justifyContent: 'center'
        }}>
          <h1>Something went wrong</h1>
          <p>Monday encountered an error and needs to restart.</p>
          <button 
            className="btn"
            onClick={() => window.location.reload()}
            style={{ margin: '1rem auto', display: 'block' }}
          >
            Restart Monday
          </button>
        </div>
      )
    }

    return this.props.children
  }
}

export default ErrorBoundary 