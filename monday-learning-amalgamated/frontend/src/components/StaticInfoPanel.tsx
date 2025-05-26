import React from 'react'

interface Citation {
  id: string
  url: string
  title: string
  snippet: string
  publishedDate?: string
  domain?: string
}

interface StaticInfoPanelProps {
  id: string
  title: string
  content: string
  citations?: Citation[]
  model?: string
  isVisible?: boolean
  onClose?: () => void
}

const StaticInfoPanel: React.FC<StaticInfoPanelProps> = ({
  id,
  title,
  content,
  citations = [],
  model = 'sonar',
  isVisible = false,
  onClose
}) => {
  if (!isVisible) return null

  return (
    <div 
      className="static-info-panel"
      style={{
        position: 'fixed',
        top: '50%',
        left: '50%',
        transform: 'translate(-50%, -50%)',
        width: '600px',
        maxWidth: '90vw',
        maxHeight: '80vh',
        backgroundColor: '#FBFAF4', // Paper White
        border: '2px solid #20808D', // True Turquoise
        borderRadius: '12px',
        padding: '24px',
        color: '#091717', // Offblack
        fontFamily: 'Arial, sans-serif',
        fontSize: '14px',
        lineHeight: '1.6',
        overflowY: 'auto',
        zIndex: 1000,
        boxShadow: '0 8px 32px rgba(32, 128, 141, 0.3)',
        backdropFilter: 'blur(10px)'
      }}
    >
      {/* Header */}
      <div style={{ 
        display: 'flex', 
        justifyContent: 'space-between', 
        alignItems: 'center',
        marginBottom: '20px',
        borderBottom: '2px solid #20808D',
        paddingBottom: '12px'
      }}>
        <h3 style={{ 
          margin: 0, 
          color: '#20808D', // True Turquoise
          fontSize: '18px',
          fontWeight: 'bold'
        }}>
          {title}
        </h3>
        {onClose && (
          <button
            onClick={onClose}
            style={{
              background: 'none',
              border: 'none',
              color: '#20808D',
              fontSize: '24px',
              cursor: 'pointer',
              padding: '4px',
              width: '32px',
              height: '32px',
              borderRadius: '50%',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center'
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.backgroundColor = 'rgba(32, 128, 141, 0.1)'
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.backgroundColor = 'transparent'
            }}
          >
            ×
          </button>
        )}
      </div>

      {/* Model Indicator */}
      <div style={{
        fontSize: '12px',
        color: '#20808D',
        marginBottom: '20px',
        textAlign: 'right',
        fontWeight: '500'
      }}>
        Powered by: {model.toUpperCase()}
      </div>

      {/* Main Content */}
      <div style={{
        marginBottom: citations.length > 0 ? '24px' : '0'
      }}>
        <div style={{
          whiteSpace: 'pre-wrap',
          wordWrap: 'break-word',
          backgroundColor: 'rgba(32, 128, 141, 0.05)',
          padding: '20px',
          borderRadius: '8px',
          border: '1px solid rgba(32, 128, 141, 0.2)',
          fontSize: '15px',
          lineHeight: '1.7'
        }}>
          {content}
        </div>
      </div>

      {/* Citations */}
      {citations.length > 0 && (
        <div>
          <h4 style={{
            color: '#20808D',
            fontSize: '16px',
            marginBottom: '12px',
            borderTop: '1px solid rgba(32, 128, 141, 0.3)',
            paddingTop: '20px',
            fontWeight: 'bold'
          }}>
            Sources ({citations.length})
          </h4>
          <div style={{ fontSize: '13px' }}>
            {citations.map((citation, index) => (
              <div 
                key={citation.id || index}
                style={{
                  marginBottom: '12px',
                  padding: '12px',
                  backgroundColor: 'rgba(32, 128, 141, 0.08)',
                  borderRadius: '8px',
                  borderLeft: '4px solid #20808D'
                }}
              >
                <div style={{ fontWeight: 'bold', marginBottom: '6px', color: '#091717' }}>
                  [{index + 1}] {citation.title || 'Untitled Source'}
                </div>
                {citation.snippet && (
                  <div style={{ 
                    color: '#666', 
                    marginBottom: '6px',
                    fontSize: '12px',
                    fontStyle: 'italic'
                  }}>
                    {citation.snippet}
                  </div>
                )}
                {citation.url && (
                  <div style={{ fontSize: '11px', color: '#20808D' }}>
                    <a 
                      href={citation.url} 
                      target="_blank" 
                      rel="noopener noreferrer"
                      style={{ 
                        color: '#20808D', 
                        textDecoration: 'none',
                        fontWeight: '500'
                      }}
                    >
                      {citation.domain || citation.url}
                    </a>
                  </div>
                )}
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  )
}

export default StaticInfoPanel 