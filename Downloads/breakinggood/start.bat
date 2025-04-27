@echo off
echo Starting Breaking Good Application...
echo.

rem Set the Anthropic API key directly
set ANTHROPIC_API_KEY=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

echo Starting Backend Server...
start cmd /k "cd backend && npm start"

echo Waiting for backend to start...
timeout /t 5 /nobreak

echo Starting Frontend Server...
start cmd /k "cd frontend && npm start"

echo.
echo All servers started! The application should open in your browser shortly.
echo If not, navigate to: http://localhost:3000
echo.
echo Press any key to close this window...
pause > nul 