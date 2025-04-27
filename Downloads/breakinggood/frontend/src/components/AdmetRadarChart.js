import React from 'react';
import { Radar } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  RadialLinearScale,
  PointElement,
  LineElement,
  Filler,
  Tooltip,
  Legend,
} from 'chart.js';
import { Paper, Typography, makeStyles } from '@material-ui/core';

ChartJS.register(
  RadialLinearScale,
  PointElement,
  LineElement,
  Filler,
  Tooltip,
  Legend
);

const useStyles = makeStyles((theme) => ({
  paper: {
    padding: theme.spacing(2),
    height: '100%',
  },
  chartContainer: {
    height: '300px', // Adjust height as needed
    position: 'relative',
  },
}));

const AdmetRadarChart = ({ data }) => {
  const classes = useStyles();

  if (!data) {
    return (
      <Paper className={classes.paper}>
        <Typography>No ADMET data available.</Typography>
      </Paper>
    );
  }

  // --- Data Transformation --- 
  // The mock data is nested. We need to extract relevant scores for the radar chart.
  // Let's pick key ADMET properties. Adjust these labels and data extraction as needed.
  const labels = [
    'Oral Absorp.', 
    'BBB Perm.', 
    'P-gp Substrate',
    'CYP3A4 Substrate',
    'hERG Risk',
    'Hepatotoxicity'
  ];
  
  // Helper to get score safely, defaulting to 0
  const getScore = (obj, defaultValue = 0) => (obj && typeof obj.score === 'number' ? obj.score : defaultValue);
  // Helper to get risk score (assuming higher score = higher risk, needs inversion for chart)
  const getRiskScore = (obj, maxRiskScore = 100, defaultValue = 50) => {
    const score = obj && typeof obj.score === 'number' ? obj.score : defaultValue;
    // Invert risk: Lower score = better profile (higher value on chart)
    return maxRiskScore - score; 
  };

  const chartValues = [
    getScore(data.absorption?.oral), // Higher score is better
    getScore(data.absorption?.bbb),  // Higher score is better
    getRiskScore(data.absorption?.pgpSubstrate, 100, 50), // Lower score (less likely substrate) is better
    getRiskScore(data.metabolism?.cyp450Substrates?.CYP3A4, 100, 50), // Assuming higher score means substrate (less desirable)
    getRiskScore(data.toxicity?.herg), // Lower score (less risk) is better
    getRiskScore(data.toxicity?.hepatotoxicity) // Lower score (less risk) is better
  ];

  const chartData = {
    labels: labels,
    datasets: [
      {
        label: 'ADMET Profile (Higher is Better)',
        data: chartValues,
        backgroundColor: 'rgba(153, 102, 255, 0.2)',
        borderColor: 'rgba(153, 102, 255, 1)',
        borderWidth: 1,
        pointBackgroundColor: 'rgba(153, 102, 255, 1)',
        pointBorderColor: '#fff',
        pointHoverBackgroundColor: '#fff',
        pointHoverBorderColor: 'rgba(153, 102, 255, 1)'
      },
    ],
  };

  const options = {
    responsive: true,
    maintainAspectRatio: false,
    plugins: {
      legend: {
        position: 'top',
      },
      title: {
        display: true,
        text: 'ADMET Properties Overview (Normalized)',
      },
      tooltip: {
         callbacks: {
           label: function(context) {
             let label = context.dataset.label || '';
             if (label) {
               label += ': ';
             }
             const propertyName = context.label;
             // Find original data point to show classification/details if needed
             // This requires mapping back from the label, which is complex
             // For simplicity, just show the normalized score for now.
             if (context.parsed.r !== null) {
               label += context.parsed.r.toFixed(0) + ' / 100'; 
             }
             return label;
           }
         }
      }
    },
    scales: {
      r: {
        beginAtZero: true,
        max: 100,
        ticks: {
           stepSize: 20
        },
        pointLabels: {
          font: {
            size: 10 // Smaller font for labels if they overlap
          }
        }
      },
    },
  };

  return (
    <Paper className={classes.paper}>
      <div className={classes.chartContainer}>
        <Radar options={options} data={chartData} />
      </div>
    </Paper>
  );
};

export default AdmetRadarChart; 