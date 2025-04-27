import React from 'react';
import { Bar } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
} from 'chart.js';
import { Paper, Typography, makeStyles } from '@material-ui/core';

ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
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

const BindingAffinityBarChart = ({ data }) => {
  const classes = useStyles();

  if (!data) {
    return (
      <Paper className={classes.paper}>
        <Typography>No binding data available.</Typography>
      </Paper>
    );
  }

  const chartData = {
    labels: Object.keys(data),
    datasets: [
      {
        label: 'Binding Score (out of 100)',
        data: Object.values(data).map(item => item.score),
        backgroundColor: 'rgba(75, 192, 192, 0.6)',
        borderColor: 'rgba(75, 192, 192, 1)',
        borderWidth: 1,
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
        text: 'Receptor Binding Affinity Scores',
      },
      tooltip: {
        callbacks: {
          label: function(context) {
            let label = context.dataset.label || '';
            if (label) {
              label += ': ';
            }
            if (context.parsed.y !== null) {
              label += context.parsed.y;
              const receptorName = context.label;
              if (data[receptorName] && data[receptorName].classification) {
                 label += ` (${data[receptorName].classification})`;
              }
            }
            return label;
          }
        }
      }
    },
    scales: {
      y: {
        beginAtZero: true,
        max: 100, // Assuming score is out of 100
        title: {
          display: true,
          text: 'Score'
        }
      },
      x: {
        title: {
          display: true,
          text: 'Receptor'
        }
      }
    },
  };

  return (
    <Paper className={classes.paper}>
      <div className={classes.chartContainer}>
        <Bar options={options} data={chartData} />
      </div>
    </Paper>
  );
};

export default BindingAffinityBarChart; 