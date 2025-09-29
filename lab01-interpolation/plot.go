package main

import (
	"fmt"
	"math"
	"os"
	"strings"
)

// generateHTML создает HTML файл с графиками
func generateHTML(uniformData, chebyshevData *interpolationData, testFunc func(float64) float64, filename string) error {
	spline := newCubicSpline(uniformData)

	// Генерируем данные для графиков
	numPoints := 200
	step := (uniformData.b - uniformData.a) / float64(numPoints)

	var xValues, originalValues, lagrangeUniformValues, lagrangeChebyshevValues, splineValues []float64
	var lagrangeUniformErrors, lagrangeChebyshevErrors, splineErrors []float64

	for i := 0; i <= numPoints; i++ {
		x := uniformData.a + float64(i)*step
		original := testFunc(x)
		lagrangeUniform := lagrangeInterpolation(uniformData, x)
		lagrangeChebyshev := lagrangeInterpolation(chebyshevData, x)
		splineVal := spline.evaluate(x)

		xValues = append(xValues, x)
		originalValues = append(originalValues, original)
		lagrangeUniformValues = append(lagrangeUniformValues, lagrangeUniform)
		lagrangeChebyshevValues = append(lagrangeChebyshevValues, lagrangeChebyshev)
		splineValues = append(splineValues, splineVal)
		lagrangeUniformErrors = append(lagrangeUniformErrors, math.Abs(original-lagrangeUniform))
		lagrangeChebyshevErrors = append(lagrangeChebyshevErrors, math.Abs(original-lagrangeChebyshev))
		splineErrors = append(splineErrors, math.Abs(original-splineVal))
	}

	// Конвертируем данные в JSON формат
	xValuesStr := floatSliceToJS(xValues)
	originalValuesStr := floatSliceToJS(originalValues)
	lagrangeUniformValuesStr := floatSliceToJS(lagrangeUniformValues)
	lagrangeChebyshevValuesStr := floatSliceToJS(lagrangeChebyshevValues)
	splineValuesStr := floatSliceToJS(splineValues)
	lagrangeUniformErrorsStr := floatSliceToJS(lagrangeUniformErrors)
	lagrangeChebyshevErrorsStr := floatSliceToJS(lagrangeChebyshevErrors)
	splineErrorsStr := floatSliceToJS(splineErrors)

	// Данные узлов (равномерные)
	var uniformNodesX, uniformNodesY []float64
	for _, p := range uniformData.points {
		uniformNodesX = append(uniformNodesX, p.x)
		uniformNodesY = append(uniformNodesY, p.y)
	}
	uniformNodesXStr := floatSliceToJS(uniformNodesX)
	uniformNodesYStr := floatSliceToJS(uniformNodesY)

	// Данные узлов (Чебышева)
	var chebyshevNodesX, chebyshevNodesY []float64
	for _, p := range chebyshevData.points {
		chebyshevNodesX = append(chebyshevNodesX, p.x)
		chebyshevNodesY = append(chebyshevNodesY, p.y)
	}
	chebyshevNodesXStr := floatSliceToJS(chebyshevNodesX)
	chebyshevNodesYStr := floatSliceToJS(chebyshevNodesY)

	htmlContent := fmt.Sprintf(`<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Результаты интерполяции</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/3.9.1/chart.min.js"></script>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 1600px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }
        h1 {
            text-align: center;
            color: #333;
        }
        .charts-container {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-bottom: 20px;
        }
        .chart-container {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .full-width {
            grid-column: 1 / -1;
        }
        canvas {
            max-width: 100%%;
            height: 400px !important;
        }
        h2 {
            margin-top: 0;
            color: #555;
        }
    </style>
</head>
<body>
    <h1>Результаты интерполяции (N = %d узлов)</h1>
    
    <div class="charts-container">
        <div class="chart-container full-width">
            <h2>Сравнение методов интерполяции</h2>
            <canvas id="interpolationChart"></canvas>
        </div>
        
        <div class="chart-container">
            <h2>Равномерные узлы</h2>
            <canvas id="uniformNodesChart"></canvas>
        </div>
        
        <div class="chart-container">
            <h2>Узлы Чебышева</h2>
            <canvas id="chebyshevNodesChart"></canvas>
        </div>
        
        <div class="chart-container full-width">
            <h2>Сравнение ошибок интерполяции</h2>
            <canvas id="errorChart"></canvas>
        </div>
    </div>

    <script>
        // График интерполяции
        const ctx1 = document.getElementById('interpolationChart').getContext('2d');
        new Chart(ctx1, {
            type: 'line',
            data: {
                labels: %s,
                datasets: [{
                    label: 'Исходная функция',
                    data: %s,
                    borderColor: 'rgb(75, 192, 192)',
                    borderWidth: 3,
                    pointRadius: 0,
                    tension: 0.1
                }, {
                    label: 'Лагранж (равномерные узлы)',
                    data: %s,
                    borderColor: 'rgb(255, 99, 132)',
                    borderWidth: 2,
                    borderDash: [5, 5],
                    pointRadius: 0,
                    tension: 0.1
                }, {
                    label: 'Лагранж (узлы Чебышева)',
                    data: %s,
                    borderColor: 'rgb(153, 102, 255)',
                    borderWidth: 2,
                    borderDash: [10, 5],
                    pointRadius: 0,
                    tension: 0.1
                }, {
                    label: 'Кубический сплайн',
                    data: %s,
                    borderColor: 'rgb(54, 162, 235)',
                    borderWidth: 2,
                    borderDash: [2, 2],
                    pointRadius: 0,
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { position: 'top' }
                },
                scales: {
                    x: { title: { display: true, text: 'x' } },
                    y: { title: { display: true, text: 'f(x)' } }
                }
            }
        });

        // График равномерных узлов
        const ctx2 = document.getElementById('uniformNodesChart').getContext('2d');
        new Chart(ctx2, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Равномерные узлы',
                    data: %s.map((x, i) => ({x: x, y: %s[i]})),
                    borderColor: 'rgb(255, 99, 132)',
                    backgroundColor: 'rgba(255, 99, 132, 0.8)',
                    pointRadius: 6
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { position: 'top' }
                },
                scales: {
                    x: { title: { display: true, text: 'x' } },
                    y: { title: { display: true, text: 'f(x)' } }
                }
            }
        });

        // График узлов Чебышева
        const ctx3 = document.getElementById('chebyshevNodesChart').getContext('2d');
        new Chart(ctx3, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Узлы Чебышева',
                    data: %s.map((x, i) => ({x: x, y: %s[i]})),
                    borderColor: 'rgb(153, 102, 255)',
                    backgroundColor: 'rgba(153, 102, 255, 0.8)',
                    pointRadius: 6
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { position: 'top' }
                },
                scales: {
                    x: { title: { display: true, text: 'x' } },
                    y: { title: { display: true, text: 'f(x)' } }
                }
            }
        });

        // График ошибок
        const ctx4 = document.getElementById('errorChart').getContext('2d');
        new Chart(ctx4, {
            type: 'line',
            data: {
                labels: %s,
                datasets: [{
                    label: 'Ошибка Лагранжа (равномерные)',
                    data: %s,
                    borderColor: 'rgb(255, 99, 132)',
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0.1
                }, {
                    label: 'Ошибка Лагранжа (Чебышев)',
                    data: %s,
                    borderColor: 'rgb(153, 102, 255)',
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0.1
                }, {
                    label: 'Ошибка сплайна',
                    data: %s,
                    borderColor: 'rgb(54, 162, 235)',
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { position: 'top' }
                },
                scales: {
                    x: { title: { display: true, text: 'x' } },
                    y: { 
                        type: 'logarithmic',
                        title: { display: true, text: 'Ошибка (log)' } 
                    }
                }
            }
        });
    </script>
</body>
</html>`, uniformData.n, xValuesStr, originalValuesStr, lagrangeUniformValuesStr, lagrangeChebyshevValuesStr,
		splineValuesStr, uniformNodesXStr, uniformNodesYStr, chebyshevNodesXStr, chebyshevNodesYStr,
		xValuesStr, lagrangeUniformErrorsStr, lagrangeChebyshevErrorsStr, splineErrorsStr)

	return os.WriteFile(filename, []byte(htmlContent), 0644)
}

// floatSliceToJS конвертирует срез float64 в JavaScript массив
func floatSliceToJS(values []float64) string {
	var result strings.Builder
	result.WriteString("[")
	for i, v := range values {
		if i > 0 {
			result.WriteString(",")
		}
		result.WriteString(fmt.Sprintf("%.6f", v))
	}
	result.WriteString("]")
	return result.String()
}
