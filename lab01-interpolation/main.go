package main

import (
	"fmt"
	"math"
	"strings"
)

// Point представляет точку (x, y)
type Point struct {
	X, Y float64
}

// InterpolationData содержит исходные данные для интерполяции
type InterpolationData struct {
	Points []Point // Узлы интерполяции
	A, B   float64 // Интервал [A, B]
	N      int     // Количество узлов
}

// TestFunction - тестовая функция x * log10(x + 1) - 1
func TestFunction(x float64) float64 {
	return x*math.Log10(x+1) - 1
}

// CreateGrid создает равномерную сетку точек
func CreateGrid(a, b float64, n int, f func(float64) float64) *InterpolationData {
	h := (b - a) / float64(n)
	points := make([]Point, n+1)

	for i := 0; i <= n; i++ {
		x := a + float64(i)*h
		y := f(x)
		points[i] = Point{X: x, Y: y}
	}

	return &InterpolationData{
		Points: points,
		A:      a,
		B:      b,
		N:      n,
	}
}

// LagrangeInterpolation вычисляет значение интерполяционного полинома Лагранжа в точке x
func LagrangeInterpolation(data *InterpolationData, x float64) float64 {
	n := len(data.Points)
	result := 0.0

	for i := 0; i < n; i++ {
		// Вычисляем полином Лагранжа Li(x)
		li := 1.0
		for j := 0; j < n; j++ {
			if i != j {
				li *= (x - data.Points[j].X) / (data.Points[i].X - data.Points[j].X)
			}
		}
		result += data.Points[i].Y * li
	}

	return result
}

// Matrix представляет матрицу для решения системы линейных уравнений
type Matrix struct {
	data [][]float64
	rows int
	cols int
}

// NewMatrix создает новую матрицу
func NewMatrix(rows, cols int) *Matrix {
	data := make([][]float64, rows)
	for i := range data {
		data[i] = make([]float64, cols)
	}
	return &Matrix{data: data, rows: rows, cols: cols}
}

// Set устанавливает значение элемента матрицы
func (m *Matrix) Set(i, j int, val float64) {
	m.data[i][j] = val
}

// Get получает значение элемента матрицы
func (m *Matrix) Get(i, j int) float64 {
	return m.data[i][j]
}

// SolveLinearSystem решает систему линейных уравнений Ax = b методом Гаусса
func SolveLinearSystem(A *Matrix, b []float64) []float64 {
	n := A.rows

	// Создаем расширенную матрицу
	augmented := NewMatrix(n, n+1)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			augmented.Set(i, j, A.Get(i, j))
		}
		augmented.Set(i, n, b[i])
	}

	// Прямой ход
	for i := 0; i < n; i++ {
		// Приведение к верхнетреугольному виду
		for k := i + 1; k < n; k++ {
			factor := augmented.Get(k, i) / augmented.Get(i, i)
			for j := i; j <= n; j++ {
				augmented.Set(k, j, augmented.Get(k, j)-factor*augmented.Get(i, j))
			}
		}
	}

	// Обратный ход
	solution := make([]float64, n) // массив решений x
	for i := n - 1; i >= 0; i-- {
		solution[i] = augmented.Get(i, n)
		for j := i + 1; j < n; j++ {
			solution[i] -= augmented.Get(i, j) * solution[j]
		}
		solution[i] /= augmented.Get(i, i)
	}

	return solution
}

// SplineCoefficient представляет коэффициенты одного сегмента сплайна
type SplineCoefficient struct {
	A, B, C, D float64
}

// CubicSpline представляет кубический сплайн
type CubicSpline struct {
	Points       []Point
	Coefficients []SplineCoefficient
}

// NewCubicSpline создает кубический сплайн с естественными граничными условиями
func NewCubicSpline(data *InterpolationData) *CubicSpline {
	points := data.Points
	n := len(points)

	// Извлекаем x и y координаты
	x := make([]float64, n)
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = points[i].X
		y[i] = points[i].Y
	}

	// Вычисляем h[i] = x[i] - x[i-1]
	h := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = x[i+1] - x[i]
	}

	// Создаем матрицу A и вектор B для системы уравнений
	A := NewMatrix(n, n)
	B := make([]float64, n)

	// Заполняем систему уравнений для внутренних точек
	for i := 1; i < n-1; i++ {
		A.Set(i, i-1, h[i-1])
		A.Set(i, i, 2*(h[i-1]+h[i]))
		A.Set(i, i+1, h[i])
		B[i] = 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
	}

	// Граничные условия для естественного сплайна (вторые производные на концах равны нулю)
	A.Set(0, 0, 1)
	A.Set(n-1, n-1, 1)
	B[0] = 0
	B[n-1] = 0

	// Решаем систему для вторых производных
	secondDerivatives := SolveLinearSystem(A, B)

	// Вычисляем коэффициенты полиномов
	coefficients := make([]SplineCoefficient, n-1)
	for i := 0; i < n-1; i++ {
		a := (secondDerivatives[i+1] - secondDerivatives[i]) / (6 * h[i])
		b := secondDerivatives[i] / 2
		c := (y[i+1]-y[i])/h[i] - (2*h[i]*secondDerivatives[i]+h[i]*secondDerivatives[i+1])/6
		d := y[i]

		coefficients[i] = SplineCoefficient{A: a, B: b, C: c, D: d}
	}

	return &CubicSpline{
		Points:       points,
		Coefficients: coefficients,
	}
}

// Evaluate вычисляет значение сплайна в точке t
func (cs *CubicSpline) Evaluate(t float64) float64 {
	n := len(cs.Points)

	// Находим подходящий интервал
	i := 0
	for i < n-1 {
		if t >= cs.Points[i].X && t <= cs.Points[i+1].X {
			break
		}
		i++
	}

	// Вычисляем значение полинома
	coeff := cs.Coefficients[i]
	deltaX := t - cs.Points[i].X

	return coeff.A*deltaX*deltaX*deltaX + coeff.B*deltaX*deltaX + coeff.C*deltaX + coeff.D
}

// PrintTable выводит таблицу исходных данных
func PrintTable(data *InterpolationData) {
	fmt.Println("Таблица исходных данных:")
	fmt.Printf("%-10s %-15s\n", "xi", "f(xi)")
	fmt.Println(strings.Repeat("-", 25))

	for _, point := range data.Points {
		fmt.Printf("%-10.4f %-15.6f\n", point.X, point.Y)
	}
	fmt.Println()
}

// CompareInterpolations сравнивает методы интерполяции
func CompareInterpolations(data *InterpolationData, testFunc func(float64) float64) {
	fmt.Println("Сравнение методов интерполяции:")
	fmt.Printf("%-10s %-15s %-15s %-15s %-15s %-15s %-15s\n", "x", "Исходная f(x)", "Лагранж", "Ошибка Лагр", "Сплайн", "Ошибка Спл", "Кто точнее")
	fmt.Println(strings.Repeat("-", 110))

	spline := NewCubicSpline(data)

	for i := 0; i < 20; i++ {
		x := data.A + float64(i)*(data.B-data.A)/19.0

		original := testFunc(x)
		lagrange := LagrangeInterpolation(data, x)
		splineVal := spline.Evaluate(x)

		errorLagrange := math.Abs(original - lagrange)
		errorSpline := math.Abs(original - splineVal)

		var moreAccurate string
		if errorLagrange < errorSpline {
			moreAccurate = "Лагр"
		} else if errorSpline < errorLagrange {
			moreAccurate = "Спл"
		} else {
			moreAccurate = "Одинаково"
		}

		fmt.Printf("%-10.4f %-15.6f %-15.6f %-15.6f %-15.6f %-15.6f %-15s\n", x, original, lagrange, errorLagrange, splineVal, errorSpline, moreAccurate)
	}
	fmt.Println()
}

func main() {
	fmt.Printf("=== Лабораторная работа №1: Интерполяция ===\n")

	// Параметры для интерполяции
	a, b := 1.0, 6.0

	// Тестирование с разным количеством узлов
	nValues := []int{5, 10, 15}

	for _, n := range nValues {
		fmt.Printf("=== Тестирование с N = %d узлами ===\n", n)

		// Создаем данные с тестовой функцией x * log10(x + 1) - 1
		fmt.Printf("Функция: f(x) = x * log10(x + 1) - 1, интервал [%.1f, %.1f]\n\n", a, b)

		data := CreateGrid(a, b, n, TestFunction)
		PrintTable(data)

		// Сравниваем методы интерполяции
		CompareInterpolations(data, TestFunction)
	}
}
