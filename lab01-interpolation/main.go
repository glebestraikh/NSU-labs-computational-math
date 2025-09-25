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

// CubicSpline представляет кубический сплайн
type CubicSpline struct {
	Points     []Point
	A, B, C, D []float64 // Коэффициенты сплайна
}

// NewCubicSpline создает кубический сплайн с естественными граничными условиями
func NewCubicSpline(data *InterpolationData) *CubicSpline {
	n := len(data.Points) - 1
	points := data.Points

	// Шаги между узлами
	h := make([]float64, n)
	for i := 0; i < n; i++ {
		h[i] = points[i+1].X - points[i].X
	}

	// Коэффициенты для системы уравнений
	alpha := make([]float64, n)
	for i := 1; i < n; i++ {
		alpha[i] = 3.0/h[i]*(points[i+1].Y-points[i].Y) - 3.0/h[i-1]*(points[i].Y-points[i-1].Y)
	}

	// Решение системы уравнений методом прогонки
	l := make([]float64, n+1)
	mu := make([]float64, n+1)
	z := make([]float64, n+1)

	l[0] = 1.0  // l_i - диагональные элементы после преобразования
	mu[0] = 0.0 // μ_i - наддиагональные элементы
	z[0] = 0.0  // z_i - правая часть после преобразования

	for i := 1; i < n; i++ {
		l[i] = 2.0*(points[i+1].X-points[i-1].X) - h[i-1]*mu[i-1]
		mu[i] = h[i] / l[i]
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
	}

	l[n] = 1.0
	z[n] = 0.0

	// Обратный ход
	c := make([]float64, n+1)
	b := make([]float64, n)
	d := make([]float64, n)
	a := make([]float64, n)

	c[n] = 0.0

	for j := n - 1; j >= 0; j-- {
		c[j] = z[j] - mu[j]*c[j+1]
		b[j] = (points[j+1].Y-points[j].Y)/h[j] - h[j]*(c[j+1]+2.0*c[j])/3.0
		d[j] = (c[j+1] - c[j]) / (3.0 * h[j])
		a[j] = points[j].Y
	}

	return &CubicSpline{
		Points: points,
		A:      a,
		B:      b,
		C:      c,
		D:      d,
	}
}

// Evaluate вычисляет значение сплайна в точке x
func (cs *CubicSpline) Evaluate(x float64) float64 {
	// Находим интервал, содержащий x
	n := len(cs.Points) - 1
	j := 0

	for i := 0; i < n; i++ {
		if x >= cs.Points[i].X && x <= cs.Points[i+1].X {
			j = i
			break
		}
	}

	// Если x вне интервала, используем ближайший сегмент
	if x < cs.Points[0].X {
		j = 0
	} else if x > cs.Points[n].X {
		j = n - 1
	}

	dx := x - cs.Points[j].X
	return cs.A[j] + cs.B[j]*dx + cs.C[j]*dx*dx + cs.D[j]*dx*dx*dx
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
	fmt.Printf("%-10s %-15s %-15s %-15s %-15s %-15s\n", "x", "Исходная f(x)", "Лагранж", "Ошибка Лагр", "Сплайн", "Ошибка Спл")
	fmt.Println(strings.Repeat("-", 90))

	spline := NewCubicSpline(data)

	for i := 0; i < 20; i++ {
		x := data.A + float64(i)*(data.B-data.A)/19.0

		original := testFunc(x)
		lagrange := LagrangeInterpolation(data, x)
		splineVal := spline.Evaluate(x)

		errorLagrange := math.Abs(original - lagrange)
		errorSpline := math.Abs(original - splineVal)

		fmt.Printf("%-10.4f %-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", x, original, lagrange, errorLagrange, splineVal, errorSpline)
	}
	fmt.Println()
}

func main() {
	fmt.Printf("=== Лабораторная работа №1: Интерполяция ===\n")

	// Параметры для интерполяции
	a, b := 1.0, 6.0

	// Тестирование с разным количеством узлов
	nValues := []int{5}

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
