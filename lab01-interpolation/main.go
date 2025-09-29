package main

import (
	"fmt"
	"math"
	"strings"
)

// point представляет точку (x, y)
type point struct {
	x, y float64
}

// interpolationData содержит исходные данные для интерполяции
type interpolationData struct {
	points []point // Узлы интерполяции
	a, b   float64 // Интервал [a, b]
	n      int     // Количество узлов
}

// testFunction - тестовая функция x * log10(x + 1) - 1
func testFunction(x float64) float64 {
	return x*math.Log10(x+1) - 1
}

// moduleFunction - тестовая функция модуля
func moduleFunction(x float64) float64 {
	return math.Abs(x)
}

// createGrid создает равномерную сетку точек
func createGrid(a, b float64, n int, f func(float64) float64) *interpolationData {
	h := (b - a) / float64(n)
	points := make([]point, n+1)

	for i := 0; i <= n; i++ {
		x := a + float64(i)*h
		y := f(x)
		points[i] = point{x: x, y: y}
	}

	return &interpolationData{
		points: points,
		a:      a,
		b:      b,
		n:      n,
	}
}

// createChebyshevGrid создает сетку точек на основе узлов Чебышева
func createChebyshevGrid(a, b float64, n int, f func(float64) float64) *interpolationData {
	points := make([]point, n+1)

	for i := 0; i <= n; i++ {
		// Узлы Чебышева на интервале [-1, 1]
		ti := math.Cos(math.Pi * float64(2*i+1) / float64(2*(n+1)))

		// Преобразование в интервал [a, b]
		x := (a+b)/2 + (b-a)/2*ti
		y := f(x)
		points[i] = point{x: x, y: y}
	}

	return &interpolationData{
		points: points,
		a:      a,
		b:      b,
		n:      n,
	}
}

// lagrangeInterpolation вычисляет значение интерполяционного полинома Лагранжа в точке x
func lagrangeInterpolation(data *interpolationData, x float64) float64 {
	n := len(data.points)
	result := 0.0

	for i := 0; i < n; i++ {
		// Вычисляем полином Лагранжа Li(x)
		li := 1.0
		for j := 0; j < n; j++ {
			if i != j {
				li *= (x - data.points[j].x) / (data.points[i].x - data.points[j].x)
			}
		}
		result += data.points[i].y * li
	}

	return result
}

type matrix struct {
	data [][]float64
	rows int
	cols int
}

// newMatrix создает новую матрицу
func newMatrix(rows, cols int) *matrix {
	data := make([][]float64, rows)
	for i := range data {
		data[i] = make([]float64, cols)
	}
	return &matrix{data: data, rows: rows, cols: cols}
}

func (m *matrix) set(i, j int, val float64) {
	m.data[i][j] = val
}

func (m *matrix) get(i, j int) float64 {
	return m.data[i][j]
}

// solveLinearSystem решает систему линейных уравнений Ax = b методом Гаусса
func solveLinearSystem(a *matrix, b []float64) []float64 {
	n := a.rows

	// Создаем расширенную матрицу
	augmented := newMatrix(n, n+1)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			augmented.set(i, j, a.get(i, j))
		}
		augmented.set(i, n, b[i])
	}

	// Прямой ход метода Гаусса
	for i := 0; i < n; i++ {
		// Приведение к верхнетреугольному виду
		for k := i + 1; k < n; k++ {
			if math.Abs(augmented.get(i, i)) < 1e-12 {
				continue
			}
			factor := augmented.get(k, i) / augmented.get(i, i)
			for j := i; j <= n; j++ {
				augmented.set(k, j, augmented.get(k, j)-factor*augmented.get(i, j))
			}
		}
	}

	// Обратный ход
	solution := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		solution[i] = augmented.get(i, n)
		for j := i + 1; j < n; j++ {
			solution[i] -= augmented.get(i, j) * solution[j]
		}
		if math.Abs(augmented.get(i, i)) > 1e-12 {
			solution[i] /= augmented.get(i, i)
		}
	}

	return solution
}

// cubicSpline представляет кубический сплайн с прямым вычислением по формуле
type cubicSpline struct {
	points            []point
	secondDerivatives []float64
	h                 []float64
}

// newCubicSpline создает кубический сплайн с естественными граничными условиями
func newCubicSpline(data *interpolationData) *cubicSpline {
	points := data.points
	n := len(points)

	// Извлекаем x и y координаты
	x := make([]float64, n)
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = points[i].x
		y[i] = points[i].y
	}

	// Вычисляем h[i] = x[i+1] - x[i]
	h := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = x[i+1] - x[i]
	}

	// Создаем матрицу a и вектор b для системы уравнений
	a := newMatrix(n, n)
	b := make([]float64, n)

	// Заполняем систему уравнений для внутренних точек
	for i := 1; i < n-1; i++ {
		a.set(i, i-1, h[i-1])
		a.set(i, i, 2*(h[i-1]+h[i]))
		a.set(i, i+1, h[i])
		b[i] = 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
	}

	// Граничные условия для естественного сплайна (вторые производные на концах равны нулю)
	a.set(0, 0, 1)
	a.set(n-1, n-1, 1)
	b[0] = 0
	b[n-1] = 0

	// Решаем систему для вторых производных
	secondDerivatives := solveLinearSystem(a, b)

	return &cubicSpline{
		points:            points,
		secondDerivatives: secondDerivatives,
		h:                 h,
	}
}

// Evaluate вычисляет значение сплайна в точке x по формуле (2.61)
func (cs *cubicSpline) evaluate(x float64) float64 {
	n := len(cs.points)

	// Находим интервал, содержащий точку x
	i := 0
	for i < n-1 {
		if x >= cs.points[i].x && x <= cs.points[i+1].x {
			break
		}
		i++
	}

	// формула (2.61)
	xi := cs.points[i].x
	xi1 := cs.points[i+1].x
	yi := cs.points[i].y
	yi1 := cs.points[i+1].y
	hi1 := cs.h[i]
	gammai := cs.secondDerivatives[i]
	gammai1 := cs.secondDerivatives[i+1]

	term1 := yi * (xi1 - x) / hi1
	term2 := yi1 * (x - xi) / hi1

	xi1minusx := xi1 - x
	xminusxi := x - xi

	term3 := gammai * (xi1minusx*xi1minusx*xi1minusx - hi1*hi1*xi1minusx) / (6 * hi1)
	term4 := gammai1 * (xminusxi*xminusxi*xminusxi - hi1*hi1*xminusxi) / (6 * hi1)

	return term1 + term2 + term3 + term4
}

// printTable выводит таблицу исходных данных
func printTable(data *interpolationData, title string) {
	fmt.Printf("Таблица исходных данных (%s):\n", title)
	fmt.Printf("%-10s %-15s\n", "xi", "f(xi)")
	fmt.Println(strings.Repeat("-", 25))

	for _, point := range data.points {
		fmt.Printf("%-10.4f %-15.6f\n", point.x, point.y)
	}
	fmt.Println()
}

// compareInterpolations сравнивает методы интерполяции
func compareInterpolations(uniformData, chebyshevData *interpolationData, testFunc func(float64) float64) {
	fmt.Println("Сравнение методов интерполяции:")
	fmt.Printf("%-10s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n",
		"x", "f(x)", "пол Лагр", "Ош Лагр", "узлы Чеб", "Ош Чеб", "Сплайн", "Ош Спл")
	fmt.Println(strings.Repeat("-", 110))

	spline := newCubicSpline(uniformData)

	for i := 0; i < 20; i++ {
		x := uniformData.a + float64(i)*(uniformData.b-uniformData.a)/19.0

		original := testFunc(x)
		lagrangeUniform := lagrangeInterpolation(uniformData, x)
		lagrangeChebyshev := lagrangeInterpolation(chebyshevData, x)
		splineVal := spline.evaluate(x)

		errorUniform := math.Abs(original - lagrangeUniform)
		errorChebyshev := math.Abs(original - lagrangeChebyshev)
		errorSpline := math.Abs(original - splineVal)

		fmt.Printf("%-10.4f %-12.6f %-12.6f %-12.6e %-12.6f %-12.6e %-12.6f %-12.6e\n",
			x, original, lagrangeUniform, errorUniform, lagrangeChebyshev, errorChebyshev, splineVal, errorSpline)
	}
	fmt.Println()

	maxErrorUniform := 0.0
	maxErrorChebyshev := 0.0
	maxErrorSpline := 0.0

	for i := 0; i < 100; i++ {
		x := uniformData.a + float64(i)*(uniformData.b-uniformData.a)/99.0
		original := testFunc(x)

		errorUniform := math.Abs(original - lagrangeInterpolation(uniformData, x))
		errorChebyshev := math.Abs(original - lagrangeInterpolation(chebyshevData, x))
		errorSpline := math.Abs(original - spline.evaluate(x))

		if errorUniform > maxErrorUniform {
			maxErrorUniform = errorUniform
		}
		if errorChebyshev > maxErrorChebyshev {
			maxErrorChebyshev = errorChebyshev
		}
		if errorSpline > maxErrorSpline {
			maxErrorSpline = errorSpline
		}
	}

	fmt.Println("Максимальные ошибки:")
	fmt.Printf("  Лагранж (равномерные узлы): %.6e\n", maxErrorUniform)
	fmt.Printf("  Лагранж (узлы Чебышева):    %.6e\n", maxErrorChebyshev)
	fmt.Printf("  Кубический сплайн:          %.6e\n", maxErrorSpline)
	fmt.Println()
}

func main() {
	fmt.Printf("=== Лабораторная работа №1: Интерполяция ===\n")

	// Параметры для интерполяции
	a, b := 1.0, 5.0

	// Тестирование с разным количеством узлов
	nValues := []int{10}

	for _, n := range nValues {
		fmt.Printf("\n=== Тестирование с N = %d узлами ===\n\n", n)

		// Создаем равномерную сетку
		uniformData := createGrid(a, b, n, testFunction)
		printTable(uniformData, "равномерные узлы")

		// Создаем сетку Чебышева
		chebyshevData := createChebyshevGrid(a, b, n, testFunction)
		printTable(chebyshevData, "узлы Чебышева")

		// Сравниваем методы интерполяции
		compareInterpolations(uniformData, chebyshevData, testFunction)

		// Генерируем HTML файл с графиками
		filename := fmt.Sprintf("interpolation_n%d.html", n)
		err := generateHTML(uniformData, chebyshevData, testFunction, filename)
		if err != nil {
			fmt.Printf("Ошибка при создании HTML файла: %v\n", err)
		} else {
			fmt.Printf("✓ График сохранен в файл: %s\n\n", filename)
		}
	}

	fmt.Println("Все графики созданы! Откройте HTML файлы в браузере для просмотра.")
}
