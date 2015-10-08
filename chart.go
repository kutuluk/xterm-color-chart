package main

import (
	"fmt"
	"math"

	"github.com/kutuluk/xterm-color-chart/color"
	"github.com/nsf/termbox-go"
)

// обертка над termbox'ом
type Screen struct {
	Buffer []termbox.Cell
	Width  int
	Height int
}

func (s *Screen) reScreen() {
	s.Buffer = termbox.CellBuffer()
	s.Width, s.Height = termbox.Size()
}

func (s *Screen) Init() error {
	err := termbox.Init()
	if err == nil {
		s.reScreen()
	}
	return err
}

func (s *Screen) Flush() {
	termbox.Flush()
	s.reScreen()
}

func (s *Screen) Clear(fg, bg termbox.Attribute) {
	termbox.Clear(fg, bg)
	s.reScreen()
}

var screen Screen

func termNative(c int) termbox.Attribute {
	return termbox.Attribute(c + 1)
}

func printString(x, y int, s string) {
	// TODO срабатывает перенос на следующую строку при достижении правого края
	// и ошибка при выходе из границ буфера снизу
	offsetX := 0
	offsetY := 0
	for _, char := range s {
		if char == '\n' {
			offsetX = 0
			offsetY++
		} else {
			screen.Buffer[x+offsetX+(y+offsetY)*screen.Width].Ch = char
			offsetX++
		}
	}
}

func printLine(x, y int, s string) {
	printLineAttr(x, y, s, termbox.ColorWhite, background)
}

func printLineAttr(x, y int, s string, fg, bg termbox.Attribute) {
	offsetX := 0
	offsetY := 0
	for _, char := range s {
		if char == '\n' {
			offsetX = 0
			offsetY++
		} else {
			termbox.SetCell(x+offsetX, y+offsetY, char, fg, bg)
			offsetX++
		}
	}
}

func DrawColor(x, y int, c int) {

	Grayscale := []byte{
		16, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243,
		244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 231,
	}

	for xx := 0; xx < 13; xx++ {
		termbox.SetCell(x+xx, y, ' ', termbox.ColorBlack, termbox.Attribute(c))
	}

	// round(36 * (r * 5) + 6 * (g * 5) + (b * 5) + 16)
	// https://gist.github.com/MicahElliott/719710

	v := int(0.2125*float64(color.XtermRGBPalette[c].R) + 0.7154*float64(color.XtermRGBPalette[c].G) + 0.0721*float64(color.XtermRGBPalette[c].B))

	fg := termbox.ColorRed
	//	if v < 128 {
	//		fg = termbox.ColorWhite
	//fg = termbox.Attribute(231)
	delta := 64
	center := 128
	if (v > (center - delta)) && (v <= center) {
		fg = termbox.Attribute(Grayscale[25])
	} else if (v > center) && (v < (center + delta)) {
		fg = termbox.Attribute(Grayscale[0])
	} else {
		//fg = termbox.Attribute(Grayscale[12])
		fg = termbox.Attribute(Grayscale[25-int(v/11)])
	}

	printLineAttr(x+1, y, fmt.Sprintf("%03d #%s", int(c), color.XtermRGBPalette[c]), fg, termbox.Attribute(c))
}

func DrawBox(left, top, d, l int) {
	var approxer string
	switch currentComparer {
	case 0:
		approxer = "CIE76"
	case 1:
		approxer = "CIE94"
	case 2:
		approxer = "CIE2000"
	}

	r := d / 2
	for y := 0; y < d; y++ {
		for x := 0; x < d; x++ {
			a := 100.0 * float64(x-r) / float64(r)
			b := -100.0 * float64(y-r) / float64(r)
			_, c := color.Approximate(color.LabColor{float64(l), a, b}, comparers[currentComparer])
			//termbox.SetCell(x*2+left, y+top, ' ', background, termbox.Attribute(color))
			//termbox.SetCell(x*2+left+1, y+top, ' ', background, termbox.Attribute(color))
			fg := int(background)
			if y == 1 || y == 2 || y == d-2 {
				fgLab := color.XtermRGBPalette[c-17].Lab()
				if l > 50 {
					fgLab.L -= 50
					//fgLab.L = 20
				} else {
					fgLab.L += 50
					//fgLab.L = 70
				}
				_, fg = color.Approximate(fgLab, color.DeltaCIE2000{})
			}
			termbox.SetCell(x*2+left, y+top, ' ', termbox.Attribute(fg), termbox.Attribute(c))
			termbox.SetCell(x*2+left+1, y+top, ' ', termbox.Attribute(fg), termbox.Attribute(c))
		}

	}

	labLabel := "CIE L*a*b color space"
	if len(labLabel)+4 < d*2 {
		printString(left+2, top+1, labLabel)
		printString(left+2, top+2, fmt.Sprintf("Lightness = %d%%", l))
	}
	//printString(left+2, top+4, fmt.Sprintln("Press Up/Down key for change"))

	approxerLabel := fmt.Sprintf("Approximate by %s algorithm", approxer)
	if len(approxerLabel)+4 <= d*2 {
		printString(left+2, top+d-2, approxerLabel)
	}
	//printString(left+2, top+7, fmt.Sprintln("Press Left/Right key for change"))
}

func DrawCircle(left, top, r int, l float64) {

	var ch [2]rune
	point := func(xb, yb int) {
		xl := xb / r
		yl := yb / r
		distance, color := color.Approximate(color.LChColor{l, math.Sqrt(float64(xl*xl) + float64(yl*yl)), 180.0 * math.Atan2(float64(xl), float64(yl))}, comparers[currentComparer])
		//distance, color := color.Palette(color.LabColor{l, float64(xb * 120.0 / r), float64(-yb * 120.0 / r)}, color.DeltaCIE94{})
		if distance < 50 {
			ch[0] = ' '
			ch[1] = ' '
		} else {
			//ch[0] = '¤'
			//ch[1] = '¤'
			ch[0] = '('
			ch[1] = ')'
		}
		termbox.SetCell((xb+r)*2+left, yb+r+top, ch[0], termbox.Attribute(0), termbox.Attribute(color))
		termbox.SetCell((xb+r)*2+left+1, yb+r+top, ch[1], termbox.Attribute(0), termbox.Attribute(color))
	}

	// алгоритм Брезенхэма
	line := func(x, y1, y2 int) {
		for yf := y1; yf <= y2; yf++ {
			point(x, yf)
		}
	}

	x := 0
	y := r
	d := 3 - 2*r
	for x <= y {
		line(x, -y, y)
		line(-x, -y, y)
		line(y, -x, x)
		line(-y, -x, x)
		if d < 0 {
			d += 4*x + 6
		} else {
			d += 4*(x-y) + 10
			y--
		}
		x++
	}
}

func DrawColorCell(x, y, c int) {
	termbox.SetCell(x, y, ' ', termbox.Attribute(1), termbox.Attribute(c))
	termbox.SetCell(x+1, y, ' ', termbox.Attribute(1), termbox.Attribute(c))
}

// left, top - левый верхний угол отрисовки палитры
// row - количество кубов 6x6 по горизонтали
// num - признак отрисовки номеров цветов
func DrawPalette(left, top, row int, num bool) {
	var numOffset int
	if num {
		numOffset = 4
	}

	//bh := [2]termbox.Attribute{237, background}
	bh := [2]termbox.Attribute{236, background}
	//th := [2]termbox.Attribute{232, 241}
	th := [2]termbox.Attribute{243, 243}

	cc := 1

	printLineAttr(left, top, fmt.Sprintln("16 system colors"), text, background)

	printLineAttr(left, top+2, fmt.Sprintf("%3d ", cc-1), th[1], bh[1])
	printLineAttr(left+numOffset+16*2, top+2, fmt.Sprintf(" %-3d", cc+14), th[1], bh[1])

	for x := 0; x < 16; x++ {
		DrawColorCell(left+numOffset+x*2, top+2, cc)
		cc++
	}

	padding := 1
	if num {
		//		padding = 0
	}

	printLineAttr(left, top+4, fmt.Sprintln("216 colors of color cube 6x6x6"), text, background)
	col := int(6 / row)        // ширина
	for l := 0; l < row; l++ { // l - вертикальный счетчик
		for k := 0; k < col; k++ { // k - горизонтальный счетчик
			for y := 0; y < 6; y++ {
				h := y % 2
				if h != 0 {
					h = 1
				}
				printLineAttr(left+k*numOffset*2+k*(6+padding)*2, top+l*(6+padding)+y+6, fmt.Sprintf("%3d ", cc-1), th[h], bh[h])
				printLineAttr(left+k*numOffset*2+k*(6+padding)*2+16, top+l*(6+padding)+y+6, fmt.Sprintf(" %-3d", cc+4), th[h], bh[h])
				for x := 0; x < 6; x++ {
					DrawColorCell(left+numOffset+k*(6+padding)*2+k*numOffset*2+x*2, top+l*(6+padding)+y+6, cc)
					cc++
				}
			}
		}
	}

	printLineAttr(left, top+row*7+6, fmt.Sprintln("24 grayscale colors"), text, background)
	for y := 0; y < 2; y++ {
		for x := 0; x < 12; x++ {
			DrawColorCell(left+numOffset+x*2, top+row*7+y+8, cc)
			cc++
		}
	}

	printLineAttr(left, top+row*7+8, fmt.Sprintf("%3d ", 232), th[1], bh[1])
	printLineAttr(left+28, top+row*7+8, fmt.Sprintf(" %-3d", 243), th[1], bh[1])
	printLineAttr(left, top+row*7+9, fmt.Sprintf("%3d ", 244), th[1], bh[1])
	printLineAttr(left+28, top+row*7+9, fmt.Sprintf(" %-3d ", 255), th[1], bh[1])
}

func Fill(x1, y1, x2, y2 int, fg, bg termbox.Attribute) {
	for y := y1; y <= y2; y++ {
		for x := x1; x <= x2; x++ {
			termbox.SetCell(x, y, ' ', fg, bg)
		}
	}
}

var (
	text, background termbox.Attribute
	lightness        int
	helpVisible      bool
)

func DrawChart() {
	screen.Clear(text, background)
	//Fill(0, 0, screen.Width-1, screen.Height-1, text, background)

	Fill(0, 0, screen.Width-1, 0, background, background+8)
	printString(2, 0, "XTerm 256 color palette chart")
	printString(screen.Width-19, 0, "F1 Help  F10 Exit")
	printLineAttr(screen.Width-19, 0, "F1", background+20, background+8)
	printLineAttr(screen.Width-10, 0, "F10", background+20, background+8)

	paletteSize := 42

	boxSize := screen.Height - 3
	w := int((screen.Width - paletteSize - 6) / 2)
	if w < boxSize {
		boxSize = w
	}

	DrawBox(2, 2, boxSize, lightness)
	// http://snag.gy/XMro8.jpg - картинка для сравнения

	helpText := []string{
		"Control keys:",
		"",
		"F1 - show/hide this help",
		"Up/Down - change lightness",
		"Left/Right - change approximation method",
		"",
		"CIE76 - fastest",
		"CIE94 - middle",
		"CIE2000 - slowest",
		"",
		"",
	}

	xc := boxSize*2 + 4
	if xc < 2 {
		xc = 2
	}
	yc := 2

	if helpVisible {
		for i, line := range helpText {
			printString(xc, yc+i, line)
		}
	} else {
		DrawPalette(xc, yc, 3, true)
	}

	screen.Flush()
}

var comparers []color.Comparer
var currentComparer int

func main() {

	err := screen.Init()
	if err != nil {
		panic(err)
	}
	defer termbox.Close()
	termbox.SetOutputMode(termbox.Output256)
	//termbox.SetInputMode(termbox.InputEsc)

	background = termNative(234)
	text = background + 12
	lightness = 100

	comparers = append(comparers, color.DeltaCIE76{}, color.DeltaCIE94{}, color.DeltaCIE2000{})

	//size := int(h/2 - 3)
	//DrawCircle(2, 1, size, 1.0)

	DrawChart()

	for loop := true; loop; {
		switch ev := termbox.PollEvent(); ev.Type {
		case termbox.EventKey:
			switch ev.Key {
			case termbox.KeyCtrlQ, termbox.KeyF10:
				loop = false
			case termbox.KeyCtrlH, termbox.KeyF1:
				helpVisible = !helpVisible
				DrawChart()
			case termbox.KeyArrowUp:
				if lightness < 100 {
					lightness = lightness + 10.0
					DrawChart()
				}
			case termbox.KeyArrowDown:
				if lightness > 0 {
					lightness = lightness - 10.0
					DrawChart()
				}
			case termbox.KeyArrowLeft:
				if currentComparer > 0 {
					currentComparer--
				} else {
					currentComparer = len(comparers) - 1
				}
				DrawChart()
			case termbox.KeyArrowRight:
				if currentComparer < 2 {
					currentComparer++
				} else {
					currentComparer = 0
				}
				DrawChart()
			}
		case termbox.EventResize:
			DrawChart()
		}
	}
}
