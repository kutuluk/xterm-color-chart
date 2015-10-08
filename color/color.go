package color

import (
	"fmt"
	"math"
)

const (
	epsilon float64 = 216.0 / 24389.0
	kappa   float64 = 24389.0 / 27.0
)

var (
	// reference white in XYZ coordinates
	// цветовые координаты излучателя белого света в пространстве XYZ
	whitePoints = map[string][3]float64{
		"A":   {1.09850, 1.0, 0.35585},
		"B":   {0.99072, 1.0, 0.85223},
		"C":   {0.98074, 1.0, 1.18232},
		"D50": {0.96422, 1.0, 0.82521},
		"D55": {0.95682, 1.0, 0.92149},
		"D65": {0.95047, 1.0, 1.08883},
		"D75": {0.94972, 1.0, 1.22638},
		"E":   {1.00000, 1.0, 1.00000},
		"F2":  {0.99186, 1.0, 0.67393},
		"F7":  {0.95041, 1.0, 1.08747},
		"F11": {1.00962, 1.0, 0.64350},
	}

	// sRGB to XYZ conversion matrix
	convertMatrix = [3][3]float64{
		{0.4124564, 0.3575761, 0.1804375},
		{0.2126729, 0.7151522, 0.0721750},
		{0.0193339, 0.1191920, 0.9503041},
	}

	// XYZ to sRGB conversion matrix
	inverseMatrix = [3][3]float64{
		{3.2404542, -1.5371385, -0.4985314},
		{-0.9692660, 1.8760108, 0.0415560},
		{0.0556434, -0.2040259, 1.0572252},
	}
)

func sqr(v float64) float64 {
	return v * v
}

type RGBColor struct {
	R, G, B int // RGB-компоненты [0..255]
}

func (rgb RGBColor) String() string {
	return fmt.Sprintf("RGB:{%d %d %d}", rgb.R, rgb.G, rgb.B)
}

func (rgb RGBColor) RGB() RGBColor {
	return rgb
}

func (rgb RGBColor) XYZ() XYZColor {

	f := func(n float64) float64 {
		if n > 0.04045 {
			// math.Pow(base, exponent) = math.Exp(math.Log(base) * exponent)
			//return math.Exp(2.4 * math.Log((n+0.055)/1.055))
			return math.Pow(((n + 0.055) / 1.055), 2.4)
		} else {
			return n / 12.92
		}
	}

	r := f(float64(rgb.R) / 255.0)
	g := f(float64(rgb.G) / 255.0)
	b := f(float64(rgb.B) / 255.0)

	//
	// [X]   [Xr Xg Xb]   [R]
	// [Y] = [Yr Yg Yb] x [G]
	// [Z]   [Zr Zg Zb]   [B]
	//
	//Observer. = 2°, Illuminant = D65
	return XYZColor{
		r*0.4124 + g*0.3576 + b*0.1805,
		r*0.2126 + g*0.7152 + b*0.0722,
		r*0.0193 + g*0.1192 + b*0.9505,
	}
}

func (rgb RGBColor) Lab() LabColor {
	return rgb.XYZ().Lab()
}

type XYZColor struct {
	X, Y, Z float64
}

func (xyz XYZColor) String() string {
	return fmt.Sprintf("XYZ:{%4.7f %4.7f %4.7f}", xyz.X, xyz.Y, xyz.Z)
}

func (xyz XYZColor) RGB() RGBColor {
	// TODO!!!
	return RGBColor{0, 0, 0}
}

func (xyz XYZColor) XYZ() XYZColor {
	return xyz
}

// ok
func (xyz XYZColor) Lab() LabColor {

	f := func(n float64) float64 {
		if n > epsilon {
			return math.Pow(n, 1.0/3.0)
		} else {
			return (kappa*n + 16.0) / 116.0
		}
	}

	//  Observer= 2°, Illuminant= D65
	whitePoint := whitePoints["D65"]

	x := f(xyz.X / whitePoint[0])
	y := f(xyz.Y / whitePoint[1])
	z := f(xyz.Z / whitePoint[2])

	return LabColor{
		(116.0 * y) - 16.0,
		500.0 * (x - y),
		200.0 * (y - z),
	}
}

type LabColor struct {
	L, A, B float64
}

func (lab LabColor) String() string {
	return fmt.Sprintf("Lab:{%4.7f %4.7f %4.7f}", lab.L, lab.A, lab.B)
}

func (lab LabColor) RGB() RGBColor {
	return lab.XYZ().RGB()
}

// ok
func (lab LabColor) XYZ() XYZColor {

	f := func(n float64) float64 {
		if n3 := n * n * n; n3 > epsilon {
			return n3
		} else {
			return (n*116.0 - 16.0) / kappa
		}
	}

	y := (lab.L + 16.0) / 116.0
	x := 0.002*lab.A + y
	z := y - 0.005*lab.B

	//  Observer= 2°, Illuminant= D65
	whitePoint := whitePoints["D65"]

	return XYZColor{
		f(x) * whitePoint[0],
		f(y) * whitePoint[1],
		f(z) * whitePoint[2],
	}
}

func (lab LabColor) Lab() LabColor {
	return lab
}

// ok
func (lab LabColor) LCh() LChColor {
	c := math.Sqrt(sqr(lab.A) + sqr(lab.B))
	h := 180.0 * math.Atan2(lab.B, lab.A) / math.Pi
	if h < 0.0 {
		h += 360.0
	}
	return LChColor{lab.L, c, h}
}

type LChColor struct {
	L, C, H float64
}

func (lch LChColor) RGB() RGBColor {
	return lch.Lab().XYZ().RGB()
}

func (lch LChColor) XYZ() XYZColor {
	return lch.Lab().XYZ()
}

// ok
func (lch LChColor) Lab() LabColor {
	a := lch.C * math.Cos(lch.H*math.Pi/180.0)
	b := lch.C * math.Sin(lch.H*math.Pi/180.0)
	return LabColor{lch.L, a, b}
}

type Colorer interface {
	RGB() RGBColor
	XYZ() XYZColor
	Lab() LabColor
}

type Comparer interface {
	Compare(c1, c2 Colorer) float64
}

type DeltaCIE76 struct{}

// DeltaCIE76.Compare computes the CIE76 color difference.
// This is just Euclidean distance in Lab space, and therefore quite fast,
// though it exhibits perceptual uniformity issues especially in the blue and desaturated regions.
func (DeltaCIE76) Compare(c1, c2 Colorer) float64 {
	lab1 := c1.Lab()
	lab2 := c2.Lab()
	return math.Sqrt(sqr(lab1.L-lab2.L) + sqr(lab1.A-lab2.A) + sqr(lab1.B-lab2.B))
}

type DeltaCIE94 struct{}

// KLCh94 is a struct for weighting factors for CIE94 ΔE calculation.
type KLCh94 struct {
	KL, KC, Kh, K1, K2 float64
}

// KLCH94GraphicArts are the weighting factors for CIE94 used for most uses except textiles.
var KLCH94GraphicArts = KLCh94{1, 1, 1, 0.045, 0.015}

// KLCH94Textiles are the weighting factors for CIE94 used for textiles.
var KLCH94Textiles = KLCh94{2, 1, 1, 0.048, 0.014}

//var KLCH94 = &KLCH94GraphicArts

// DeltaCIE94.Compare computes the CIE94 color difference of two L*a*b* colors.
// This is a distance calculation with the addition of weighting factors specified by klch.
func (DeltaCIE94) Compare(c1, c2 Colorer) float64 {
	lab1 := c1.Lab()
	lab2 := c2.Lab()

	dL := sqr(lab1.L - lab2.L)

	c1ab := math.Sqrt(sqr(lab1.A) + sqr(lab1.B))
	c2ab := math.Sqrt(sqr(lab2.A) + sqr(lab2.B))

	dC := sqr(c1ab - c2ab)

	dH := sqr(lab1.A-lab2.A) + sqr(lab1.B-lab2.B) - dC

	KLCH94 := &KLCH94GraphicArts

	sC := 1.0 + KLCH94.K1*c1ab
	sH := 1.0 + KLCH94.K2*c1ab

	return math.Sqrt(
		dL/sqr(KLCH94.KL) +
			dC/sqr(KLCH94.KC*sC) +
			dH/sqr(KLCH94.Kh*sH))
}

type DeltaCIE2000 struct{}

type KLCh struct {
	KL, KC, Kh float64
}

// KLCHDefault is the most commonly used set of weighting parameters for CIEDE2000
var KLCH = KLCh{1, 1, 1}

// CIE2000 computes the CIEDE2000 delta-E for two L*a*b* space color coordinates
// klch is for configuring the weighting factors, but this almost always should be KLCHDefault
// Note that this implementation will exhibit slightly different behavior around the discontinuities
// of the function (these are grey colors) compared to Java and most C runtimes. The golang atan
// function has different accuracy characteristics compared to most Unix platforms and Java Strict math
func (DeltaCIE2000) Compare(col1, col2 Colorer) float64 {
	lab1 := col1.Lab()
	lab2 := col2.Lab()

	lBarPrime := (lab1.L + lab2.L) * 0.5
	c1 := math.Sqrt(lab1.A*lab1.A + lab1.B*lab1.B)
	c2 := math.Sqrt(lab2.A*lab2.A + lab2.B*lab2.B)
	cBar := (c1 + c2) * 0.5

	cBar7 := cBar * cBar * cBar
	cBar7 *= cBar7 * cBar
	g := 0.5 * (1.0 - math.Sqrt(cBar7/(cBar7+6103515625.0))) // 25**7

	a1Prime := (1.0 + g) * lab1.A
	a2Prime := (1.0 + g) * lab2.A

	c1Prime := math.Sqrt(a1Prime*a1Prime + lab1.B*lab1.B)
	c2Prime := math.Sqrt(a2Prime*a2Prime + lab2.B*lab2.B)

	cBarPrime := (c1Prime + c2Prime) * 0.5

	h1Prime := math.Atan2(lab1.B, a1Prime)
	if h1Prime < 0 {
		h1Prime += 2 * math.Pi

	}
	h2Prime := math.Atan2(lab2.B, a2Prime)
	if h2Prime < 0 {
		h2Prime += 2 * math.Pi

	}

	hBarPrime := (h1Prime + h2Prime) * 0.5
	dhPrime := h2Prime - h1Prime
	if math.Abs(dhPrime) > math.Pi {
		hBarPrime += math.Pi
		if h2Prime <= h1Prime {
			dhPrime += 2 * math.Pi

		} else {
			dhPrime -= 2 * math.Pi

		}

	}

	t := 1.0 -
		0.17*math.Cos(hBarPrime-math.Pi/6) +
		0.24*math.Cos(2.0*hBarPrime) +
		0.32*math.Cos(3.0*hBarPrime+math.Pi/30) -
		0.20*math.Cos(4.0*hBarPrime-63.0*math.Pi/180)

	dLPrime := lab2.L - lab1.L
	dCPrime := c2Prime - c1Prime
	dHPrime := 2.0 * math.Sqrt(c1Prime*c2Prime) * math.Sin(dhPrime/2.0)

	lBarPrimeM50Sqr := lBarPrime - 50.0
	lBarPrimeM50Sqr *= lBarPrimeM50Sqr
	sL := 1.0 + (0.015*lBarPrimeM50Sqr)/math.Sqrt(20.0+lBarPrimeM50Sqr)
	sC := 1.0 + 0.045*cBarPrime
	sH := 1.0 + 0.015*cBarPrime*t

	hBarPrimeM := (180/math.Pi*hBarPrime - 275.0) / 25.0
	dTheta := math.Pi / 6 * math.Exp(-hBarPrimeM*hBarPrimeM)
	cBarPrime7 := cBarPrime * cBarPrime * cBarPrime
	cBarPrime7 *= cBarPrime7 * cBarPrime
	rC := math.Sqrt(cBarPrime7 / (cBarPrime7 + 6103515625.0))
	rT := -2.0 * rC * math.Sin(2.0*dTheta)

	return math.Sqrt(
		sqr(dLPrime/(KLCH.KL*sL)) +
			sqr(dCPrime/(KLCH.KC*sC)) +
			sqr(dHPrime/(KLCH.Kh*sH)) +
			(dCPrime/(KLCH.KC*sC))*(dHPrime/(KLCH.Kh*sH))*rT)

}

func Approximate(color Colorer, comparer Comparer) (float64, int) {
	best := 0
	bestdist := 10000000.0
	for i, applicant := range XtermLabPalette {
		if dist := comparer.Compare(color, applicant); dist < bestdist {
			best, bestdist = i, dist
		}
	}
	return bestdist, best + 17
}

//color = round(36 * (r * 5) + 6 * (g * 5) + (b * 5) + 16)
func HackApproximate(color Colorer) int {
	c := color.RGB()
	r := float64(c.R / 255)
	g := float64(c.B / 255)
	b := float64(c.G / 255)
	return int(36*(r*5)+6*(g*5)+b*5) + 17
}

var (
	/*
	   {
	   		{R: 0x00, G: 0x00, B: 0x00}, // "#000000" 16 index of palette
	   		{R: 0x00, G: 0x00, B: 0x5F}, // "#00005f"
	   		{R: 0x00, G: 0x00, B: 0x87}, // "#000087"
	   		{R: 0x00, G: 0x00, B: 0xAF}, // "#0000af"
	   		{R: 0x00, G: 0x00, B: 0xD7}, // "#0000d7"
	   		{R: 0x00, G: 0x00, B: 0xFF}, // "#0000ff"
	   		{R: 0x00, G: 0x5F, B: 0x00}, // "#005f00"
	   		{R: 0x00, G: 0x5F, B: 0x5F}, // "#005f5f"
	   		{R: 0x00, G: 0x5F, B: 0x87}, // "#005f87"
	   		{R: 0x00, G: 0x5F, B: 0xAF}, // "#005faf"
	   		{R: 0x00, G: 0x5F, B: 0xD7}, // "#005fd7"
	   		{R: 0x00, G: 0x5F, B: 0xFF}, // "#005fff"
	   		{R: 0x00, G: 0x87, B: 0x00}, // "#008700"
	   		{R: 0x00, G: 0x87, B: 0x5F}, // "#00875f"
	   		{R: 0x00, G: 0x87, B: 0x87}, // "#008787"
	   		{R: 0x00, G: 0x87, B: 0xAF}, // "#0087af"
	   		{R: 0x00, G: 0x87, B: 0xD7}, // "#0087d7"
	   		{R: 0x00, G: 0x87, B: 0xFF}, // "#0087ff"
	   		{R: 0x00, G: 0xAF, B: 0x00}, // "#00af00"
	   		{R: 0x00, G: 0xAF, B: 0x5F}, // "#00af5f"
	   		{R: 0x00, G: 0xAF, B: 0x87}, // "#00af87"
	   		{R: 0x00, G: 0xAF, B: 0xAF}, // "#00afaf"
	   		{R: 0x00, G: 0xAF, B: 0xD7}, // "#00afd7"
	   		{R: 0x00, G: 0xAF, B: 0xFF}, // "#00afff"
	   		{R: 0x00, G: 0xD7, B: 0x00}, // "#00d700"
	   		{R: 0x00, G: 0xD7, B: 0x5F}, // "#00d75f"
	   		{R: 0x00, G: 0xD7, B: 0x87}, // "#00d787"
	   		{R: 0x00, G: 0xD7, B: 0xAF}, // "#00d7af"
	   		{R: 0x00, G: 0xD7, B: 0xD7}, // "#00d7d7"
	   		{R: 0x00, G: 0xD7, B: 0xFF}, // "#00d7ff"
	   		{R: 0x00, G: 0xFF, B: 0x00}, // "#00ff00"
	   		{R: 0x00, G: 0xFF, B: 0x5F}, // "#00ff5f"
	   		{R: 0x00, G: 0xFF, B: 0x87}, // "#00ff87"
	   		{R: 0x00, G: 0xFF, B: 0xAF}, // "#00ffaf"
	   		{R: 0x00, G: 0xFF, B: 0xD7}, // "#00ffd7"
	   		{R: 0x00, G: 0xFF, B: 0xFF}, // "#00ffff"
	   		{R: 0x5F, G: 0x00, B: 0x00}, // "#5f0000"
	   		{R: 0x5F, G: 0x00, B: 0x5F}, // "#5f005f"
	   		{R: 0x5F, G: 0x00, B: 0x87}, // "#5f0087"
	   		{R: 0x5F, G: 0x00, B: 0xAF}, // "#5f00af"
	   		{R: 0x5F, G: 0x00, B: 0xD7}, // "#5f00d7"
	   		{R: 0x5F, G: 0x00, B: 0xFF}, // "#5f00ff"
	   		{R: 0x5F, G: 0x5F, B: 0x00}, // "#5f5f00"
	   		{R: 0x5F, G: 0x5F, B: 0x5F}, // "#5f5f5f"
	   		{R: 0x5F, G: 0x5F, B: 0x87}, // "#5f5f87"
	   		{R: 0x5F, G: 0x5F, B: 0xAF}, // "#5f5faf"
	   		{R: 0x5F, G: 0x5F, B: 0xD7}, // "#5f5fd7"
	   		{R: 0x5F, G: 0x5F, B: 0xFF}, // "#5f5fff"
	   		{R: 0x5F, G: 0x87, B: 0x00}, // "#5f8700"
	   		{R: 0x5F, G: 0x87, B: 0x5F}, // "#5f875f"
	   		{R: 0x5F, G: 0x87, B: 0x87}, // "#5f8787"
	   		{R: 0x5F, G: 0x87, B: 0xAF}, // "#5f87af"
	   		{R: 0x5F, G: 0x87, B: 0xD7}, // "#5f87d7"
	   		{R: 0x5F, G: 0x87, B: 0xFF}, // "#5f87ff"
	   		{R: 0x5F, G: 0xAF, B: 0x00}, // "#5faf00"
	   		{R: 0x5F, G: 0xAF, B: 0x5F}, // "#5faf5f"
	   		{R: 0x5F, G: 0xAF, B: 0x87}, // "#5faf87"
	   		{R: 0x5F, G: 0xAF, B: 0xAF}, // "#5fafaf"
	   		{R: 0x5F, G: 0xAF, B: 0xD7}, // "#5fafd7"
	   		{R: 0x5F, G: 0xAF, B: 0xFF}, // "#5fafff"
	   		{R: 0x5F, G: 0xD7, B: 0x00}, // "#5fd700"
	   		{R: 0x5F, G: 0xD7, B: 0x5F}, // "#5fd75f"
	   		{R: 0x5F, G: 0xD7, B: 0x87}, // "#5fd787"
	   		{R: 0x5F, G: 0xD7, B: 0xAF}, // "#5fd7af"
	   		{R: 0x5F, G: 0xD7, B: 0xD7}, // "#5fd7d7"
	   		{R: 0x5F, G: 0xD7, B: 0xFF}, // "#5fd7ff"
	   		{R: 0x5F, G: 0xFF, B: 0x00}, // "#5fff00"
	   		{R: 0x5F, G: 0xFF, B: 0x5F}, // "#5fff5f"
	   		{R: 0x5F, G: 0xFF, B: 0x87}, // "#5fff87"
	   		{R: 0x5F, G: 0xFF, B: 0xAF}, // "#5fffaf"
	   		{R: 0x5F, G: 0xFF, B: 0xD7}, // "#5fffd7"
	   		{R: 0x5F, G: 0xFF, B: 0xFF}, // "#5fffff"
	   		{R: 0x87, G: 0x00, B: 0x00}, // "#870000"
	   		{R: 0x87, G: 0x00, B: 0x5F}, // "#87005f"
	   		{R: 0x87, G: 0x00, B: 0x87}, // "#870087"
	   		{R: 0x87, G: 0x00, B: 0xAF}, // "#8700af"
	   		{R: 0x87, G: 0x00, B: 0xD7}, // "#8700d7"
	   		{R: 0x87, G: 0x00, B: 0xFF}, // "#8700ff"
	   		{R: 0x87, G: 0x5F, B: 0x00}, // "#875f00"
	   		{R: 0x87, G: 0x5F, B: 0x5F}, // "#875f5f"
	   		{R: 0x87, G: 0x5F, B: 0x87}, // "#875f87"
	   		{R: 0x87, G: 0x5F, B: 0xAF}, // "#875faf"
	   		{R: 0x87, G: 0x5F, B: 0xD7}, // "#875fd7"
	   		{R: 0x87, G: 0x5F, B: 0xFF}, // "#875fff"
	   		{R: 0x87, G: 0x87, B: 0x00}, // "#878700"
	   		{R: 0x87, G: 0x87, B: 0x5F}, // "#87875f"
	   		{R: 0x87, G: 0x87, B: 0x87}, // "#878787"
	   		{R: 0x87, G: 0x87, B: 0xAF}, // "#8787af"
	   		{R: 0x87, G: 0x87, B: 0xD7}, // "#8787d7"
	   		{R: 0x87, G: 0x87, B: 0xFF}, // "#8787ff"
	   		{R: 0x87, G: 0xAF, B: 0x00}, // "#87af00"
	   		{R: 0x87, G: 0xAF, B: 0x5F}, // "#87af5f"
	   		{R: 0x87, G: 0xAF, B: 0x87}, // "#87af87"
	   		{R: 0x87, G: 0xAF, B: 0xAF}, // "#87afaf"
	   		{R: 0x87, G: 0xAF, B: 0xD7}, // "#87afd7"
	   		{R: 0x87, G: 0xAF, B: 0xFF}, // "#87afff"
	   		{R: 0x87, G: 0xD7, B: 0x00}, // "#87d700"
	   		{R: 0x87, G: 0xD7, B: 0x5F}, // "#87d75f"
	   		{R: 0x87, G: 0xD7, B: 0x87}, // "#87d787"
	   		{R: 0x87, G: 0xD7, B: 0xAF}, // "#87d7af"
	   		{R: 0x87, G: 0xD7, B: 0xD7}, // "#87d7d7"
	   		{R: 0x87, G: 0xD7, B: 0xFF}, // "#87d7ff"
	   		{R: 0x87, G: 0xFF, B: 0x00}, // "#87ff00"
	   		{R: 0x87, G: 0xFF, B: 0x5F}, // "#87ff5f"
	   		{R: 0x87, G: 0xFF, B: 0x87}, // "#87ff87"
	   		{R: 0x87, G: 0xFF, B: 0xAF}, // "#87ffaf"
	   		{R: 0x87, G: 0xFF, B: 0xD7}, // "#87ffd7"
	   		{R: 0x87, G: 0xFF, B: 0xFF}, // "#87ffff"
	   		{R: 0xAF, G: 0x00, B: 0x00}, // "#af0000"
	   		{R: 0xAF, G: 0x00, B: 0x5F}, // "#af005f"
	   		{R: 0xAF, G: 0x00, B: 0x87}, // "#af0087"
	   		{R: 0xAF, G: 0x00, B: 0xAF}, // "#af00af"
	   		{R: 0xAF, G: 0x00, B: 0xD7}, // "#af00d7"
	   		{R: 0xAF, G: 0x00, B: 0xFF}, // "#af00ff"
	   		{R: 0xAF, G: 0x5F, B: 0x00}, // "#af5f00"
	   		{R: 0xAF, G: 0x5F, B: 0x5F}, // "#af5f5f"
	   		{R: 0xAF, G: 0x5F, B: 0x87}, // "#af5f87"
	   		{R: 0xAF, G: 0x5F, B: 0xAF}, // "#af5faf"
	   		{R: 0xAF, G: 0x5F, B: 0xD7}, // "#af5fd7"
	   		{R: 0xAF, G: 0x5F, B: 0xFF}, // "#af5fff"
	   		{R: 0xAF, G: 0x87, B: 0x00}, // "#af8700"
	   		{R: 0xAF, G: 0x87, B: 0x5F}, // "#af875f"
	   		{R: 0xAF, G: 0x87, B: 0x87}, // "#af8787"
	   		{R: 0xAF, G: 0x87, B: 0xAF}, // "#af87af"
	   		{R: 0xAF, G: 0x87, B: 0xD7}, // "#af87d7"
	   		{R: 0xAF, G: 0x87, B: 0xFF}, // "#af87ff"
	   		{R: 0xAF, G: 0xAF, B: 0x00}, // "#afaf00"
	   		{R: 0xAF, G: 0xAF, B: 0x5F}, // "#afaf5f"
	   		{R: 0xAF, G: 0xAF, B: 0x87}, // "#afaf87"
	   		{R: 0xAF, G: 0xAF, B: 0xAF}, // "#afafaf"
	   		{R: 0xAF, G: 0xAF, B: 0xD7}, // "#afafd7"
	   		{R: 0xAF, G: 0xAF, B: 0xFF}, // "#afafff"
	   		{R: 0xAF, G: 0xD7, B: 0x00}, // "#afd700"
	   		{R: 0xAF, G: 0xD7, B: 0x5F}, // "#afd75f"
	   		{R: 0xAF, G: 0xD7, B: 0x87}, // "#afd787"
	   		{R: 0xAF, G: 0xD7, B: 0xAF}, // "#afd7af"
	   		{R: 0xAF, G: 0xD7, B: 0xD7}, // "#afd7d7"
	   		{R: 0xAF, G: 0xD7, B: 0xFF}, // "#afd7ff"
	   		{R: 0xAF, G: 0xFF, B: 0x00}, // "#afff00"
	   		{R: 0xAF, G: 0xFF, B: 0x5F}, // "#afff5f"
	   		{R: 0xAF, G: 0xFF, B: 0x87}, // "#afff87"
	   		{R: 0xAF, G: 0xFF, B: 0xAF}, // "#afffaf"
	   		{R: 0xAF, G: 0xFF, B: 0xD7}, // "#afffd7"
	   		{R: 0xAF, G: 0xFF, B: 0xFF}, // "#afffff"
	   		{R: 0xD7, G: 0x00, B: 0x00}, // "#d70000"
	   		{R: 0xD7, G: 0x00, B: 0x5F}, // "#d7005f"
	   		{R: 0xD7, G: 0x00, B: 0x87}, // "#d70087"
	   		{R: 0xD7, G: 0x00, B: 0xAF}, // "#d700af"
	   		{R: 0xD7, G: 0x00, B: 0xD7}, // "#d700d7"
	   		{R: 0xD7, G: 0x00, B: 0xFF}, // "#d700ff"
	   		{R: 0xD7, G: 0x5F, B: 0x00}, // "#d75f00"
	   		{R: 0xD7, G: 0x5F, B: 0x5F}, // "#d75f5f"
	   		{R: 0xD7, G: 0x5F, B: 0x87}, // "#d75f87"
	   		{R: 0xD7, G: 0x5F, B: 0xAF}, // "#d75faf"
	   		{R: 0xD7, G: 0x5F, B: 0xD7}, // "#d75fd7"
	   		{R: 0xD7, G: 0x5F, B: 0xFF}, // "#d75fff"
	   		{R: 0xD7, G: 0x87, B: 0x00}, // "#d78700"
	   		{R: 0xD7, G: 0x87, B: 0x5F}, // "#d7875f"
	   		{R: 0xD7, G: 0x87, B: 0x87}, // "#d78787"
	   		{R: 0xD7, G: 0x87, B: 0xAF}, // "#d787af"
	   		{R: 0xD7, G: 0x87, B: 0xD7}, // "#d787d7"
	   		{R: 0xD7, G: 0x87, B: 0xFF}, // "#d787ff"
	   		{R: 0xD7, G: 0xAF, B: 0x00}, // "#d7af00"
	   		{R: 0xD7, G: 0xAF, B: 0x5F}, // "#d7af5f"
	   		{R: 0xD7, G: 0xAF, B: 0x87}, // "#d7af87"
	   		{R: 0xD7, G: 0xAF, B: 0xAF}, // "#d7afaf"
	   		{R: 0xD7, G: 0xAF, B: 0xD7}, // "#d7afd7"
	   		{R: 0xD7, G: 0xAF, B: 0xFF}, // "#d7afff"
	   		{R: 0xD7, G: 0xD7, B: 0x00}, // "#d7d700"
	   		{R: 0xD7, G: 0xD7, B: 0x5F}, // "#d7d75f"
	   		{R: 0xD7, G: 0xD7, B: 0x87}, // "#d7d787"
	   		{R: 0xD7, G: 0xD7, B: 0xAF}, // "#d7d7af"
	   		{R: 0xD7, G: 0xD7, B: 0xD7}, // "#d7d7d7"
	   		{R: 0xD7, G: 0xD7, B: 0xFF}, // "#d7d7ff"
	   		{R: 0xD7, G: 0xFF, B: 0x00}, // "#d7ff00"
	   		{R: 0xD7, G: 0xFF, B: 0x5F}, // "#d7ff5f"
	   		{R: 0xD7, G: 0xFF, B: 0x87}, // "#d7ff87"
	   		{R: 0xD7, G: 0xFF, B: 0xAF}, // "#d7ffaf"
	   		{R: 0xD7, G: 0xFF, B: 0xD7}, // "#d7ffd7"
	   		{R: 0xD7, G: 0xFF, B: 0xFF}, // "#d7ffff"
	   		{R: 0xFF, G: 0x00, B: 0x00}, // "#ff0000"
	   		{R: 0xFF, G: 0x00, B: 0x5F}, // "#ff005f"
	   		{R: 0xFF, G: 0x00, B: 0x87}, // "#ff0087"
	   		{R: 0xFF, G: 0x00, B: 0xAF}, // "#ff00af"
	   		{R: 0xFF, G: 0x00, B: 0xD7}, // "#ff00d7"
	   		{R: 0xFF, G: 0x00, B: 0xFF}, // "#ff00ff"
	   		{R: 0xFF, G: 0x5F, B: 0x00}, // "#ff5f00"
	   		{R: 0xFF, G: 0x5F, B: 0x5F}, // "#ff5f5f"
	   		{R: 0xFF, G: 0x5F, B: 0x87}, // "#ff5f87"
	   		{R: 0xFF, G: 0x5F, B: 0xAF}, // "#ff5faf"
	   		{R: 0xFF, G: 0x5F, B: 0xD7}, // "#ff5fd7"
	   		{R: 0xFF, G: 0x5F, B: 0xFF}, // "#ff5fff"
	   		{R: 0xFF, G: 0x87, B: 0x00}, // "#ff8700"
	   		{R: 0xFF, G: 0x87, B: 0x5F}, // "#ff875f"
	   		{R: 0xFF, G: 0x87, B: 0x87}, // "#ff8787"
	   		{R: 0xFF, G: 0x87, B: 0xAF}, // "#ff87af"
	   		{R: 0xFF, G: 0x87, B: 0xD7}, // "#ff87d7"
	   		{R: 0xFF, G: 0x87, B: 0xFF}, // "#ff87ff"
	   		{R: 0xFF, G: 0xAF, B: 0x00}, // "#ffaf00"
	   		{R: 0xFF, G: 0xAF, B: 0x5F}, // "#ffaf5f"
	   		{R: 0xFF, G: 0xAF, B: 0x87}, // "#ffaf87"
	   		{R: 0xFF, G: 0xAF, B: 0xAF}, // "#ffafaf"
	   		{R: 0xFF, G: 0xAF, B: 0xD7}, // "#ffafd7"
	   		{R: 0xFF, G: 0xAF, B: 0xFF}, // "#ffafff"
	   		{R: 0xFF, G: 0xD7, B: 0x00}, // "#ffd700"
	   		{R: 0xFF, G: 0xD7, B: 0x5F}, // "#ffd75f"
	   		{R: 0xFF, G: 0xD7, B: 0x87}, // "#ffd787"
	   		{R: 0xFF, G: 0xD7, B: 0xAF}, // "#ffd7af"
	   		{R: 0xFF, G: 0xD7, B: 0xD7}, // "#ffd7d7"
	   		{R: 0xFF, G: 0xD7, B: 0xFF}, // "#ffd7ff"
	   		{R: 0xFF, G: 0xFF, B: 0x00}, // "#ffff00"
	   		{R: 0xFF, G: 0xFF, B: 0x5F}, // "#ffff5f"
	   		{R: 0xFF, G: 0xFF, B: 0x87}, // "#ffff87"
	   		{R: 0xFF, G: 0xFF, B: 0xAF}, // "#ffffaf"
	   		{R: 0xFF, G: 0xFF, B: 0xD7}, // "#ffffd7"
	   		{R: 0xFF, G: 0xFF, B: 0xFF}, // "#ffffff"
	   		{R: 0x08, G: 0x08, B: 0x08}, // "#080808"
	   		{R: 0x12, G: 0x12, B: 0x12}, // "#121212"
	   		{R: 0x1C, G: 0x1C, B: 0x1C}, // "#1c1c1c"
	   		{R: 0x26, G: 0x26, B: 0x26}, // "#262626"
	   		{R: 0x30, G: 0x30, B: 0x30}, // "#303030"
	   		{R: 0x3A, G: 0x3A, B: 0x3A}, // "#3a3a3a"
	   		{R: 0x44, G: 0x44, B: 0x44}, // "#444444"
	   		{R: 0x4E, G: 0x4E, B: 0x4E}, // "#4e4e4e"
	   		{R: 0x58, G: 0x58, B: 0x58}, // "#585858"
	   		{R: 0x60, G: 0x60, B: 0x60}, // "#606060"
	   		{R: 0x66, G: 0x66, B: 0x66}, // "#666666"
	   		{R: 0x76, G: 0x76, B: 0x76}, // "#767676"
	   		{R: 0x80, G: 0x80, B: 0x80}, // "#808080"
	   		{R: 0x8A, G: 0x8A, B: 0x8A}, // "#8a8a8a"
	   		{R: 0x94, G: 0x94, B: 0x94}, // "#949494"
	   		{R: 0x9E, G: 0x9E, B: 0x9E}, // "#9e9e9e"
	   		{R: 0xA8, G: 0xA8, B: 0xA8}, // "#a8a8a8"
	   		{R: 0xB2, G: 0xB2, B: 0xB2}, // "#b2b2b2"
	   		{R: 0xBC, G: 0xBC, B: 0xBC}, // "#bcbcbc"
	   		{R: 0xC6, G: 0xC6, B: 0xC6}, // "#c6c6c6"
	   		{R: 0xD0, G: 0xD0, B: 0xD0}, // "#d0d0d0"
	   		{R: 0xDA, G: 0xDA, B: 0xDA}, // "#dadada"
	   		{R: 0xE4, G: 0xE4, B: 0xE4}, // "#e4e4e4"
	   		{R: 0xEE, G: 0xEE, B: 0xEE}, // "#eeeeee"
	   	}
	*/
	XtermRGBPalette [240]RGBColor
	XtermLabPalette [240]LabColor
)

func init() {
	// calculate xterm 240 color RGB palette

	// calculate 6x6x6 color cube
	cubeLevels := [6]int{0x00, 0x5f, 0x87, 0xaf, 0xd7, 0xff}
	for r := 0; r < 6; r++ {
		for g := 0; g < 6; g++ {
			for b := 0; b < 6; b++ {
				XtermRGBPalette[r*36+g*6+b] = RGBColor{cubeLevels[r], cubeLevels[g], cubeLevels[b]}
			}
		}
	}
	// calculate grayscale ramp
	for i := 0; i < 24; i++ {
		XtermRGBPalette[i+216] = RGBColor{0x08 + i*0xA, 0x08 + i*0xA, 0x08 + i*0xA}
	}

	// calculate xterm 240 color Lab palette
	for i, color := range XtermRGBPalette {
		XtermLabPalette[i] = color.Lab()
	}
}
