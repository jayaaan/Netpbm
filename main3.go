package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// PPM represents a PPM image.
type PPM struct {
	data          [][]Pixel
	width, height int
	magicNumber   string
	max           int
}

// Pixel represents a color pixel with red, green, and blue components.
type Pixel struct {
	R, G, B uint8
}

// ReadPPM reads a PPM image from a file and returns a struct that represents the image.
func ReadPPM(filename string) (*PPM, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read magic number
	scanner.Scan()
	magicNumber := scanner.Text()

	// Skip comments
	for strings.HasPrefix(scanner.Text(), "#") {
		scanner.Scan()
	}

	// Read width, height, and max
	scanner.Scan()
	dimensions := strings.Fields(scanner.Text())
	width, _ := strconv.Atoi(dimensions[0])
	height, _ := strconv.Atoi(dimensions[1])
	max, _ := strconv.Atoi(dimensions[2])

	// Read pixel data
	data := make([][]Pixel, height)
	for i := 0; i < height; i++ {
		data[i] = make([]Pixel, width)
		for j := 0; j < width; j++ {
			scanner.Scan()
			values := strings.Fields(scanner.Text())
			r, _ := strconv.Atoi(values[0])
			g, _ := strconv.Atoi(values[1])
			b, _ := strconv.Atoi(values[2])
			data[i][j] = Pixel{uint8(r), uint8(g), uint8(b)}
		}
	}

	return &PPM{
		data:        data,
		width:       width,
		height:      height,
		magicNumber: magicNumber,
		max:         max,
	}, nil
}

// Size returns the width and height of the image.
func (ppm *PPM) Size() (int, int) {
	return ppm.width, ppm.height
}

// At returns the value of the pixel at (x, y).
func (ppm *PPM) At(x, y int) Pixel {
	return ppm.data[y][x]
}

// Set sets the value of the pixel at (x, y).
func (ppm *PPM) Set(x, y int, value Pixel) {
	ppm.data[y][x] = value
}

// Save saves the PPM image to a file and returns an error if there was a problem.
func (ppm *PPM) Save(filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// Write magic number
	_, err = writer.WriteString(ppm.magicNumber + "\n")
	if err != nil {
		return err
	}

	// Write width, height, and max
	_, err = writer.WriteString(fmt.Sprintf("%d %d\n%d\n", ppm.width, ppm.height, ppm.max))
	if err != nil {
		return err
	}

	// Write pixel data
	for _, row := range ppm.data {
		for _, pixel := range row {
			_, err = writer.WriteString(fmt.Sprintf("%d %d %d ", pixel.R, pixel.G, pixel.B))
		}
		_, err = writer.WriteString("\n")
	}

	writer.Flush()

	return nil
}

// Invert inverts the colors of the PPM image.
func (ppm *PPM) Invert() {
	for y, row := range ppm.data {
		for x := range row {
			ppm.data[y][x] = Pixel{255 - row[x].R, 255 - row[x].G, 255 - row[x].B}
		}
	}
}

// Flip flips the PPM image horizontally.
func (ppm *PPM) Flip() {
	for i := 0; i < ppm.height; i++ {
		for j := 0; j < ppm.width/2; j++ {
			ppm.data[i][j], ppm.data[i][ppm.width-j-1] = ppm.data[i][ppm.width-j-1], ppm.data[i][j]
		}
	}
}

// Flop flops the PPM image vertically.
func (ppm *PPM) Flop() {
	for i, j := 0, ppm.height-1; i < j; i, j = i+1, j-1 {
		ppm.data[i], ppm.data[j] = ppm.data[j], ppm.data[i]
	}
}

// SetMagicNumber sets the magic number of the PPM image.
func (ppm *PPM) SetMagicNumber(magicNumber string) {
	ppm.magicNumber = magicNumber
}

// SetMaxValue sets the max value of the PPM image.
func (ppm *PPM) SetMaxValue(maxValue uint8) {
	ppm.max = int(maxValue)
}

// Rotate90CW rotates the PPM image 90° clockwise.
func (ppm *PPM) Rotate90CW() {
	rotatedData := make([][]Pixel, ppm.width)
	for i := 0; i < ppm.width; i++ {
		rotatedData[i] = make([]Pixel, ppm.height)
		for j := 0; j < ppm.height; j++ {
			rotatedData[i][j] = ppm.data[ppm.height-j-1][i]
		}
	}
	ppm.data = rotatedData
	ppm.width, ppm.height = ppm.height, ppm.width
}

// ToPGM converts the PPM image to PGM.
func (ppm *PPM) ToPGM() *PGM {
	pgmData := make([][]uint8, ppm.height)
	for i := 0; i < ppm.height; i++ {
		pgmData[i] = make([]uint8, ppm.width)
		for j := 0; j < ppm.width; j++ {
			// Convert RGB to grayscale using the luminosity method
			gray := uint8(0.299*float64(ppm.data[i][j].R) + 0.587*float64(ppm.data[i][j].G) + 0.114*float64(ppm.data[i][j].B))
			pgmData[i][j] = gray
		}
	}
	return &PGM{
		data:        pgmData,
		width:       ppm.width,
		height:      ppm.height,
		magicNumber: "P2",
		max:         255,
	}
}

// ToPBM converts the PPM image to PBM.
func (ppm *PPM) ToPBM() *PBM {
	pbmData := make([][]bool, ppm.height)
	threshold := uint8(ppm.max / 2) // You can adjust the threshold as needed

	for i := 0; i < ppm.height; i++ {
		pbmData[i] = make([]bool, ppm.width)
		for j := 0; j < ppm.width; j++ {
			// Convert RGB to binary using a threshold
			gray := (uint16(ppm.data[i][j].R) + uint16(ppm.data[i][j].G) + uint16(ppm.data[i][j].B)) / 3
			pbmData[i][j] = gray > uint16(threshold)
		}
	}

	return &PBM{
		data:        pbmData,
		width:       ppm.width,
		height:      ppm.height,
		magicNumber: "P1",
	}
}

// Point represents a point in the image.
type Point struct {
    X, Y int
}

// DrawLine draws a line between two points.
func (ppm *PPM) DrawLine(p1, p2 Point, color Pixel) {
    dx := p2.X - p1.X
    dy := p2.Y - p1.Y

    steps := max(abs(dx), abs(dy))

    for i := 0; i <= steps; i++ {
        x := p1.X + i*dx/steps
        y := p1.Y + i*dy/steps

        ppm.Set(x, y, color)
    }
}

// DrawRectangle draws a rectangle.
func (ppm *PPM) DrawRectangle(p1 Point, width, height int, color Pixel) {
    for y := p1.Y; y < p1.Y+height; y++ {
        for x := p1.X; x < p1.X+width; x++ {
            ppm.Set(x, y, color)
        }
    }
}

// DrawFilledRectangle draws a filled rectangle.
func (ppm *PPM) DrawFilledRectangle(p1 Point, width, height int, color Pixel) {
    for y := p1.Y; y < p1.Y+height; y++ {
        for x := p1.X; x < p1.X+width; x++ {
            ppm.Set(x, y, color)
        }
    }
}

// DrawCircle draws a circle.
func (ppm *PPM) DrawCircle(center Point, radius int, color Pixel) {
    for angle := 0; angle < 360; angle++ {
        x := int(float64(center.X) + float64(radius)*cosine(angle))
        y := int(float64(center.Y) + float64(radius)*sine(angle))

        ppm.Set(x, y, color)
    }
}

// DrawFilledCircle draws a filled circle.
func (ppm *PPM) DrawFilledCircle(center Point, radius int, color Pixel) {
    for y := center.Y - radius; y <= center.Y+radius; y++ {
        for x := center.X - radius; x <= center.X+radius; x++ {
            if (x-center.X)*(x-center.X)+(y-center.Y)*(y-center.Y) <= radius*radius {
                ppm.Set(x, y, color)
            }
        }
    }
}

// DrawTriangle draws a triangle.
func (ppm *PPM) DrawTriangle(p1, p2, p3 Point, color Pixel) {
    ppm.DrawLine(p1, p2, color)
    ppm.DrawLine(p2, p3, color)
    ppm.DrawLine(p3, p1, color)
}

// DrawFilledTriangle draws a filled triangle.
func (ppm *PPM) DrawFilledTriangle(p1, p2, p3 Point, color Pixel) {
    // Sorting vertices by y-coordinate
    vertices := []Point{p1, p2, p3}
    sort.Slice(vertices, func(i, j int) bool {
        return vertices[i].Y < vertices[j].Y
    })

    p1, p2, p3 = vertices[0], vertices[1], vertices[2]

    dx1, dx2 := float64(p2.X-p1.X)/float64(p2.Y-p1.Y), float64(p3.X-p1.X)/float64(p3.Y-p1.Y)
    x1, x2 := float64(p1.X), float64(p1.X)

    for y := p1.Y; y <= p2.Y; y++ {
        ppm.DrawLine(Point{int(x1), y}, Point{int(x2), y}, color)
        x1 += dx1
        x2 += dx2
    }

    dx1, dx2 = float64(p3.X-p2.X)/float64(p3.Y-p2.Y), float64(p3.X-p1.X)/float64(p3.Y-p1.Y)
    x1, x2 = float64(p2.X), float64(p1.X)

    for y := p2.Y; y <= p3.Y; y++ {
        ppm.DrawLine(Point{int(x1), y}, Point{int(x2), y}, color)
        x1 += dx1
        x2 += dx2
    }
}

// DrawPolygon draws a polygon.
func (ppm *PPM) DrawPolygon(points []Point, color Pixel) {
    n := len(points)
    for i := 0; i < n-1; i++ {
        ppm.DrawLine(points[i], points[i+1], color)
    }
    ppm.DrawLine(points[n-1], points[0], color)
}

// DrawFilledPolygon draws a filled polygon.
func (ppm *PPM) DrawFilledPolygon(points []Point, color Pixel) {
    n := len(points)
    maxY := points[0].Y
    minY := points[0].Y

    for i := 1; i < n; i++ {
        if points[i].Y > maxY {
            maxY = points[i].Y
        }
        if points[i].Y < minY {
            minY = points[i].Y
        }
    }

    for y := minY; y <= maxY; y++ {
        intersections := findIntersections(points, y)
        sort.Ints(intersections)

        for i := 0; i < len(intersections)-1; i += 2 {
            ppm.DrawLine(Point{intersections[i], y}, Point{intersections[i+1], y}, color)
        }
    }
}

// DrawKochSnowflake draws a Koch snowflake.
func (ppm *PPM) DrawKochSnowflake(center Point, radius int, color Pixel) {
    p1 := Point{center.X - radius, center.Y + radius}
    p2 := Point{center.X + radius, center.Y + radius}
    p3 := Point{center.X, center.Y - radius}

    ppm.DrawKochCurve(p1, p2, 4, color)
    ppm.DrawKochCurve(p2, p3, 4, color)
    ppm.DrawKochCurve(p3, p1, 4, color)
}

// DrawKochCurve recursively draws a Koch curve.
func (ppm *PPM) DrawKochCurve(p1, p2 Point, depth int, color Pixel) {
    if depth == 0 {
        ppm.DrawLine(p1, p2, color)
    } else {
        // Calculate intermediate points
        p3 := Point{
            X: (2*p1.X + p2.X) / 3,
            Y: (2*p1.Y + p2.Y) / 3,
        }

        p4 := Point{
            X: (p1.X + 2*p2.X) / 3,
            Y: (p1.Y + 2*p2.Y) / 3,
        }

        angle := math.Pi / 3.0
        dx := p4.X - p3.X
        dy := p4.Y - p3.Y

        // Calculate the third point of an equilateral triangle
        p5 := Point{
            X: p3.X + int(math.Cos(angle)*dx-math) 
		}
	} 
} 
