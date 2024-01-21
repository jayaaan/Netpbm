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

func (ppm *PPM) DrawKochSnowflake(n int, start Point, width int, color Pixel) {
	// Draw the initial triangle
	ppm.drawKochCurve(start, Point{start.X + width, start.Y}, n, color)
	ppm.drawKochCurve(Point{start.X + width, start.Y}, Point{start.X + width/2, start.Y - int(math.Sqrt(3)*float64(width/2))}, n, color)
	ppm.drawKochCurve(Point{start.X + width/2, start.Y - int(math.Sqrt(3)*float64(width/2))}, start, n, color)
}

// Helper function to draw a Koch curve
func (ppm *PPM) drawKochCurve(p1, p2 Point, n int, color Pixel) {
	if n == 0 {
		// Draw a line segment
		ppm.drawLine(p1, p2, color)
		return
	}

	// Calculate the intermediate points
	p3 := calculateIntermediatePoint(p1, p2, 1.0/3.0)
	p4 := calculateIntermediatePoint(p1, p2, 2.0/3.0)
	p5 := calculateKochVertex(p1, p2)

	// Recursively draw the four sub-segments
	ppm.drawKochCurve(p1, p3, n-1, color)
	ppm.drawKochCurve(p3, p5, n-1, color)
	ppm.drawKochCurve(p5, p4, n-1, color)
	ppm.drawKochCurve(p4, p2, n-1, color)
}

// Helper function to draw a line segment
func (ppm *PPM) drawLine(p1, p2 Point, color Pixel) {
	// You can implement the line drawing logic here
	// Update the image or perform other necessary actions
}

// Helper function to calculate the intermediate point between two points
func calculateIntermediatePoint(p1, p2 Point, ratio float64) Point {
	return Point{
		X: int(float64(p1.X)*(1-ratio) + float64(p2.X)*ratio),
		Y: int(float64(p1.Y)*(1-ratio) + float64(p2.Y)*ratio),
	}
}

// Helper function to calculate the Koch vertex
func calculateKochVertex(p1, p2 Point) Point {
	return Point{
		X: p1.X + (p2.X-p1.X)/2,
		Y: p1.Y + int(math.Sqrt(3)*float64(p2.X-p1.X)/2),
	}
}

func (ppm *PPM) DrawSierpinskiTriangle(n int, start Point, width int, color Pixel) {
	// Define the three vertices of the initial equilateral triangle
	p1 := start
	p2 := Point{start.X + width, start.Y}
	p3 := Point{start.X + width / 2, start.Y - int(float64(width) * (3 ** 0.5) / 2)}

	// Call the recursive function to draw the Sierpinski triangle
	ppm.drawSierpinskiTriangleRecursive(n, p1, p2, p3, color)
}

// Helper function to draw a Sierpinski triangle recursively
func (ppm *PPM) drawSierpinskiTriangleRecursive(n int, p1, p2, p3 Point, color Pixel) {
	if n == 0 {
		// Draw the triangle
		ppm.drawLine(p1, p2, color)
		ppm.drawLine(p2, p3, color)
		ppm.drawLine(p3, p1, color)
		return
	}

	// Calculate midpoints of the edges
	mid1 := calculateMidpoint(p1, p2)
	mid2 := calculateMidpoint(p2, p3)
	mid3 := calculateMidpoint(p3, p1)
 
	// Recursively draw three smaller triangles
	ppm.drawSierpinskiTriangleRecursive(n-1, p1, mid1, mid3, color)
	ppm.drawSierpinskiTriangleRecursive(n-1, mid1, p2, mid2, color)
	ppm.drawSierpinskiTriangleRecursive(n-1, mid3, mid2, p3, color)
}

// Helper function to draw a line segment
func (ppm *PPM) drawLine(p1, p2 Point, color Pixel) {
	// You can implement the line drawing logic here
	// Update the image or perform other necessary actions
}

// Helper function to calculate the midpoint of two points
func calculateMidpoint(p1, p2 Point) Point {
	return Point{
		X: (p1.X + p2.X) / 2,
		Y: (p1.Y + p2.Y) / 2,
	}
} 

func (ppm *PPM) DrawPerlinNoise(color1 Pixel, color2 Pixel) {
	// Create a Perlin noise generator
	noiseGenerator := perlin.NewPerlin(1, 1, 2, rand.Int63())

	// Iterate over each pixel and set its color based on Perlin noise
	for y := 0; y < ppm.Height; y++ {
		for x := 0; x < ppm.Width; x++ {
			// Generate Perlin noise value
			noiseValue := noiseGenerator.Noise2D(float64(x)/50, float64(y)/50)

			// Interpolate between color1 and color2 based on the noise value
			interpolatedColor := interpolateColor(color1, color2, (noiseValue+1)/2)

			// Set the pixel color
			ppm.Pixels[y][x] = interpolatedColor
		}
	}
}

// Helper function to interpolate between two colors based on a factor
func interpolateColor(color1 Pixel, color2 Pixel, factor float64) Pixel {
	return Pixel{
		R: uint8(float64(color1.R)*(1-factor) + float64(color2.R)*factor),
		G: uint8(float64(color1.G)*(1-factor) + float64(color2.G)*factor),
		B: uint8(float64(color1.B)*(1-factor) + float64(color2.B)*factor),
	}
}

func (ppm *PPM) KNearestNeighbors(newWidth, newHeight int) {
	// Create a new image with the desired dimensions
	resizedImage := PPM{
		Width:  newWidth,
		Height: newHeight,
		Pixels: make([][]Pixel, newHeight),
	}

	for i := range resizedImage.Pixels {
		resizedImage.Pixels[i] = make([]Pixel, newWidth)
	}

	// Resize the image using k-nearest neighbors algorithm
	for y := 0; y < newHeight; y++ {
		for x := 0; x < newWidth; x++ {
			// Calculate the corresponding position in the original image
			originalX := int(float64(x) * float64(ppm.Width-1) / float64(newWidth-1))
			originalY := int(float64(y) * float64(ppm.Height-1) / float64(newHeight-1))

			// Perform k-nearest neighbors averaging
			kNeighbors := ppm.getKNearestNeighbors(originalX, originalY, 3)
			averageIntensity := ppm.calculateAverageIntensity(kNeighbors)

			// Set the pixel intensity in the resized image
			resizedImage.Pixels[y][x].Intensity = averageIntensity
		}
	}

	// Update the original image with the resized image
	ppm.Width = newWidth
	ppm.Height = newHeight
	ppm.Pixels = resizedImage.Pixels
}

// Helper function to get k-nearest neighbors for a given pixel position
func (ppm *PPM) getKNearestNeighbors(x, y, k int) []Pixel {
	neighbors := make([]Pixel, 0)

	for i := -k/2; i <= k/2; i++ {
		for j := -k/2; j <= k/2; j++ {
			// Avoid going out of bounds
			neighborX := x + i
			neighborY := y + j

			if neighborX >= 0 && neighborX < ppm.Width && neighborY >= 0 && neighborY < ppm.Height {
				neighbors = append(neighbors, ppm.Pixels[neighborY][neighborX])
			}
		}
	}

	return neighbors
}

// Helper function to calculate the average intensity of pixels
func (ppm *PPM) calculateAverageIntensity(neighbors []Pixel) uint8 {
	var sum uint16
	for _, neighbor := range neighbors {
		sum += uint16(neighbor.Intensity)
	}
	return uint8(sum / uint16(len(neighbors)))
} 
