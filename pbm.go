package main

import (
	"bufio"
	"fmt"
	"os"
)

type PBM struct {
	data          [][]bool
	width, height int
	magicNumber   string
}

// ReadPBM reads a PBM image from a file and returns a struct that represents the image.
func ReadPBM(filename string) (*PBM, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read magic number (P1 or P4)
	scanner.Scan()
	magicNumber := scanner.Text()

	// Determine the format based on the magic number
	var isP1 bool
	if magicNumber == "P1" {
		isP1 = true
	} else if magicNumber == "P4" {
		isP1 = false
	} else {
		return nil, fmt.Errorf("unsupported PBM format: %s", magicNumber)
	}

	// Read width and height
	scanner.Scan()
	width, height := 0, 0
	fmt.Sscanf(scanner.Text(), "%d %d", &width, &height)

	// Read image data
	var data [][]bool
	for i := 0; i < height; i++ {
		scanner.Scan()
		line := scanner.Text()
		var row []bool

		if isP1 {
			// P1 format
			for _, char := range line {
				if char == '0' {
					row = append(row, false)
				} else if char == '1' {
					row = append(row, true)
				}
			}
		} else {
			// P4 format
			for _, char := range line {
				for j := 7; j >= 0; j-- {
					bit := (char >> uint(j)) & 1
					row = append(row, bit == 1)
				}
			}
		}

		data = append(data, row)
	}

	return &PBM{
		data:        data,
		width:       width,
		height:      height,
		magicNumber: magicNumber,
	}, nil
}

// Size returns the width and height of the image.
func (pbm *PBM) Size() (int, int) {
	return pbm.width, pbm.height
}

// At returns the value of the pixel at (x, y).
func (pbm *PBM) At(x, y int) bool {
	return pbm.data[y][x]
}

// Set sets the value of the pixel at (x, y).
func (pbm *PBM) Set(x, y int, value bool) {
	pbm.data[y][x] = value
}

// Save saves the PBM image to a file and returns an error if there was a problem.
func (pbm *PBM) Save(filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// Write magic number, width, and height
	fmt.Fprintf(writer, "%s\n%d %d\n", pbm.magicNumber, pbm.width, pbm.height)

	// Write image data
	for _, row := range pbm.data {
		if pbm.magicNumber == "P1" {
			// P1 format
			for _, pixel := range row {
				if pixel {
					fmt.Fprint(writer, "1 ")
				} else {
					fmt.Fprint(writer, "0 ")
				}
			}
		} else {
			// P4 format
			for i := 0; i < len(row); i += 8 {
				var byteValue byte
				for j := 0; j < 8; j++ {
					if i+j < len(row) && row[i+j] {
						byteValue |= 1 << uint(7-j)
					}
				}
				fmt.Fprintf(writer, "%c", byteValue)
			}
		}
		fmt.Fprintln(writer)
	}

	return writer.Flush()
}

// Invert inverts the colors of the PBM image.
func (pbm *PBM) Invert() {
	for y := 0; y < pbm.height; y++ {
		for x := 0; x < pbm.width; x++ {
			pbm.data[y][x] = !pbm.data[y][x]
		}
	}
}

// Flip flips the PBM image horizontally.
func (pbm *PBM) Flip() {
	for y := 0; y < pbm.height; y++ {
		for x := 0; x < pbm.width/2; x++ {
			pbm.data[y][x], pbm.data[y][pbm.width-x-1] = pbm.data[y][pbm.width-x-1], pbm.data[y][x]
		}
	}
}

// Flop flops the PBM image vertically.
func (pbm *PBM) Flop() {
	for y := 0; y < pbm.height/2; y++ {
		pbm.data[y], pbm.data[pbm.height-y-1] = pbm.data[pbm.height-y-1], pbm.data[y]
	}
}

// SetMagicNumber sets the magic number of the PBM image.
func (pbm *PBM) SetMagicNumber(magicNumber string) {
	pbm.magicNumber = magicNumber
}

// Exemple d'usage
func main() {
	filename := "example.pbm"
	pbm, err := ReadPBM(filename)
	if err != nil {
		fmt.Println("Error:", err)
		return
	}

	fmt.Println("PBM Image:")
	fmt.Println("Magic Number:", pbm.magicNumber)
	fmt.Println("Width:", pbm.width)
	fmt.Println("Height:", pbm.height)
	fmt.Println("Data:", pbm.data)

	// Example usage of other functions
	width, height := pbm.Size()
	fmt.Printf("Image Size: %d x %d\n", width, height)

	value := pbm.At(2, 3)
	fmt.Printf("Value at (2, 3): %t\n", value)

	pbm.Set(2, 3, true)
	fmt.Println("After setting value at (2, 3) to true:", pbm.data)

	err = pbm.Save("output.pbm")
	if err != nil {
		fmt.Println("Error saving PBM image:", err)
		return
	}

	fmt.Println("Image saved successfully.")
}
