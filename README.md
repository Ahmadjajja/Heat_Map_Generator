# Montana Heat Map Generator

A Python desktop application for generating hexagonal heat maps of Montana based on geographical data points.

## Features

- Load Excel files containing latitude and longitude data points
- Generate customizable hexagonal grid overlays
- Interactive color range configuration
- Real-time map visualization
- High-resolution TIFF export
- User-friendly GUI interface

## Requirements

- Python 3.8 or higher
- Required packages listed in `requirements.txt`

## Installation

1. Clone this repository
2. Create a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Run the application:
   ```bash
   python montana_heatmap.py
   ```

2. Load your Excel file:
   - Click "Browse" to select your Excel file
   - File must contain 'lat' and 'long' columns

3. Configure the heat map:
   - Set the number of hexagons to divide the map into
   - Adjust color ranges and corresponding colors
   - Colors can be specified by name or hex code

4. Generate and export:
   - Click "Generate Heat Map" to create the visualization
   - Click "Download Heat Map" to save as TIFF

## Input Data Format

Your Excel file should contain the following columns:
- `lat`: Latitude values
- `long`: Longitude values
- Optional: `county`, `year`, etc. for future filtering capabilities

## Directory Structure

```
montana_heatmap/
├── montana_heatmap.py    # Main application file
├── requirements.txt      # Python dependencies
├── README.md            # This file
└── shapefiles/          # Directory for Montana boundary files
    └── montana.shp      # Montana state boundary shapefile
```

## Notes

- The application requires Montana boundary shapefiles to be present in the `shapefiles` directory
- For best results, ensure your Excel data points fall within Montana's boundaries
- The application automatically adjusts to window resizing

## License

MIT License 