# Montana Heat Map Generator

A powerful desktop application for visualizing and analyzing spatial data in Montana using a hexagon-based heat map. Built with Python, Tkinter, GeoPandas, Matplotlib, and Pandas, this tool is designed for researchers, ecologists, and anyone working with geospatial point data (e.g., bee sampling, biodiversity, field surveys).

---

## Features

- **Load Excel Data:** Supports `.xlsx` files with latitude/longitude in either decimal or DMS (degrees-minutes-seconds) format, and direction fields (`lat_dir`, `long_dir`).
- **Automatic Coordinate Parsing:** Handles both decimal and DMS formats, and applies N/S/E/W direction logic.
- **Montana-Only Filtering:** Only points that fall within the actual Montana state polygon are counted and visualized.
- **Hexagonal Grid Generation:**
  - User specifies the number of hexagons to tile Montana.
  - Grid is perfectly tessellated and clipped to the Montana boundary.
  - Border hexagons are included if they touch the state.
- **Point Counting and Coloring:**
  - Each hexagon is colored based on the number of points it contains, using user-defined color ranges (e.g., white for 0, yellow for 1–15, etc.).
  - Only points inside Montana are counted, even for border hexagons.
- **Interactive GUI:**
  - Modern, resizable interface with left-side controls and right-side live map preview.
  - Color ranges and hexagon count are fully customizable.
  - Toast notifications for user feedback.
- **Export:**
  - Download the generated map as a high-resolution TIFF file.
  - File is automatically saved to your Downloads folder with a timestamped, meaningful filename.
- **Robust Data Handling:**
  - Skips and warns about invalid or out-of-bounds coordinates.
  - Handles large datasets efficiently.

---

## Requirements

- Python 3.8+
- See `requirements.txt` for all dependencies:
  - pandas
  - geopandas
  - matplotlib
  - shapely
  - numpy
  - openpyxl
  - pillow

---

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/Ahmadjajja/Heat_Map_Generator.git
   cd Heat_Map_Generator
   ```
2. (Recommended) Create and activate a virtual environment:
   ```bash
   python -m venv venv
   venv\Scripts\activate  # On Windows
   # or
   source venv/bin/activate  # On Mac/Linux
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. Ensure you have the Montana county shapefile in the `shapefiles/` directory (see below).

---

## Usage

1. **Run the application:**
   ```bash
   python montana_heatmap.py
   ```
2. **Load your Excel file:**
   - Click "Browse" and select your `.xlsx` file.
   - Required columns: `lat`, `lat_dir`, `long`, `long_dir` (decimal or DMS format supported).
3. **Configure the map:**
   - Enter the number of hexagons for the grid.
   - Adjust color ranges and colors as desired.
   - Click "Preview Grid" to see the hex grid overlay.
4. **Generate the heat map:**
   - Click "Generate Heat Map" to color hexagons based on point density.
5. **Export:**
   - Click "Download Heat Map" to save the map as a TIFF in your Downloads folder. The filename will include the current date and time.

---

## Input Data Format

Your Excel file must include:
- `lat`: Latitude (decimal or DMS, e.g., `44.695` or `44°41.576'`)
- `lat_dir`: 'N' or 'S'
- `long`: Longitude (decimal or DMS, e.g., `-110.456` or `110°27.360'`)
- `long_dir`: 'E' or 'W'

Other columns (e.g., `county`, `year`, `species`) are ignored for mapping but can be present.

---

## How It Works

1. **Coordinate Parsing:**
   - All coordinates are parsed and converted to decimal degrees, with direction applied.
2. **Montana Filtering:**
   - Only points inside the actual Montana polygon are kept.
3. **Hex Grid Generation:**
   - The state is tiled with the specified number of hexagons, perfectly tessellated and clipped to the border.
4. **Point Counting:**
   - For each hexagon, the number of Montana points inside is counted.
5. **Color Assignment:**
   - Each hexagon is colored according to the user's color range settings.
6. **Export:**
   - The map can be saved as a TIFF with a timestamped filename.

---

## Example Output

- A clean, publication-ready map of Montana with a honeycomb grid, colored by point density.
- Only Montana data is visualized; out-of-state points are ignored.
- Border hexagons are included if they touch Montana.

---

## Shapefile Requirement

- Place the Montana county shapefile (e.g., `cb_2021_us_county_5m.shp` and related files) in a `shapefiles/` directory in your project root.
- The software will automatically extract the Montana boundary from this file.

---

## Troubleshooting

- **Invalid coordinates:** The app will warn and skip rows with unparseable or missing coordinates.
- **Points outside Montana:** The app will warn and skip points outside the state polygon.
- **No points in Montana:** If your data has no valid Montana points, the map will not be generated.
- **Git integration:** See the repo for version control and collaboration.

---

## License

MIT License

---

## Author

[Ahmadjajja/Heat_Map_Generator](https://github.com/Ahmadjajja/Heat_Map_Generator) 