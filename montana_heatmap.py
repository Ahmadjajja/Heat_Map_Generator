import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point, box
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import os
from typing import Dict, List, Tuple, Optional
import re

class SplashScreen:
    def __init__(self, parent):
        self.parent = parent
        self.splash = tk.Toplevel(parent)
        self.splash.title("Montana Heat Map Generator")
        
        # Get screen dimensions
        screen_width = self.splash.winfo_screenwidth()
        screen_height = self.splash.winfo_screenheight()
        
        # Calculate position
        width = 400
        height = 200
        x = (screen_width - width) // 2
        y = (screen_height - height) // 2
        
        self.splash.geometry(f"{width}x{height}+{x}+{y}")
        self.splash.overrideredirect(True)
        self.splash.configure(bg='white')
        
        # Add loading text
        self.status_label = tk.Label(
            self.splash,
            text="Initializing...",
            bg='white',
            font=('Arial', 12)
        )
        self.status_label.pack(pady=20)
        
        # Add progress bar
        self.progress = ttk.Progressbar(
            self.splash,
            length=300,
            mode='determinate'
        )
        self.progress.pack(pady=20)
        
        self.splash.update()

    def update_status(self, message: str, progress: int = None):
        self.status_label.config(text=message)
        if progress is not None:
            self.progress['value'] = progress
        self.splash.update()

    def destroy(self):
        self.splash.destroy()

class ToastNotification:
    def __init__(self, parent):
        self.parent = parent
        
    def show_toast(self, message: str, duration: int = 3000, error: bool = False):
        toast = tk.Toplevel(self.parent)
        toast.overrideredirect(True)
        
        # Position toast at bottom right
        toast.geometry(f"+{self.parent.winfo_screenwidth() - 310}+{self.parent.winfo_screenheight() - 100}")
        
        # Configure toast appearance
        bg_color = '#ff4444' if error else '#44aa44'
        frame = tk.Frame(toast, bg=bg_color, padx=10, pady=5)
        frame.pack(fill='both', expand=True)
        
        tk.Label(
            frame,
            text=message,
            bg=bg_color,
            fg='white',
            wraplength=250,
            font=('Arial', 10)
        ).pack()
        
        toast.after(duration, toast.destroy)

class MainApplication:
    def __init__(self):
        self.root = tk.Tk()
        self.root.withdraw()  # Hide main window initially
        
        # Show splash screen
        self.splash = SplashScreen(self.root)
        self.splash.update_status("Loading application...", 0)
        
        # Initialize variables
        self.excel_data = None
        self.montana_gdf = None
        self.hexagons = None
        self.current_map = None
        
        # Configure main window
        self.root.title("Montana Heat Map Generator")
        self.root.state('zoomed')  # Start maximized
        
        # Initialize notification system
        self.toast = ToastNotification(self.root)
        
        # Set up the GUI
        self.initialize_gui()
        
        # Destroy splash screen and show main window
        self.splash.destroy()
        self.root.deiconify()

    def initialize_gui(self):
        # Configure style
        style = ttk.Style()
        style.configure('TFrame', background='white')
        style.configure('TLabel', background='white')
        style.configure('TButton', padding=5)
        
        # Main container
        self.main_container = ttk.Frame(self.root)
        self.main_container.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Left panel (inputs)
        self.left_panel = ttk.Frame(self.main_container, style='TFrame')
        self.left_panel.pack(side='left', fill='y', padx=(0, 10))
        
        # Right panel (map display)
        self.right_panel = ttk.Frame(self.main_container, style='TFrame')
        self.right_panel.pack(side='right', fill='both', expand=True)
        
        self._setup_input_fields()
        self._setup_map_display()
        
        # Bind resize event
        self.root.bind('<Configure>', self.on_window_resize)

    def _setup_input_fields(self):
        # File selection
        ttk.Label(self.left_panel, text="Excel File:").pack(anchor='w', pady=(0, 5))
        self.file_frame = ttk.Frame(self.left_panel)
        self.file_frame.pack(fill='x', pady=(0, 20))
        
        self.file_path_var = tk.StringVar()
        ttk.Entry(self.file_frame, textvariable=self.file_path_var, state='readonly').pack(side='left', fill='x', expand=True)
        ttk.Button(self.file_frame, text="Browse", command=self.load_excel).pack(side='right', padx=(5, 0))
        
        # Hexagon count with preview button
        ttk.Label(self.left_panel, text="Number of Hexagons:").pack(anchor='w', pady=(0, 5))
        hex_frame = ttk.Frame(self.left_panel)
        hex_frame.pack(fill='x', pady=(0, 20))
        
        self.hex_count_var = tk.StringVar(value="100")
        ttk.Entry(hex_frame, textvariable=self.hex_count_var).pack(side='left', fill='x', expand=True)
        ttk.Button(hex_frame, text="Preview Grid", command=self.preview_grid).pack(side='right', padx=(5, 0))
        
        # Color ranges
        self.color_ranges = []
        default_ranges = [
            (0, 0, "white"),
            (1, 15, "yellow"),
            (16, 100, "orange"),
            (101, 105, "red"),
            (106, 120, "purple"),
            (121, float('inf'), "black")
        ]
        
        ttk.Label(self.left_panel, text="Color Ranges:").pack(anchor='w', pady=(0, 5))
        
        for i, (min_val, max_val, color) in enumerate(default_ranges):
            range_frame = ttk.Frame(self.left_panel)
            range_frame.pack(fill='x', pady=(0, 10))
            
            min_var = tk.StringVar(value=str(min_val))
            max_var = tk.StringVar(value="∞" if max_val == float('inf') else str(max_val))
            color_var = tk.StringVar(value=color)
            
            ttk.Entry(range_frame, textvariable=min_var, width=8).pack(side='left', padx=(0, 5))
            ttk.Label(range_frame, text="-").pack(side='left', padx=5)
            ttk.Entry(range_frame, textvariable=max_var, width=8).pack(side='left', padx=(5, 10))
            
            color_entry = ttk.Entry(range_frame, textvariable=color_var)
            color_entry.pack(side='left', fill='x', expand=True)
            
            self.color_ranges.append((min_var, max_var, color_var))
        
        # Action buttons
        ttk.Button(self.left_panel, text="Generate Heat Map", command=self.generate_map).pack(fill='x', pady=(20, 5))
        ttk.Button(self.left_panel, text="Download Heat Map", command=self.download_map).pack(fill='x', pady=(5, 0))

    def _setup_map_display(self):
        self.figure = Figure(figsize=(10, 8))
        self.ax = self.figure.add_subplot(111)
        # Remove the box from initial display
        self.ax.set_frame_on(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.right_panel)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)

    def load_excel(self):
        file_path = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
        )
        if not file_path:
            return
            
        try:
            self.excel_data = pd.read_excel(file_path)
            required_columns = ['lat', 'long']
            if not all(col in self.excel_data.columns for col in required_columns):
                raise ValueError("Excel file must contain 'lat' and 'long' columns")
                
            self.file_path_var.set(file_path)
            # print("Excel file loaded successfully")
            # print("self.excel_data: ", self.excel_data)
            print("self.excel_data.columns: ", self.excel_data.columns)
            self.toast.show_toast("Excel file loaded successfully")
        except Exception as e:
            self.toast.show_toast(f"Error loading file: {str(e)}", error=True)

    def generate_hexagonal_grid(self, bounds: Tuple[float, float, float, float], n_hexagons: int) -> gpd.GeoDataFrame:
        """Generate a hexagonal grid covering Montana."""
        xmin, ymin, xmax, ymax = bounds
        
        # Calculate the width and height of the area
        width = xmax - xmin
        height = ymax - ymin
        
        # Add small padding to ensure we catch edge hexagons
        padding_x = width * 0.05  # 5% padding
        padding_y = height * 0.05
        
        # Adjust bounds with padding
        xmin -= padding_x
        xmax += padding_x
        ymin -= padding_y
        ymax += padding_y
        
        # Recalculate width and height with padding
        width = xmax - xmin
        height = ymax - ymin
        
        # Calculate the target area for each hexagon
        total_area = width * height
        hex_target_area = total_area / n_hexagons
        
        # Calculate hexagon size (radius) based on target area
        # Area of a hexagon = 2 * sqrt(3) * r^2
        # where r is the radius (distance from center to vertex)
        hex_size = np.sqrt(hex_target_area / (2 * np.sqrt(3)))
        
        # Calculate horizontal and vertical spacing
        w = hex_size * 2  # Width of a hexagon (point to point)
        h = w * np.sqrt(3)/2  # Height of a hexagon (flat to flat)
        
        # Calculate number of hexagons needed in each direction to ensure coverage
        nx = int(np.ceil(width / (w * 0.866))) + 2  # Add buffer
        ny = int(np.ceil(height / (h * 0.866))) + 2  # Add buffer
        
        # Adjust the starting position to center the grid
        x_start = xmin - (nx * w * 0.866 - width) / 2
        y_start = ymin - (ny * h * 0.866 - height) / 2
        
        hexagons = []
        for row in range(ny):
            for col in range(nx):
                # Calculate center coordinates
                center_x = x_start + w * 0.866 * col  # cos(30°) = 0.866, for perfect hexagon spacing
                center_y = y_start + (h * 0.866) * row  # Using 0.866 for perfect vertical alignment
                
                # Offset even rows by half the width
                if row % 2 == 0:
                    center_x += w * 0.433  # Half of 0.866 for perfect offset
                
                # Generate vertices starting from rightmost point, going counterclockwise
                vertices = []
                for angle in range(0, 360, 60):
                    # Start from 30 degrees to point hexagon up
                    rad = np.radians(angle + 30)
                    vx = center_x + hex_size * np.cos(rad)
                    vy = center_y + hex_size * np.sin(rad)
                    vertices.append((vx, vy))
                vertices.append(vertices[0])  # Close the polygon
                
                hex_polygon = Polygon(vertices)
                
                # Only add hexagons that intersect with Montana's boundary
                if hex_polygon.intersects(self.montana_gdf.iloc[0].geometry):
                    hexagons.append(hex_polygon)
        
        hex_gdf = gpd.GeoDataFrame(geometry=hexagons, crs=self.montana_gdf.crs)
        return hex_gdf

    def preview_grid(self):
        try:
            # Convert hex count to int
            n_hexagons = int(self.hex_count_var.get())
            if n_hexagons <= 0:
                raise ValueError("Number of hexagons must be positive")
            
            # Load Montana boundary if not already loaded
            if self.montana_gdf is None:
                # Load US counties and filter for Montana
                all_counties = gpd.read_file("shapefiles/cb_2021_us_county_5m.shp")
                self.montana_gdf = all_counties[all_counties['STATEFP'] == '30']  # Montana's FIPS code is 30
                # Convert to a better projection for Montana (NAD 83 / Montana State Plane)
                self.montana_gdf = self.montana_gdf.to_crs("EPSG:32100")
                # Dissolve counties to get state boundary
                self.montana_gdf = self.montana_gdf.dissolve()
            
            # Generate hexagonal grid
            bounds = self.montana_gdf.total_bounds
            self.hexagons = self.generate_hexagonal_grid(bounds, n_hexagons)
            
            # Plot the preview
            self.ax.clear()
            # Remove the box
            self.ax.set_frame_on(False)
            
            # Plot state boundary with no fill and thin black line
            self.montana_gdf.boundary.plot(ax=self.ax, color='black', linewidth=0.5)
            # Plot hexagons with white fill and light gray edges
            self.hexagons.boundary.plot(ax=self.ax, color='gray', linewidth=0.5)
            
            # Set proper aspect ratio and limits
            self.ax.set_aspect('equal')
            
            # Add padding to the bounds
            bounds = self.montana_gdf.total_bounds
            padding = (bounds[2] - bounds[0]) * 0.15  # 15% padding
            self.ax.set_xlim([bounds[0] - padding, bounds[2] + padding])
            self.ax.set_ylim([bounds[1] - padding, bounds[3] + padding])
            
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.figure.tight_layout()
            self.canvas.draw()
            
            self.toast.show_toast(f"Preview grid with {n_hexagons} hexagons generated")
            
        except Exception as e:
            self.toast.show_toast(f"Error generating preview: {str(e)}", error=True)

    def dms_to_decimal(self, coord):
        """
        Convert a coordinate in DMS format (e.g., '44°41.576'') to decimal degrees.
        Handles both unicode and ascii degree/minute/second symbols.
        """
        if isinstance(coord, float) or isinstance(coord, int):
            return float(coord)
        if not isinstance(coord, str):
            return float('nan')
        # Remove unwanted characters and normalize
        coord = coord.replace("'", "'").replace("″", '"').replace("""", '"').replace(""", '"')
        dms_pattern = r"(\d+)[°\s]+(\d+(?:\.\d+)?)[\'′]?\s*(\d*(?:\.\d+)?)[\"″]?"
        match = re.match(dms_pattern, coord.strip())
        if match:
            deg = float(match.group(1))
            min_ = float(match.group(2))
            sec = float(match.group(3)) if match.group(3) else 0.0
            return deg + min_ / 60 + sec / 3600
        try:
            return float(coord)
        except Exception:
            return float('nan')

    def convert_coordinates(self, row):
        """Convert coordinates taking into account direction (N/S, E/W) and DMS/decimal formats"""
        try:
            lat = self.dms_to_decimal(row['lat'])
            long = self.dms_to_decimal(row['long'])
            
            # Convert direction values to string and handle potential NaN/float values
            lat_dir = str(row['lat_dir']).strip().upper() if pd.notna(row['lat_dir']) else 'N'
            long_dir = str(row['long_dir']).strip().upper() if pd.notna(row['long_dir']) else 'W'
            
            # Validate direction values
            if lat_dir not in ['N', 'S']:
                print(f"Invalid latitude direction: {lat_dir}, defaulting to 'N'")
                lat_dir = 'N'
            if long_dir not in ['E', 'W']:
                print(f"Invalid longitude direction: {long_dir}, defaulting to 'W'")
                long_dir = 'W'
            
            # Adjust for direction
            if lat_dir == 'S':  # If Southern hemisphere
                lat = -lat
            if long_dir == 'W':  # If Western hemisphere
                long = -long
            
            # Montana is roughly between 44°N to 49°N and 104°W to 116°W
            # Validate the coordinates are somewhat reasonable
            if not (44 <= abs(lat) <= 49 and 104 <= abs(long) <= 116):
                print(f"Warning: Coordinates ({lat}, {long}) might be outside Montana's bounds")
            
            return Point(long, lat)
        except Exception as e:
            print(f"Error converting coordinates: {str(e)}")
            # Return a point outside Montana's bounds which will be filtered out
            return Point(0, 0)

    def generate_map(self):
        if self.excel_data is None:
            self.toast.show_toast("Please load an Excel file first", error=True)
            return
            
        try:
            # Convert hex count to int
            n_hexagons = int(self.hex_count_var.get())
            if n_hexagons <= 0:
                raise ValueError("Number of hexagons must be positive")
                
            # Validate required columns
            required_columns = ['lat', 'lat_dir', 'long', 'long_dir']
            if not all(col in self.excel_data.columns for col in required_columns):
                raise ValueError("Excel file must contain 'lat', 'lat_dir', 'long', and 'long_dir' columns")
            
            # Create points from Excel data with direction consideration
            geometries = self.excel_data.apply(self.convert_coordinates, axis=1)
            points = gpd.GeoDataFrame(
                self.excel_data,
                geometry=geometries,
                crs="EPSG:4326"
            )
            
            # Filter points to only those inside the actual Montana polygon
            montana_poly = self.montana_gdf.to_crs("EPSG:4326").geometry.iloc[0]
            points = points[points.geometry.within(montana_poly)]
            
            if len(points) == 0:
                self.toast.show_toast("No points found within Montana's boundaries", error=True)
                return
            
            # Convert points to the same CRS as Montana
            points = points.to_crs(self.montana_gdf.crs)
            
            # Generate hexagonal grid if not already generated
            if self.hexagons is None:
                # Load Montana boundary if not already loaded
                if self.montana_gdf is None:
                    # Load US counties and filter for Montana
                    all_counties = gpd.read_file("shapefiles/cb_2021_us_county_5m.shp")
                    self.montana_gdf = all_counties[all_counties['STATEFP'] == '30']  # Montana's FIPS code is 30
                    # Convert to a better projection for Montana (NAD 83 / Montana State Plane)
                    self.montana_gdf = self.montana_gdf.to_crs("EPSG:32100")
                    # Dissolve counties to get state boundary
                    self.montana_gdf = self.montana_gdf.dissolve()
                
                # Generate hexagonal grid
                bounds = self.montana_gdf.total_bounds
                self.hexagons = self.generate_hexagonal_grid(bounds, n_hexagons)
            
            # Initialize point_count column with zeros
            self.hexagons['point_count'] = 0
            
            # Count points in each hexagon
            for idx, hexagon in self.hexagons.iterrows():
                # Count points that fall within this hexagon
                # print("idx: ", idx)
                # print("hexagon: ", hexagon)
                # print("points: ", points)
                # print("hexagon.geometry: ", hexagon.geometry)
                # print("points.within(hexagon.geometry): ", points.within(hexagon.geometry))
                # print("points.within(hexagon.geometry).count(): ", points.within(hexagon.geometry).count())
                points_in_hex = points[points.within(hexagon.geometry)]
                self.hexagons.at[idx, 'point_count'] = len(points_in_hex)
            
            # Assign colors based on ranges
            self.hexagons['color'] = None  # Initialize with None
            
            # Sort ranges by min value to ensure proper order of application
            ranges = []
            for min_var, max_var, color_var in self.color_ranges:
                min_val = float(min_var.get())
                max_val = float('inf') if max_var.get() == "∞" else float(max_var.get())
                ranges.append((min_val, max_val, color_var.get()))
            
            # Sort ranges by minimum value
            ranges.sort(key=lambda x: x[0])
            
            # Apply colors based on point counts
            for min_val, max_val, color in ranges:
                # Create mask for hexagons within this range
                mask = (self.hexagons['point_count'] >= min_val) & (self.hexagons['point_count'] <= max_val)
                self.hexagons.loc[mask, 'color'] = color
            
            # Plot the map
            self.ax.clear()
            # Remove the box
            self.ax.set_frame_on(False)
            
            # Plot state boundary with no fill and thin black line
            self.montana_gdf.boundary.plot(ax=self.ax, color='black', linewidth=0.5)
            # Plot hexagons with colors and very thin edges
            for idx, hexagon in self.hexagons.iterrows():
                color = hexagon['color']
                if color:
                    self.ax.fill(hexagon.geometry.exterior.xy[0], 
                               hexagon.geometry.exterior.xy[1],
                               facecolor=color,
                               edgecolor='gray',
                               linewidth=0.1,
                               alpha=0.7)
            
            # Add legend
            legend_elements = []
            for min_var, max_var, color_var in self.color_ranges:
                min_val = min_var.get()
                max_val = "∞" if max_var.get() == "∞" else max_var.get()
                label = f"{min_val}-{max_val}"
                legend_elements.append(plt.Rectangle((0, 0), 1, 1, fc=color_var.get()))
            
            self.ax.legend(legend_elements,
                         [f"{min_var.get()}-{max_var.get()}" for min_var, max_var, _ in self.color_ranges],
                         title="Point Count Ranges",
                         loc='center left',
                         bbox_to_anchor=(1, 0.5))
            
            # Set proper aspect ratio and limits
            self.ax.set_aspect('equal')
            
            # Add padding to the bounds
            bounds = self.montana_gdf.total_bounds
            padding = (bounds[2] - bounds[0]) * 0.15  # 15% padding
            self.ax.set_xlim([bounds[0] - padding, bounds[2] + padding])
            self.ax.set_ylim([bounds[1] - padding, bounds[3] + padding])
            
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.figure.tight_layout()
            self.canvas.draw()
            
            self.toast.show_toast("Map generated successfully")
            
        except Exception as e:
            self.toast.show_toast(f"Error generating map: {str(e)}", error=True)

    def download_map(self):
        if self.hexagons is None:
            self.toast.show_toast("Please generate a map first", error=True)
            return
            
        try:
            import datetime
            from pathlib import Path
            import os
            
            # Get Downloads folder path
            downloads_path = str(Path.home() / "Downloads")
            
            # Get current date and time in the desired format
            now = datetime.datetime.now()
            timestamp = now.strftime("%I_%M_%p_%m_%d_%Y")  # e.g., 12_49_PM_6_12_2025
            
            # Create a meaningful filename
            filename = f"MontanaHeatMap_{timestamp}.tiff"
            file_path = os.path.join(downloads_path, filename)
            
            # Save the figure
            self.figure.savefig(file_path, format="tiff", dpi=300, bbox_inches='tight')
            
            # Show toast notification
            self.toast.show_toast(f"Map saved as {filename}")
            
            print(f"✅ TIFF map saved as '{file_path}'")
            
        except Exception as e:
            messagebox.showerror("Error", 
                f"Error saving file:\n{str(e)}\n\n"
                "Please try again."
            )

    def on_window_resize(self, event=None):
        # Update the figure size to match the panel size
        w = self.right_panel.winfo_width() / 100
        h = self.right_panel.winfo_height() / 100
        self.figure.set_size_inches(w, h)
        self.canvas.draw()

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = MainApplication()
    app.run() 