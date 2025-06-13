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

# [Previous classes remain the same until MainApplication]

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
        
        # Add variables for species selection
        self.selected_family = tk.StringVar()
        self.selected_genus = tk.StringVar()
        self.selected_species = tk.StringVar()
        
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

    def calculate_hexagon_dimensions(self, hex_gdf):
        """Calculate the dimensions of a hexagon in miles and kilometers."""
        # Get the first hexagon's geometry
        hex_geom = hex_gdf.iloc[0].geometry
        
        # Calculate width (point to point)
        width = hex_geom.length / 6  # Since it's a regular hexagon, each side is equal
        
        # Calculate height (flat to flat)
        height = width * np.sqrt(3)
        
        # Convert to miles and kilometers
        width_miles = width * 0.000621371  # Convert meters to miles
        height_miles = height * 0.000621371
        width_km = width / 1000  # Convert meters to kilometers
        height_km = height / 1000
        
        return {
            'miles': {'width': width_miles, 'height': height_miles},
            'kilometers': {'width': width_km, 'height': height_km}
        }

    def preview_grid(self):
        try:
            loading = LoadingIndicator(self.root, "Generating preview grid...")
            n_hexagons = int(self.hex_count_var.get())
            if n_hexagons <= 0:
                loading.destroy()
                raise ValueError("Number of hexagons must be positive")
            if self.montana_gdf is None:
                loading.update_message("Loading Montana boundary...")
                all_counties = gpd.read_file("shapefiles/cb_2021_us_county_5m.shp")
                self.montana_gdf = all_counties[all_counties['STATEFP'] == '30']
                self.montana_gdf = self.montana_gdf.to_crs("EPSG:32100")
                self.montana_gdf = self.montana_gdf.dissolve()
            loading.update_message("Generating hexagonal grid...")
            bounds = self.montana_gdf.total_bounds
            self.hexagons = self.generate_hexagonal_grid(bounds, n_hexagons)
            
            # Calculate hexagon dimensions
            dimensions = self.calculate_hexagon_dimensions(self.hexagons)
            
            loading.update_message("Rendering preview...")
            self.figure.clf()
            self.ax = self.figure.add_subplot(111)
            self.ax.set_frame_on(False)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.ax.set_aspect('equal')
            self.montana_gdf.boundary.plot(ax=self.ax, color='black', linewidth=0.5)
            self.hexagons.boundary.plot(ax=self.ax, color='gray', linewidth=0.5)
            bounds = self.montana_gdf.total_bounds
            padding = (bounds[2] - bounds[0]) * 0.15
            self.ax.set_xlim([bounds[0] - padding, bounds[2] + padding])
            self.ax.set_ylim([bounds[1] - padding, bounds[3] + padding])
            self.figure.subplots_adjust(bottom=0.18, top=0.98, left=0.01, right=0.99)
            
            # Add title with hexagon dimensions
            title_text = (
                f"Preview Grid ({n_hexagons} hexagons)\n"
                f"Hexagon Dimensions:\n"
                f"Width: {dimensions['miles']['width']:.2f} mi / {dimensions['kilometers']['width']:.2f} km\n"
                f"Height: {dimensions['miles']['height']:.2f} mi / {dimensions['kilometers']['height']:.2f} km"
            )
            self.ax.text(0.5, 1, title_text, transform=self.ax.transAxes, fontsize=12, 
                        fontweight='normal', color='#2c3e50', ha='center', va='bottom',
                        bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.3', alpha=0.85))
            
            self.canvas.draw()
            loading.destroy()
            self.toast.show_toast(f"Preview grid with {n_hexagons} hexagons generated")
        except Exception as e:
            if 'loading' in locals():
                loading.destroy()
            self.toast.show_toast(f"Error generating preview: {str(e)}", error=True)

    def generate_map(self):
        if self.excel_data is None:
            self.toast.show_toast("Please load an Excel file first", error=True)
            return
        try:
            loading = LoadingIndicator(self.root, "Generating heat map...")
            n_hexagons = int(self.hex_count_var.get())
            if n_hexagons <= 0:
                loading.destroy()
                raise ValueError("Number of hexagons must be positive")
            required_columns = ['lat', 'lat_dir', 'long', 'long_dir', 'family', 'genus', 'species']
            if not all(col in self.excel_data.columns for col in required_columns):
                loading.destroy()
                raise ValueError("Excel file must contain 'lat', 'lat_dir', 'long', 'long_dir', 'family', 'genus', and 'species' columns")
            fam = self.selected_family.get().strip()
            gen = self.selected_genus.get().strip()
            spec = self.selected_species.get().strip()
            if not fam or fam == "Select Family" or not gen or gen == "Select Genus" or not spec or spec == "Select Species":
                loading.destroy()
                messagebox.showerror("Missing Input", "Please select Family, Genus, and Species.")
                return
            loading.update_message("Filtering data...")
            filtered = self.excel_data
            if fam == "All":
                filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
            else:
                filtered = filtered[filtered["family"].str.lower() == fam.lower()]
            if gen == "All":
                filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
            else:
                filtered = filtered[filtered["genus"].str.lower() == gen.lower()]
            if spec == "all":
                filtered = filtered[filtered["species"].notna() & (filtered["species"].str.strip() != "")]
            else:
                filtered = filtered[filtered["species"].str.lower() == spec.lower()]
            loading.update_message("Creating points...")
            geometries = filtered.apply(self.convert_coordinates, axis=1)
            points = gpd.GeoDataFrame(
                filtered,
                geometry=geometries,
                crs="EPSG:4326"
            )
            montana_poly = self.montana_gdf.to_crs("EPSG:4326").geometry.iloc[0]
            points = points[points.geometry.within(montana_poly)]
            if len(points) == 0:
                loading.destroy()
                self.toast.show_toast("No points found within Montana's boundaries", error=True)
                return
            points = points.to_crs(self.montana_gdf.crs)
            if self.hexagons is None:
                loading.update_message("Generating hexagonal grid...")
                if self.montana_gdf is None:
                    all_counties = gpd.read_file("shapefiles/cb_2021_us_county_5m.shp")
                    self.montana_gdf = all_counties[all_counties['STATEFP'] == '30']
                    self.montana_gdf = self.montana_gdf.to_crs("EPSG:32100")
                    self.montana_gdf = self.montana_gdf.dissolve()
                bounds = self.montana_gdf.total_bounds
                self.hexagons = self.generate_hexagonal_grid(bounds, n_hexagons)
            
            # Calculate hexagon dimensions
            dimensions = self.calculate_hexagon_dimensions(self.hexagons)
            
            loading.update_message("Counting points in hexagons...")
            self.hexagons['point_count'] = 0
            for idx, hexagon in self.hexagons.iterrows():
                points_in_hex = points[points.within(hexagon.geometry)]
                self.hexagons.at[idx, 'point_count'] = len(points_in_hex)
            loading.update_message("Assigning colors...")
            self.hexagons['color'] = None
            ranges = []
            for min_var, max_var, color_var in self.color_ranges:
                min_val = float(min_var.get())
                max_val = float('inf') if max_var.get() == "∞" else float(max_var.get())
                ranges.append((min_val, max_val, color_var.get()))
            ranges.sort(key=lambda x: x[0])
            for min_val, max_val, color in ranges:
                mask = (self.hexagons['point_count'] >= min_val) & (self.hexagons['point_count'] <= max_val)
                self.hexagons.loc[mask, 'color'] = color
            loading.update_message("Rendering map...")
            self.figure.clf()
            self.ax = self.figure.add_subplot(111)
            self.ax.set_frame_on(False)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.ax.set_aspect('equal')
            self.montana_gdf.boundary.plot(ax=self.ax, color='black', linewidth=0.5)
            for idx, hexagon in self.hexagons.iterrows():
                color = hexagon['color']
                if color:
                    self.ax.fill(hexagon.geometry.exterior.xy[0], 
                               hexagon.geometry.exterior.xy[1],
                               facecolor=color,
                               edgecolor='gray',
                               linewidth=0.1,
                               alpha=0.7)
            
            bounds = self.montana_gdf.total_bounds
            padding = (bounds[2] - bounds[0]) * 0.15
            self.ax.set_xlim([bounds[0] - padding, bounds[2] + padding])
            self.ax.set_ylim([bounds[1] - padding, bounds[3] + padding])
            self.figure.subplots_adjust(bottom=0.18, top=0.93, left=0.01, right=0.99)
            
            # Add label at the top of the map with species info and hexagon dimensions
            label_text = (
                f"{fam} > {gen} > {spec}\n"
                f"Hexagon Dimensions:\n"
                f"Width: {dimensions['miles']['width']:.2f} mi / {dimensions['kilometers']['width']:.2f} km\n"
                f"Height: {dimensions['miles']['height']:.2f} mi / {dimensions['kilometers']['height']:.2f} km"
            )
            self.ax.text(0.5, 1, label_text, transform=self.ax.transAxes, fontsize=14, 
                        fontweight='normal', color='#2c3e50', ha='center', va='bottom',
                        bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.3', alpha=0.85))
            
            import matplotlib.patches as mpatches
            legend_elements = []
            legend_labels = []
            for min_var, max_var, color_var in self.color_ranges:
                min_val = min_var.get()
                max_val = "∞" if max_var.get() == "∞" else max_var.get()
                label = f"{min_val}-{max_val}"
                legend_elements.append(mpatches.Patch(facecolor=color_var.get(), edgecolor='black'))
                legend_labels.append(label)
            self.ax.legend(
                legend_elements,
                legend_labels,
                title="Point Count Ranges",
                loc='lower center',
                bbox_to_anchor=(0.5, -0.18),
                ncol=len(legend_elements),
                frameon=True,
                fancybox=True,
                shadow=False,
                borderpad=1.2
            )
            self.canvas.draw()
            loading.destroy()
            self.toast.show_toast("Map generated successfully")
        except Exception as e:
            if 'loading' in locals():
                loading.destroy()
            self.toast.show_toast(f"Error generating map: {str(e)}", error=True)

    # [Rest of the class methods remain the same]

if __name__ == "__main__":
    app = MainApplication()
    app.run() 