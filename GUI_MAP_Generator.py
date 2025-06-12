import tkinter as tk
from tkinter import ttk, StringVar, filedialog, messagebox, Canvas
import os
from pathlib import Path
import datetime
import sys
import pandas as pd

def get_screen_geometry():
    """Get the geometry of all available screens"""
    root = tk.Tk()
    root.withdraw()  # Hide the temporary window
    
    # Get primary screen dimensions
    primary_width = root.winfo_screenwidth()
    primary_height = root.winfo_screenheight()
    
    # Get all screen dimensions using _tkinter.tcl_call
    try:
        # This gets info about all screens
        screens_info = root.tk.call('wm', 'maxsize', '.')
        all_screens_width = screens_info[0]
        all_screens_height = screens_info[1]
    except:
        # Fallback to primary screen if can't get all screens
        all_screens_width = primary_width
        all_screens_height = primary_height
    
    root.destroy()
    return primary_width, primary_height, all_screens_width, all_screens_height

class SplashScreen:
    def __init__(self, parent):
        self.root = tk.Toplevel(parent)
        self.root.overrideredirect(True)
        self.root.attributes('-alpha', 0.9)
        self.root.attributes('-topmost', True)
        
        # Get screen dimensions
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        
        # Window dimensions
        width = 400
        height = 200
        x = (screen_width - width) // 2
        y = (screen_height - height) // 2
        
        self.root.geometry(f'{width}x{height}+{x}+{y}')
        
        # Create frame
        frame = ttk.Frame(self.root, relief='raised', borderwidth=2)
        frame.pack(fill='both', expand=True)
        
        # Add title
        title = ttk.Label(
            frame, 
            text="Montana Bee County Map Generator",
            font=('Helvetica', 16, 'bold'),
            foreground='dark green'
        )
        title.pack(pady=20)
        
        # Add loading text
        self.status = ttk.Label(
            frame,
            text="Initializing...",
            font=('Helvetica', 10)
        )
        self.status.pack(pady=10)
        
        # Add progress bar
        self.progress = ttk.Progressbar(
            frame,
            mode='indeterminate',
            length=300
        )
        self.progress.pack(pady=20)
        self.progress.start(10)
        
        self.root.update()
    
    def update_status(self, message):
        self.status.config(text=message)
        self.root.update()
    
    def destroy(self):
        self.root.destroy()

class ToastNotification:
    def __init__(self, parent):
        self.parent = parent
        
    def show_toast(self, message, duration=3000, error=False):
        # Create toast window
        toast = tk.Toplevel(self.parent)
        toast.withdraw()  # Hide initially for position calculation
        
        # Remove window decorations
        toast.overrideredirect(True)
        
        # Create frame with border
        frame = ttk.Frame(toast, style='Toast.TFrame')
        frame.pack(fill='both', expand=True)
        
        # Add message with icon
        msg_frame = ttk.Frame(frame)
        msg_frame.pack(pady=10, padx=20)
        
        # Icon and color based on type
        if error:
            icon_text = "✕"
            icon_color = 'red'
            msg_color = 'red'
        else:
            icon_text = "✓"
            icon_color = 'green'
            msg_color = 'black'
        
        # Icon label
        icon_label = ttk.Label(msg_frame, text=icon_text, font=('Helvetica', 14, 'bold'), foreground=icon_color)
        icon_label.pack(side='left', padx=(0, 10))
        
        # Message with color
        msg_label = ttk.Label(msg_frame, text=message, font=('Helvetica', 10), foreground=msg_color)
        msg_label.pack(side='left')
        
        # Position toast at bottom right of main window
        toast.update_idletasks()  # Update to get actual size
        main_x = self.parent.winfo_x()
        main_y = self.parent.winfo_y()
        main_width = self.parent.winfo_width()
        main_height = self.parent.winfo_height()
        
        toast_width = toast.winfo_width()
        toast_height = toast.winfo_height()
        
        x = main_x + main_width - toast_width - 20
        y = main_y + main_height - toast_height - 20
        
        toast.geometry(f"+{x}+{y}")
        toast.deiconify()  # Show the toast
        
        # Destroy after duration
        toast.after(duration, toast.destroy)

class SelectionScreen:
    def __init__(self, parent, main_app, from_analysis=False):
        self.root = tk.Toplevel(parent)
        self.root.title("Analysis Selection")
        self.main_app = main_app
        
        # Set window icon
        if getattr(sys, 'frozen', False):
            base_dir = sys._MEIPASS
        else:
            base_dir = os.path.dirname(os.path.abspath(__file__))
        icon_path = os.path.join(base_dir, "app_icon.ico")
        if os.path.exists(icon_path):
            self.root.iconbitmap(icon_path)
        
        # Allow the window to be moved to any screen
        self.root.attributes('-alpha', 1.0)
        self.root.attributes('-topmost', False)
        
        # Get screen dimensions
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        
        # Set minimum size and make window resizable
        self.root.minsize(600, 400)
        self.root.resizable(True, True)
        
        # Set initial size to 80% of screen
        initial_width = int(screen_width * 0.8)
        initial_height = int(screen_height * 0.8)
        x_position = (screen_width - initial_width) // 2
        y_position = (screen_height - initial_height) // 2
        
        # Set initial geometry but don't enforce it
        self.root.geometry(f"{initial_width}x{initial_height}+{x_position}+{y_position}")
        
        # Update window to ensure proper sizing
        self.root.update_idletasks()
        
        # Bind window state change
        self.root.bind("<Configure>", self.on_window_resize)
        
        # Configure styles
        style = ttk.Style()
        style.configure('Title.TLabel',
                       font=('Helvetica', 32, 'bold'),
                       foreground='dark green')
        
        style.configure('Subtitle.TLabel',
                       font=('Helvetica', 24),
                       foreground='dark green')
        
        style.configure('Selection.TButton',
                       font=('Helvetica', 14, 'bold'),
                       padding=20)
        
        style.configure('Description.TLabel',
                       font=('Helvetica', 12),
                       foreground='gray')
        
        # Create main container with padding
        main_container = ttk.Frame(self.root, padding="40")
        main_container.pack(fill='both', expand=True)
        
        # Create a content frame to center everything vertically
        content_frame = ttk.Frame(main_container)
        content_frame.pack(expand=True)
        
        # Add logo or icon (if available)
        try:
            logo_path = "bee_logo.png"  # You can add a logo file later
            if os.path.exists(logo_path):
                logo = tk.PhotoImage(file=logo_path)
                logo_label = ttk.Label(content_frame, image=logo)
                logo_label.image = logo  # Keep a reference
                logo_label.pack(pady=(0, 30))
        except:
            pass
        
        # Title with shadow effect
        title_frame = ttk.Frame(content_frame)
        title_frame.pack(pady=(0, 40))
        
        title_label = ttk.Label(
            title_frame,
            text="Montana Bee County Map Generator",
            style='Title.TLabel'
        )
        title_label.pack()
        
        # Subtitle
        subtitle_label = ttk.Label(
            content_frame,
            text="Select Analysis Type",
            style='Subtitle.TLabel'
        )
        subtitle_label.pack(pady=(0, 60))
        
        # Button container for consistent width
        button_container = ttk.Frame(content_frame)
        button_container.pack(expand=True)
        
        # Single Year Analysis Section
        single_year_frame = ttk.Frame(button_container)
        single_year_frame.pack(pady=(0, 30))
        
        single_year_btn = ttk.Button(
            single_year_frame,
            text="Single Year Analysis",
            style='Selection.TButton',
            command=self.start_single_year_analysis,
            width=30
        )
        single_year_btn.pack()
        
        single_year_desc = ttk.Label(
            single_year_frame,
            text="Analyze bee distribution before and after a specific year",
            style='Description.TLabel'
        )
        single_year_desc.pack(pady=(10, 0))
        
        # Dual Year Analysis Section
        dual_year_frame = ttk.Frame(button_container)
        dual_year_frame.pack()
        
        dual_year_btn = ttk.Button(
            dual_year_frame,
            text="Dual Year Analysis",
            style='Selection.TButton',
            command=self.start_dual_year_analysis,
            width=30
        )
        dual_year_btn.pack()
        
        dual_year_desc = ttk.Label(
            dual_year_frame,
            text="Analyze bee distribution across three time periods using two year boundaries",
            style='Description.TLabel'
        )
        dual_year_desc.pack(pady=(10, 0))
        
        # Configure hover effects
        style.map('Selection.TButton',
                 background=[('active', '#e1e1e1')],
                 foreground=[('active', 'dark green')])
        
        # Add version or copyright info at bottom
        footer_frame = ttk.Frame(main_container)
        footer_frame.pack(side='bottom', fill='x')
        
        current_year = datetime.datetime.now().year
        copyright_label = ttk.Label(
            footer_frame,
            text=f"© {current_year} Montana Bee Project",
            font=('Helvetica', 8),
            foreground='gray'
        )
        copyright_label.pack(side='right')
    
    def start_single_year_analysis(self):
        # Store current position and state
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        current_state = self.root.state()
        self.root.destroy()
        analysis = SingleYearAnalysis(self.root.master, self.main_app)
        # Set position first
        analysis.root.geometry(f"{analysis.root.winfo_screenwidth()}x{analysis.root.winfo_screenheight()}+{x}+{y}")
        analysis.root.update()  # Force geometry update
        # Restore the previous window state
        if current_state == 'zoomed':
            analysis.root.state('zoomed')
    
    def start_dual_year_analysis(self):
        # Store current position and state
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        current_state = self.root.state()
        self.root.destroy()
        analysis = DualYearAnalysis(self.root.master, self.main_app)
        # Set position first
        analysis.root.geometry(f"{analysis.root.winfo_screenwidth()}x{analysis.root.winfo_screenheight()}+{x}+{y}")
        analysis.root.update()  # Force geometry update
        # Restore the previous window state
        if current_state == 'zoomed':
            analysis.root.state('zoomed')
    
    def on_window_resize(self, event=None):
        # Allow free movement and resizing
        pass

class MainApplication:
    def __init__(self):
        # Create and hide the root window
        self.root = tk.Tk()
        
        # Set the application icon
        if getattr(sys, 'frozen', False):
            # If running as exe
            base_dir = sys._MEIPASS
        else:
            # If running as script
            base_dir = os.path.dirname(os.path.abspath(__file__))
        
        icon_path = os.path.join(base_dir, "app_icon.ico")
        if os.path.exists(icon_path):
            self.root.iconbitmap(icon_path)
        
        # Get screen information
        self.primary_width, self.primary_height, self.all_screens_width, self.all_screens_height = get_screen_geometry()
        
        # Allow the window to be moved to any screen
        self.root.attributes('-alpha', 1.0)  # Ensure window is visible
        self.root.attributes('-topmost', False)  # Don't force window to stay on top
        
        # Hide initially for setup
        self.root.withdraw()
        
        # Initialize variables
        self.gdf = None
        self.pd = None
        self.gpd = None
        self.plt = None
        self.FigureCanvasTkAgg = None
        
        # Show splash screen
        self.splash = SplashScreen(self.root)
        
        # Start loading sequence
        self.root.after(100, self.load_step_1)
        
        # Start the event loop
        self.root.mainloop()
    
    def load_step_1(self):
        """Load pandas"""
        try:
            self.splash.update_status("From Billings to Bozeman...")
            import pandas as pd
            self.pd = pd
            self.root.after(100, self.load_step_2)
        except Exception as e:
            self.show_error(f"Error loading pandas: {str(e)}")
    
    def load_step_2(self):
        """Load geopandas"""
        try:
            self.splash.update_status("Spanning the Big Sky Country...")
            import geopandas as gpd
            self.gpd = gpd
            self.root.after(100, self.load_step_3)
        except Exception as e:
            self.show_error(f"Error loading geopandas: {str(e)}")
    
    def load_step_3(self):
        """Load matplotlib"""
        try:
            self.splash.update_status("Mapping Montana's vast landscapes...")
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            self.plt = plt
            self.FigureCanvasTkAgg = FigureCanvasTkAgg
            self.root.after(100, self.load_step_4)
        except Exception as e:
            self.show_error(f"Error loading matplotlib: {str(e)}")
    
    def load_step_4(self):
        """Load shapefile"""
        try:
            self.splash.update_status("Connecting all 56 counties...")
            self.load_shapefile()
            self.root.after(100, self.show_selection_screen)
        except Exception as e:
            self.show_error(f"Error loading shapefile: {str(e)}")
    
    def show_selection_screen(self):
        """Show the selection screen after splash"""
        self.splash.destroy()
        SelectionScreen(self.root, self)
    
    def show_error(self, message):
        """Show error message and exit"""
        if hasattr(self, 'splash'):
            self.splash.destroy()
        messagebox.showerror("Error", message)
        self.root.quit()
        sys.exit(1)
    
    def load_shapefile(self):
        try:
            # Get the application's base directory
            if getattr(sys, 'frozen', False):
                base_dir = sys._MEIPASS
            else:
                base_dir = os.path.dirname(os.path.abspath(__file__))
            
            shapefile_path = os.path.join(base_dir, "MontanaCounties_shp", "County.shp")
            
            if not os.path.exists(shapefile_path):
                raise FileNotFoundError(
                    f"Shapefile not found at:\n{shapefile_path}\n\n"
                    "Please ensure the MontanaCounties_shp folder is in the correct location."
                )
                
            self.gdf = self.gpd.read_file(shapefile_path)
            self.gdf.columns = self.gdf.columns.str.strip()
            self.gdf["County"] = self.gdf["NAME"].str.strip().str.lower()
            self.gdf["Color"] = "white"
            
        except Exception as e:
            raise Exception(f"Error loading shapefile:\n{str(e)}\n\nPlease ensure the shapefile is not corrupted and try again.")
    
    def load_excel(self):
        path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
        if not path:
            return
        
        # Show loading dialog
        loading_window = tk.Tk()
        loading_window.title("Loading")
        screen_width = loading_window.winfo_screenwidth()
        screen_height = loading_window.winfo_screenheight()
        window_width = 300
        window_height = 100
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        loading_window.geometry(f"{window_width}x{window_height}+{x}+{y}")
        loading_window.overrideredirect(True)
        
        frame = ttk.Frame(loading_window, padding="20", relief="raised")
        frame.pack(fill='both', expand=True)
        
        loading_label = ttk.Label(frame, text="Loading file...\nPlease wait", font=('Helvetica', 10))
        loading_label.pack(pady=10)
        
        progress = ttk.Progressbar(frame, mode='indeterminate')
        progress.pack(fill='x', pady=5)
        progress.start(10)
        
        loading_window.update()
        
        try:
            # Load the Excel file
            self.df = self.pd.read_excel(path, sheet_name=0)
            
            # # Print all column names immediately after loading
            # print("\nAvailable columns in Excel file:")
            # print("--------------------------------")
            # for i, col in enumerate(self.df.columns, 1):
            #     print(f"{i}. {col}")
            # print("--------------------------------\n")
            
            # Now continue with column processing
            self.df.columns = self.df.columns.str.strip()
            
            # Get just the filename from the path
            filename = path.split('/')[-1]
            
            # Validate required columns
            required_columns = ["county", "family", "genus", "species", "year"]
            missing_columns = [col for col in required_columns if col not in self.df.columns]
            
            if missing_columns:
                progress.stop()
                loading_window.destroy()
                self.selected_file_var.set("No file selected")
                messagebox.showerror("Error", 
                    f"Missing required columns: {', '.join(missing_columns)}\n\n"
                    "The following columns are required:\n"
                    "- county: for mapping locations\n"
                    "- family: for taxonomic classification\n"
                    "- genus: for taxonomic classification\n"
                    "- species: for taxonomic classification\n"
                    "- year: for temporal analysis\n\n"
                    "Please check your Excel file and try again."
                )
                return
            
            # Process the data
            # First standardize county names
            self.df["county"] = self.standardize_county_names(self.df["county"])
            
            # Process other columns
            for col in ["family", "genus", "species"]:
                self.df[col] = self.df[col].astype(str).str.strip().str.lower()
            self.df["year"] = self.df["year"].astype(str).str.extract(r'(\d{4})').astype(float)
            
            # Get valid Montana counties from shapefile
            valid_counties = set(self.standardize_county_names(self.gdf["County"]))
            
            # Filter DataFrame to only include valid Montana counties
            montana_records = self.df[self.df["county"].isin(valid_counties)]
            
            if len(montana_records) == 0:
                progress.stop()
                loading_window.destroy()
                self.selected_file_var.set("No file selected")
                messagebox.showerror("Error", 
                    "No valid Montana county records found in the Excel file.\n\n"
                    "Please check that your data contains Montana county records."
                )
                return
            
            # Replace the main DataFrame with only Montana records
            self.df = montana_records
            
            # Calculate statistics using Montana records
            num_records = len(montana_records)
            num_families = len(montana_records["family"].unique())
            num_genera = len(montana_records["genus"].unique())
            num_species = len(montana_records["species"].unique())
            year_range = f"{int(montana_records['year'].min())} - {int(montana_records['year'].max())}"
            num_counties = len(montana_records["county"].unique())
            
            # Get valid families (non-empty/non-null values)
            valid_families = sorted(montana_records["family"].dropna().unique())
            valid_families = [f for f in valid_families if str(f).strip() and str(f).lower() != 'nan']  # Remove empty strings and 'nan'
            
            # Capitalize family names
            family_values = ["All"] + [f.title() for f in valid_families]
            
            # Update Family dropdown
            self.family_dropdown["values"] = family_values
            self.family_dropdown.set("Select Family")
            
            # Reset other dropdowns and fields
            self.genus_dropdown.set("Select Genus")
            self.genus_dropdown["values"] = []
            self.species_dropdown.set("Select Species")
            self.species_dropdown["values"] = []
            self.year_var.set("")
            
            # Update file info display
            self.selected_file_var.set(f"✓ {filename}\n{num_records:,} Montana records loaded")
            
            # Stop progress bar and close loading window
            progress.stop()
            loading_window.destroy()
            
            # Show success message with detailed statistics
            messagebox.showinfo("Success", 
                f"File loaded successfully!\n\n"
                f"Montana Dataset Summary:\n"
                f"• Total Records: {num_records:,}\n"
                f"• Unique Families: {num_families}\n"
                f"• Unique Genera: {num_genera}\n"
                f"• Unique Species: {num_species}\n"
                f"• Counties Covered: {num_counties}\n"
                f"• Year Range: {year_range}\n\n"
                "Please select a Family to continue."
            )
            
            print("✅ Excel file loaded successfully!")
            
        except Exception as e:
            progress.stop()
            loading_window.destroy()
            self.selected_file_var.set("No file selected")
            error_message = str(e)
            if "No sheet named" in error_message:
                error_message = "Invalid Excel file format. Please ensure your data is in the first sheet."
            elif "Invalid file" in error_message:
                error_message = "Invalid file format. Please ensure you're uploading a valid Excel (.xlsx) file."
            
            messagebox.showerror("Error", 
                f"Error loading file:\n{error_message}\n\n"
                "Please check your Excel file format and try again."
            )
    
    def update_genus_dropdown(self, event=None):
        family = self.selected_family.get().strip()
        
        if family == "Select Family":
            self.genus_dropdown["values"] = []
            self.genus_dropdown.set("Select Genus")
            return
        
        # Filter based on family selection
        if family == "All":
            # Get all non-empty genus values
            filtered = self.df[self.df["genus"].notna() & (self.df["genus"].str.strip() != "")]
        else:
            # Get genus for specific family (case-insensitive)
            filtered = self.df[self.df["family"].str.lower() == family.lower()]
        
        # Get valid genera (non-empty/non-null values)
        valid_genera = sorted(filtered["genus"].dropna().unique())
        valid_genera = [g for g in valid_genera if str(g).strip() and str(g).lower() != 'nan']  # Remove empty strings and 'nan'
        
        # Create genus list with special options
        genus_values = ["All"] + [g.title() for g in valid_genera]
        
        # Update Genus dropdown
        self.genus_dropdown["values"] = genus_values
        self.genus_dropdown.set("Select Genus")
        
        # Reset species dropdown
        self.species_dropdown.set("Select Species")
        self.species_dropdown["values"] = []
    
    def update_species_dropdown(self, event=None):
        family = self.selected_family.get().strip()
        genus = self.selected_genus.get().strip()
        
        if family == "Select Family" or genus == "Select Genus":
            self.species_dropdown["values"] = []
            self.species_dropdown.set("Select Species")
            return
        
        # Start with base DataFrame
        filtered = self.df
        
        # Apply family filter
        if family == "All":
            filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
        else:
            filtered = filtered[filtered["family"].str.lower() == family.lower()]
        
        # Apply genus filter
        if genus == "All":
            filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
        else:
            filtered = filtered[filtered["genus"].str.lower() == genus.lower()]
        
        # Get valid species (non-empty/non-null values)
        valid_species = sorted(filtered["species"].dropna().unique())
        valid_species = [s for s in valid_species if str(s).strip() and str(s).lower() != 'nan']  # Remove empty strings and 'nan'
        
        # Create species list with special options - note lowercase for species
        species_values = ["all"] + valid_species
        
        # Update Species dropdown
        self.species_dropdown["values"] = species_values
        self.species_dropdown.set("Select Species")
    
    def is_valid_color(self, color):
        """Validate if a color string is a valid matplotlib color"""
        try:
            # Convert color to RGB
            self.plt.matplotlib.colors.to_rgb(color)
            return True
        except ValueError:
            return False

    def validate_colors(self):
        """Validate all selected colors"""
        colors = {
            'Pre-Year': self.pre_color.get(),
            'Post-Year': self.post_color.get(),
            'All Records': self.all_color.get()
        }
        
        invalid_colors = []
        for name, color in colors.items():
            if not self.is_valid_color(color):
                invalid_colors.append(f"{name} color '{color}'")
        
        if invalid_colors:
            error_msg = "Invalid colors detected:\n" + "\n".join(invalid_colors)
            self.toast.show_toast(error_msg, duration=5000, error=True)
            return False
        return True

    def on_color_change(self, event=None):
        """Validate colors when they change"""
        self.validate_colors()

    def generate_map(self):
        if self.map_canvas:
            self.map_canvas.get_tk_widget().destroy()
        
        # Validate colors first
        if not self.validate_colors():
            self.download_button.config(state="disabled")
            return
        
        gdf_copy = self.gdf.copy()
        # Initialize all counties as white
        gdf_copy["Color"] = "white"
        
        fam = self.selected_family.get().strip()
        gen = self.selected_genus.get().strip()
        spec = self.selected_species.get().strip()
        year = self.year_var.get().strip()
        
        if not fam or fam == "Select Family" or not gen or gen == "Select Genus" or not spec or spec == "Select Species":
            messagebox.showerror("Missing Input", "Please select Family, Genus, and Species.")
            self.download_button.config(state="disabled")
            return
        
        # Start with base DataFrame
        filtered = self.df
        
        # Apply family filter
        if fam == "All":
            filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
        elif fam == "Not Specified":
            filtered = filtered[filtered["family"].isna() | (filtered["family"].str.strip() == "")]
        else:
            filtered = filtered[filtered["family"].str.lower() == fam.lower()]
        
        # Apply genus filter
        if gen == "All":
            filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
        elif gen == "Not Specified":
            filtered = filtered[filtered["genus"].isna() | (filtered["genus"].str.strip() == "")]
        else:
            filtered = filtered[filtered["genus"].str.lower() == gen.lower()]
        
        # Apply species filter
        if spec == "all":
            filtered = filtered[filtered["species"].notna() & (filtered["species"].str.strip() != "")]
        elif spec == "not specified":
            filtered = filtered[filtered["species"].isna() | (filtered["species"].str.strip() == "")]
        else:
            filtered = filtered[filtered["species"].str.lower() == spec.lower()]
        
        # Create a set of valid county names from the shapefile for quick lookup
        valid_counties = set(self.standardize_county_names(gdf_copy["County"]))
        
        # Track unmatched counties to report to user
        unmatched_counties = set()
        
        if isinstance(year, str) and year.isdigit():
            year = int(year)
            # First pass: Mark counties with post-dividing year records with post_color (only if they're white)
            for county in filtered[filtered["year"] > year]["county"].unique():
                county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
                if county_lower in valid_counties:
                    mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                    if gdf_copy.loc[mask, "Color"].iloc[0] == "white":
                        gdf_copy.loc[mask, "Color"] = self.post_color.get()
                else:
                    unmatched_counties.add(county)
            
            # Second pass: Mark counties with pre-dividing year records with pre_color (overriding any color)
            for county in filtered[filtered["year"] <= year]["county"].unique():
                county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
                if county_lower in valid_counties:
                    mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                    gdf_copy.loc[mask, "Color"] = self.pre_color.get()
                else:
                    unmatched_counties.add(county)
        else:
            # If no year specified, mark all counties with records using all_color
            for county in filtered["county"].unique():
                county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
                if county_lower in valid_counties:
                    mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                    gdf_copy.loc[mask, "Color"] = self.all_color.get()
                else:
                    unmatched_counties.add(county)
        
        # Report any unmatched counties
        if unmatched_counties:
            print("\nWarning: The following counties in your Excel file don't match the shapefile counties:")
            print("-------------------------------------------------------------------------")
            for county in sorted(unmatched_counties):
                print(f"• {county}")
            print("\nValid Montana county names:")
            print("--------------------------------")
            for county in sorted(valid_counties):
                print(f"• {county}")
            print("-------------------------------------------------------------------------\n")
            
            messagebox.showwarning("County Name Mismatch",
                f"Some counties in your Excel file don't match the shapefile counties.\n\n"
                f"Number of unmatched counties: {len(unmatched_counties)}\n\n"
                f"Please check the console output for details and ensure county names match exactly."
            )
        
        # Create figure with calculated size
        fig = self.plt.figure(figsize=(12, 11))
        
        # Create main map axis with sufficient space for title and legend
        ax = fig.add_axes([0.1, 0.2, 0.8, 0.6])  # Adjusted to leave more space at top
        gdf_copy.boundary.plot(ax=ax, linewidth=1, edgecolor="black")
        gdf_copy.plot(ax=ax, color=gdf_copy["Color"], alpha=0.6)
        
        # Add title with dynamic font size
        title = f"{fam.title()} > {gen.title()} > {spec.lower()}"
        if isinstance(year, int):
            subtitle = f"\nYear: {year}"
            title += subtitle
        
        # Calculate dynamic font size based on figure width
        title_fontsize = min(15, max(8, fig.get_figwidth() * 1.5))
        ax.set_title(title, fontsize=title_fontsize, pad=25, wrap=True)  # Increased padding
        ax.axis("off")
        
        # Create a separate axis for the legend
        legend_ax = fig.add_axes([0.2, 0.02, 0.6, 0.12])
        legend_ax.axis('off')
        
        # Add a box around the legend
        legend_box = self.plt.Rectangle((0, 0), 1, 1, 
                                      facecolor='white', edgecolor='black',
                                      transform=legend_ax.transAxes)
        legend_ax.add_patch(legend_box)
        
        if isinstance(year, int):
            # Add horizontal bars with their descriptions for year analysis
            bar_length = 0.15
            bar_height = 0.1
            
            # Pre-year color
            legend_ax.add_patch(
                self.plt.Rectangle((0.2, 0.5), bar_length, bar_height, 
                                 facecolor=self.pre_color.get(), 
                                 alpha=0.6,
                                 edgecolor='black')
            )
            legend_ax.text(0.4, 0.55, f"Records ≤ {year}", fontsize=10, va='center')
            
            # Post-year color
            legend_ax.add_patch(
                self.plt.Rectangle((0.2, 0.2), bar_length, bar_height, 
                                 facecolor=self.post_color.get(), 
                                 alpha=0.6,
                                 edgecolor='black')
            )
            legend_ax.text(0.4, 0.25, f"Records > {year}", fontsize=10, va='center')
        else:
            # Add single bar for all records
            bar_length = 0.15
            bar_height = 0.1
            legend_ax.add_patch(
                self.plt.Rectangle((0.2, 0.35), bar_length, bar_height, 
                                 facecolor=self.all_color.get(), 
                                 alpha=0.6,
                                 edgecolor='black')
            )
            legend_ax.text(0.4, 0.4, "All Records", fontsize=10, va='center')
        
        # Set the legend axis limits
        legend_ax.set_xlim(0, 1)
        legend_ax.set_ylim(0, 1)
        
        # Adjust main plot to make room for legend
        fig.subplots_adjust(bottom=0.2, top=0.9)
        
        self.map_canvas = self.FigureCanvasTkAgg(fig, master=self.right_panel)
        self.map_canvas.draw()
        self.map_canvas.get_tk_widget().pack(fill='both', expand=True)
        
        self.current_fig = fig
        self.download_button.config(state="normal")
        print("✅ Map generated successfully!")
    
    def download_map(self):
        if self.current_fig:
            try:
                # Get Downloads folder path
                downloads_path = str(Path.home() / "Downloads")
                
                # Get taxonomic info and year
                fam = self.selected_family.get().strip().title()
                gen = self.selected_genus.get().strip().title()
                spec = self.selected_species.get().strip().lower()
                year = self.year_var.get().strip()
                
                # Create timestamp with current date (YYYYMMDD_HHMM)
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
                
                # Create filename with taxonomic info
                year_info = f"_{year}" if year else ""
                filename = f"{fam}-{gen}-{spec}{year_info}_{timestamp}.tiff"
                
                # Full path for the file
                file_path = os.path.join(downloads_path, filename)
                
                # Save the figure
                self.current_fig.savefig(file_path, format="tiff", dpi=300)
                
                # Show toast notification
                self.toast.show_toast(f"Map saved as {filename}")
                
                print(f"✅ TIFF map saved as '{file_path}'")
                
            except Exception as e:
                messagebox.showerror("Error", 
                    f"Error saving file:\n{str(e)}\n\n"
                    "Please try again."
                )
    
    def initialize_gui(self):
        # Configure style
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TLabel', font=('Helvetica', 10))
        style.configure('TButton', font=('Helvetica', 10))
        style.configure('Header.TLabel', font=('Helvetica', 12, 'bold'))
        style.configure('Title.TLabel', font=('Helvetica', 24, 'bold'))  # Increased from 16 to 24
        style.configure('Back.TButton', 
                       font=('Helvetica', 10, 'bold'),
                       padding=8)
        
        # Create main container with padding
        main_container = ttk.Frame(self.root, padding="20")
        main_container.pack(fill='both', expand=True)
        
        # Title
        title_frame = ttk.Frame(main_container)
        title_frame.pack(fill='x', pady=(0, 20))
        title_label = ttk.Label(
            title_frame, 
            text="Montana Bee County Map Generator - Single Year Analysis", 
            style='Title.TLabel',
            foreground='dark green'
        )
        title_label.pack()
        
        # Calculate left and right panel widths based on screen size
        screen_width = self.root.winfo_screenwidth()
        left_panel_width = min(300, int(screen_width * 0.2))  # 20% of screen width or 300px, whichever is smaller
        
        # Create left panel for controls with scrollbar
        left_panel_container = ttk.Frame(main_container, width=left_panel_width)
        left_panel_container.pack(side='left', fill='y', padx=(0, 20))
        left_panel_container.pack_propagate(False)  # Prevent the frame from shrinking
        
        # Back button at the top of left panel container (outside scroll area)
        back_button = ttk.Button(
            left_panel_container,
            text="← Back to Selection",
            command=self.go_back,
            style='Back.TButton'
        )
        back_button.pack(fill='x', pady=(0, 10))
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(left_panel_container)
        scrollbar.pack(side='right', fill='y')
        
        # Create canvas for scrollable content with matching background
        bg_color = style.lookup('TFrame', 'background')  # Get ttk frame background color
        if not bg_color:  # Fallback if style lookup fails
            bg_color = self.left_panel_container.cget('background')
        self.scroll_canvas = Canvas(
            left_panel_container,
            yscrollcommand=scrollbar.set,
            width=left_panel_width-20,
            bg=bg_color,
            highlightthickness=0  # Remove border
        )
        self.scroll_canvas.pack(side='left', fill='y')
        
        scrollbar.config(command=self.scroll_canvas.yview)
        
        # Create frame for controls inside canvas
        left_panel = ttk.Frame(self.scroll_canvas)
        self.scroll_canvas.create_window((0, 0), window=left_panel, anchor='nw', width=left_panel_width-20)
        
        # Excel Load Section
        excel_frame = ttk.LabelFrame(left_panel, text="Data Input", padding="10")
        excel_frame.pack(fill='x', pady=(0, 20))
        
        # Create a frame for the button and file info
        file_info_frame = ttk.Frame(excel_frame)
        file_info_frame.pack(fill='x', expand=True)
        
        # Load button
        load_button = ttk.Button(
            file_info_frame, 
            text="Load Excel File", 
            command=self.load_excel, 
            style='TButton'
        )
        load_button.pack(fill='x', pady=(0, 5))
        
        # File info label
        file_label = ttk.Label(
            file_info_frame, 
            textvariable=self.selected_file_var,
            style='FileInfo.TLabel',
            wraplength=250
        )
        file_label.pack(fill='x')
        
        # Configure style for file info
        style.configure('FileInfo.TLabel', 
                       font=('Helvetica', 9),
                       foreground='dark green')
        
        # Color Selection Section
        color_frame = ttk.LabelFrame(left_panel, text="Color Settings", padding="10")
        color_frame.pack(fill='x', pady=(0, 20))
        
        # Common color options
        color_options = ["red", "blue", "grey", "black", "yellow", "purple", "orange", "pink", "brown"]
        
        # Pre Color
        ttk.Label(color_frame, text="Pre-Year Color:", style='TLabel').pack(fill='x')
        pre_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.pre_color, 
            values=color_options,
            state="normal"
        )
        pre_color_combo.pack(fill='x', pady=(0, 10))
        pre_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        pre_color_combo.bind('<Return>', self.on_color_change)
        pre_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Post Color
        ttk.Label(color_frame, text="Post-Year Color:", style='TLabel').pack(fill='x')
        post_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.post_color, 
            values=color_options,
            state="normal"
        )
        post_color_combo.pack(fill='x', pady=(0, 10))
        post_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        post_color_combo.bind('<Return>', self.on_color_change)
        post_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # All Color
        ttk.Label(color_frame, text="All Records Color:", style='TLabel').pack(fill='x')
        all_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.all_color, 
            values=color_options,
            state="normal"
        )
        all_color_combo.pack(fill='x', pady=(0, 10))
        all_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        all_color_combo.bind('<Return>', self.on_color_change)
        all_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Add helper text for custom colors
        helper_label = ttk.Label(
            color_frame, 
            text="Tip: You can type any valid color name or hex code (e.g., #FF5733)",
            style='TLabel',
            wraplength=250
        )
        helper_label.pack(fill='x', pady=(5, 0))
        
        # Year Input Section
        year_frame = ttk.LabelFrame(left_panel, text="Year Filter", padding="10")
        year_frame.pack(fill='x', pady=(0, 20))
        ttk.Label(year_frame, text="Year (Optional):", style='TLabel').pack(fill='x')
        ttk.Entry(year_frame, textvariable=self.year_var).pack(fill='x')
        
        # Species Selection Section
        species_frame = ttk.LabelFrame(left_panel, text="Species Selection", padding="10")
        species_frame.pack(fill='x', pady=(0, 20))
        
        # Family
        ttk.Label(species_frame, text="Family:", style='TLabel').pack(fill='x')
        self.family_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_family, state="readonly")
        self.family_dropdown.pack(fill='x', pady=(0, 10))
        
        # Genus
        ttk.Label(species_frame, text="Genus:", style='TLabel').pack(fill='x')
        self.genus_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_genus, state="readonly")
        self.genus_dropdown.pack(fill='x', pady=(0, 10))
        
        # Species
        ttk.Label(species_frame, text="Species:", style='TLabel').pack(fill='x')
        self.species_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_species, state="readonly")
        self.species_dropdown.pack(fill='x', pady=(0, 10))
        
        # Button Section
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill='x', pady=(0, 20))
        
        # Style configuration for buttons
        style.configure('Action.TButton', 
                       font=('Helvetica', 10, 'bold'),
                       padding=10)
        
        generate_button = ttk.Button(
            button_frame, 
            text="Generate Map", 
            command=self.generate_map, 
            style='Action.TButton'
        )
        generate_button.pack(fill='x', pady=(0, 5))
        
        self.download_button = ttk.Button(
            button_frame, 
            text="Download Map", 
            command=self.download_map, 
            state="disabled", 
            style='Action.TButton'
        )
        self.download_button.pack(fill='x', pady=(0, 5))
        
        # Create right panel for map display with dynamic width
        self.right_panel = ttk.Frame(main_container)
        self.right_panel.pack(side='left', fill='both', expand=True)
        
        # Function to update panel sizes on window resize
        def update_panel_sizes(event=None):
            if event and event.widget != self.root:
                return
                
            current_width = self.root.winfo_width()
            new_left_width = min(300, int(current_width * 0.2))
            
            # Update left panel width
            left_panel_container.configure(width=new_left_width)
            self.scroll_canvas.configure(width=new_left_width-20)
            self.scroll_canvas.itemconfig(self.scroll_canvas.find_all()[0], width=new_left_width-20)
            
            # Right panel will automatically adjust due to expand=True
            
            # If there's a map canvas, trigger a redraw
            if hasattr(self, 'map_canvas') and self.map_canvas:
                self.map_canvas.draw()
        
        # Bind resize event
        self.root.bind('<Configure>', update_panel_sizes)
        
        # Configure scrolling
        def on_configure(event):
            self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox('all'))
        
        left_panel.bind('<Configure>', on_configure)
        
        # Enable mousewheel scrolling
        def on_mousewheel(event):
            self.scroll_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        
        self.scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        # Bind dropdowns
        self.family_dropdown.bind("<<ComboboxSelected>>", self.update_genus_dropdown)
        self.genus_dropdown.bind("<<ComboboxSelected>>", self.update_species_dropdown)
    
    def go_back(self):
        """Return to the selection screen"""
        # Store current position and state
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        current_state = self.root.state()
        self.root.destroy()
        selection = SelectionScreen(self.root.master, self.main_app, from_analysis=True)
        # Set position first
        selection.root.geometry(f"{selection.root.winfo_screenwidth()}x{selection.root.winfo_screenheight()}+{x}+{y}")
        selection.root.update()  # Force geometry update
        # Restore the previous window state
        if current_state == 'zoomed':
            selection.root.state('zoomed')
    
    def on_window_resize(self, event=None):
        try:
            # Get current window dimensions
            current_width = self.root.winfo_width()
            current_height = self.root.winfo_height()
            
            # Only process if we have valid dimensions
            if current_width <= 1 or current_height <= 1:
                return
                
            # Store current dimensions for restoring from maximize
            if not hasattr(self, 'last_valid_size'):
                self.last_valid_size = {
                    'width': current_width,
                    'height': current_height,
                    'x': self.root.winfo_x(),
                    'y': self.root.winfo_y()
                }
            
            # Update panels if needed
            self.update_panel_sizes()
                
        except Exception as e:
            print(f"Warning: Resize handling error: {str(e)}")  # For debugging
            pass  # Silently handle any errors during resize

    def update_panel_sizes(self):
        try:
            current_width = self.root.winfo_width()
            
            # Calculate new left panel width (20% of window width, max 300px)
            new_left_width = min(300, max(200, int(current_width * 0.2)))
            
            # Update left panel container if it exists
            if hasattr(self, 'left_panel_container'):
                self.left_panel_container.configure(width=new_left_width)
            
            # Update scroll canvas if it exists and has content
            if hasattr(self, 'scroll_canvas'):
                try:
                    # Only update if canvas exists and has content
                    if self.scroll_canvas.find_all():
                        self.scroll_canvas.configure(width=new_left_width-20)
                        self.scroll_canvas.itemconfig(
                            self.scroll_canvas.find_all()[0],
                            width=new_left_width-20
                        )
                except tk.TclError:
                    pass  # Handle case where canvas is being destroyed
            
            # Redraw map if it exists
            if hasattr(self, 'map_canvas') and self.map_canvas:
                try:
                    self.map_canvas.draw()
                except Exception:
                    pass  # Handle any drawing errors silently
                    
            # Force geometry update
            self.root.update_idletasks()
                
        except Exception as e:
            print(f"Warning: Panel update error: {str(e)}")  # For debugging
            pass  # Silently handle any errors

class SingleYearAnalysis:
    def __init__(self, parent, main_app):
        self.root = tk.Toplevel(parent)
        self.root.title("Single Year Analysis - Montana Bee County Map Generator")
        
        # Allow the window to be moved to any screen
        self.root.attributes('-alpha', 1.0)
        self.root.attributes('-topmost', False)
        
        # Store main_app reference
        self.main_app = main_app
        
        # Get dependencies from main_app
        self.pd = main_app.pd
        self.plt = main_app.plt
        self.FigureCanvasTkAgg = main_app.FigureCanvasTkAgg
        
        # Set window icon
        if getattr(sys, 'frozen', False):
            base_dir = sys._MEIPASS
        else:
            base_dir = os.path.dirname(os.path.abspath(__file__))
        icon_path = os.path.join(base_dir, "app_icon.ico")
        if os.path.exists(icon_path):
            self.root.iconbitmap(icon_path)
        
        # Get screen dimensions
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        
        # Set initial size but don't maximize yet - that will be done by the calling method
        self.root.geometry(f"{screen_width}x{screen_height}")
        
        # Ensure window can be moved and resized
        self.root.resizable(True, True)
        self.root.minsize(800, 600)
        
        # Initialize variables
        self.map_canvas = None
        self.current_fig = None
        
        # Initialize StringVar variables
        self.pre_color = StringVar(self.root)
        self.post_color = StringVar(self.root)
        self.all_color = StringVar(self.root)
        self.year_var = StringVar(self.root)
        self.selected_family = StringVar(self.root)
        self.selected_genus = StringVar(self.root)
        self.selected_species = StringVar(self.root)
        self.selected_file_var = StringVar(self.root)
        
        # Set default values
        self.pre_color.set("grey")
        self.post_color.set("red")
        self.all_color.set("yellow")
        self.selected_file_var.set("No file selected")
        
        # Create toast notification instance
        self.toast = ToastNotification(self.root)
        
        # Initialize pandas DataFrame
        self.df = self.pd.DataFrame()
        
        # Get the shapefile data from parent
        self.gdf = main_app.gdf.copy()
        
        # Initialize GUI
        self.initialize_gui()
        
        # Bind window state change
        self.root.bind("<Configure>", self.on_window_resize)

    def standardize_county_names(self, county_series):
        """
        Standardize county names by:
        1. Converting '&' to 'and'
        2. Stripping whitespace
        3. Converting to lowercase
        """
        return county_series.str.strip().str.lower().str.replace('&', 'and')

    def load_excel(self):
        path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
        if not path:
            return
        
        # Show loading dialog
        loading_window = tk.Tk()
        loading_window.title("Loading")
        screen_width = loading_window.winfo_screenwidth()
        screen_height = loading_window.winfo_screenheight()
        window_width = 300
        window_height = 100
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        loading_window.geometry(f"{window_width}x{window_height}+{x}+{y}")
        loading_window.overrideredirect(True)
        
        frame = ttk.Frame(loading_window, padding="20", relief="raised")
        frame.pack(fill='both', expand=True)
        
        loading_label = ttk.Label(frame, text="Loading file...\nPlease wait", font=('Helvetica', 10))
        loading_label.pack(pady=10)
        
        progress = ttk.Progressbar(frame, mode='indeterminate')
        progress.pack(fill='x', pady=5)
        progress.start(10)
        
        loading_window.update()
        
        try:
            # Load the Excel file
            self.df = self.pd.read_excel(path, sheet_name=0)
            
            # # Print all column names immediately after loading
            # print("\nAvailable columns in Excel file:")
            # print("--------------------------------")
            # for i, col in enumerate(self.df.columns, 1):
            #     print(f"{i}. {col}")
            # print("--------------------------------\n")
            
            # Now continue with column processing
            self.df.columns = self.df.columns.str.strip()
            
            # Get just the filename from the path
            filename = path.split('/')[-1]
            
            # Validate required columns
            required_columns = ["county", "family", "genus", "species", "year"]
            missing_columns = [col for col in required_columns if col not in self.df.columns]
            
            if missing_columns:
                progress.stop()
                loading_window.destroy()
                self.selected_file_var.set("No file selected")
                messagebox.showerror("Error", 
                    f"Missing required columns: {', '.join(missing_columns)}\n\n"
                    "The following columns are required:\n"
                    "- county: for mapping locations\n"
                    "- family: for taxonomic classification\n"
                    "- genus: for taxonomic classification\n"
                    "- species: for taxonomic classification\n"
                    "- year: for temporal analysis\n\n"
                    "Please check your Excel file and try again."
                )
                return
            
            # Process the data
            # First standardize county names
            self.df["county"] = self.standardize_county_names(self.df["county"])
            
            # Process other columns
            for col in ["family", "genus", "species"]:
                self.df[col] = self.df[col].astype(str).str.strip().str.lower()
            self.df["year"] = self.df["year"].astype(str).str.extract(r'(\d{4})').astype(float)
            
            # Get valid Montana counties from shapefile
            valid_counties = set(self.standardize_county_names(self.gdf["County"]))
            
            # Filter DataFrame to only include valid Montana counties
            montana_records = self.df[self.df["county"].isin(valid_counties)]
            
            if len(montana_records) == 0:
                progress.stop()
                loading_window.destroy()
                self.selected_file_var.set("No file selected")
                messagebox.showerror("Error", 
                    "No valid Montana county records found in the Excel file.\n\n"
                    "Please check that your data contains Montana county records."
                )
                return
            
            # Replace the main DataFrame with only Montana records
            self.df = montana_records
            
            # Calculate statistics using Montana records
            num_records = len(montana_records)
            num_families = len(montana_records["family"].unique())
            num_genera = len(montana_records["genus"].unique())
            num_species = len(montana_records["species"].unique())
            year_range = f"{int(montana_records['year'].min())} - {int(montana_records['year'].max())}"
            num_counties = len(montana_records["county"].unique())
            
            # Get valid families (non-empty/non-null values)
            valid_families = sorted(montana_records["family"].dropna().unique())
            valid_families = [f for f in valid_families if str(f).strip() and str(f).lower() != 'nan']  # Remove empty strings and 'nan'
            
            # Capitalize family names
            family_values = ["All"] + [f.title() for f in valid_families]
            
            # Update Family dropdown
            self.family_dropdown["values"] = family_values
            self.family_dropdown.set("All")
            
            # Trigger genus dropdown update
            self.update_genus_dropdown()
            
            # Set genus to "All" and trigger species update
            self.genus_dropdown.set("All")
            self.update_species_dropdown()
            
            # Set species to "all"
            self.species_dropdown.set("all")
            
            # Reset year field
            self.year_var.set("")
            
            # Update file info display
            self.selected_file_var.set(f"✓ {filename}\n{num_records:,} Montana records loaded")
            
            # Stop progress bar and close loading window
            progress.stop()
            loading_window.destroy()
            
            # Show success message with detailed statistics
            messagebox.showinfo("Success", 
                f"File loaded successfully!\n\n"
                f"Montana Dataset Summary:\n"
                f"• Total Records: {num_records:,}\n"
                f"• Unique Families: {num_families}\n"
                f"• Unique Genera: {num_genera}\n"
                f"• Unique Species: {num_species}\n"
                f"• Counties Covered: {num_counties}\n"
                f"• Year Range: {year_range}\n\n"
                "Please select a Family to continue."
            )
            
            print("✅ Excel file loaded successfully!")
            
        except Exception as e:
            progress.stop()
            loading_window.destroy()
            self.selected_file_var.set("No file selected")
            error_message = str(e)
            if "No sheet named" in error_message:
                error_message = "Invalid Excel file format. Please ensure your data is in the first sheet."
            elif "Invalid file" in error_message:
                error_message = "Invalid file format. Please ensure you're uploading a valid Excel (.xlsx) file."
            
            messagebox.showerror("Error", 
                f"Error loading file:\n{error_message}\n\n"
                "Please check your Excel file format and try again."
            )
    
    def update_genus_dropdown(self, event=None):
        family = self.selected_family.get().strip()
        
        if family == "Select Family":
            self.genus_dropdown["values"] = []
            self.genus_dropdown.set("Select Genus")
            return
        
        # Filter based on family selection
        if family == "All":
            # Get all non-empty genus values
            filtered = self.df[self.df["genus"].notna() & (self.df["genus"].str.strip() != "")]
        else:
            # Get genus for specific family (case-insensitive)
            filtered = self.df[self.df["family"].str.lower() == family.lower()]
        
        # Get valid genera (non-empty/non-null values)
        valid_genera = sorted(filtered["genus"].dropna().unique())
        valid_genera = [g for g in valid_genera if str(g).strip() and str(g).lower() != 'nan']  # Remove empty strings and 'nan'
        
        # Create genus list with special options
        genus_values = ["All"] + [g.title() for g in valid_genera]
        
        # Update Genus dropdown
        self.genus_dropdown["values"] = genus_values
        self.genus_dropdown.set("Select Genus")
        
        # Reset species dropdown
        self.species_dropdown.set("Select Species")
        self.species_dropdown["values"] = []
    
    def update_species_dropdown(self, event=None):
        family = self.selected_family.get().strip()
        genus = self.selected_genus.get().strip()
        
        if family == "Select Family" or genus == "Select Genus":
            self.species_dropdown["values"] = []
            self.species_dropdown.set("Select Species")
            return
        
        # Start with base DataFrame
        filtered = self.df
        
        # Apply family filter
        if family == "All":
            filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
        else:
            filtered = filtered[filtered["family"].str.lower() == family.lower()]
        
        # Apply genus filter
        if genus == "All":
            filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
        else:
            filtered = filtered[filtered["genus"].str.lower() == genus.lower()]
        
        # Get valid species (non-empty/non-null values)
        valid_species = sorted(filtered["species"].dropna().unique())
        valid_species = [s for s in valid_species if str(s).strip() and str(s).lower() != 'nan']  # Remove empty strings and 'nan'
        
        # Create species list with special options - note lowercase for species
        species_values = ["all"] + valid_species
        
        # Update Species dropdown
        self.species_dropdown["values"] = species_values
        self.species_dropdown.set("Select Species")
    
    def is_valid_color(self, color):
        """Validate if a color string is a valid matplotlib color"""
        try:
            # Convert color to RGB
            self.plt.matplotlib.colors.to_rgb(color)
            return True
        except ValueError:
            return False

    def validate_colors(self):
        """Validate all selected colors"""
        colors = {
            'Pre-Year': self.pre_color.get(),
            'Post-Year': self.post_color.get(),
            'All Records': self.all_color.get()
        }
        
        invalid_colors = []
        for name, color in colors.items():
            if not self.is_valid_color(color):
                invalid_colors.append(f"{name} color '{color}'")
        
        if invalid_colors:
            error_msg = "Invalid colors detected:\n" + "\n".join(invalid_colors)
            self.toast.show_toast(error_msg, duration=5000, error=True)
            return False
        return True

    def on_color_change(self, event=None):
        """Validate colors when they change"""
        self.validate_colors()

    def generate_map(self):
        if self.map_canvas:
            self.map_canvas.get_tk_widget().destroy()
        
        # Validate colors first
        if not self.validate_colors():
            self.download_button.config(state="disabled")
            return
        
        gdf_copy = self.gdf.copy()
        # Initialize all counties as white
        gdf_copy["Color"] = "white"
        
        fam = self.selected_family.get().strip()
        gen = self.selected_genus.get().strip()
        spec = self.selected_species.get().strip()
        year = self.year_var.get().strip()
        
        if not fam or fam == "Select Family" or not gen or gen == "Select Genus" or not spec or spec == "Select Species":
            messagebox.showerror("Missing Input", "Please select Family, Genus, and Species.")
            self.download_button.config(state="disabled")
            return
        
        # Start with base DataFrame
        filtered = self.df
        
        # Apply family filter
        if fam == "All":
            filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
        elif fam == "Not Specified":
            filtered = filtered[filtered["family"].isna() | (filtered["family"].str.strip() == "")]
        else:
            filtered = filtered[filtered["family"].str.lower() == fam.lower()]
        
        # Apply genus filter
        if gen == "All":
            filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
        elif gen == "Not Specified":
            filtered = filtered[filtered["genus"].isna() | (filtered["genus"].str.strip() == "")]
        else:
            filtered = filtered[filtered["genus"].str.lower() == gen.lower()]
        
        # Apply species filter
        if spec == "all":
            filtered = filtered[filtered["species"].notna() & (filtered["species"].str.strip() != "")]
        elif spec == "not specified":
            filtered = filtered[filtered["species"].isna() | (filtered["species"].str.strip() == "")]
        else:
            filtered = filtered[filtered["species"].str.lower() == spec.lower()]
        
        # Create a set of valid county names from the shapefile for quick lookup
        valid_counties = set(self.standardize_county_names(gdf_copy["County"]))
        
        # Track unmatched counties to report to user
        unmatched_counties = set()
        
        if isinstance(year, str) and year.isdigit():
            year = int(year)
            # First pass: Mark counties with post-dividing year records with post_color (only if they're white)
            for county in filtered[filtered["year"] > year]["county"].unique():
                county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
                if county_lower in valid_counties:
                    mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                    if gdf_copy.loc[mask, "Color"].iloc[0] == "white":
                        gdf_copy.loc[mask, "Color"] = self.post_color.get()
                else:
                    unmatched_counties.add(county)
            
            # Second pass: Mark counties with pre-dividing year records with pre_color (overriding any color)
            for county in filtered[filtered["year"] <= year]["county"].unique():
                county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
                if county_lower in valid_counties:
                    mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                    gdf_copy.loc[mask, "Color"] = self.pre_color.get()
                else:
                    unmatched_counties.add(county)
        else:
            # If no year specified, mark all counties with records using all_color
            for county in filtered["county"].unique():
                county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
                if county_lower in valid_counties:
                    mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                    gdf_copy.loc[mask, "Color"] = self.all_color.get()
                else:
                    unmatched_counties.add(county)
        
        # Report any unmatched counties
        if unmatched_counties:
            print("\nWarning: The following counties in your Excel file don't match the shapefile counties:")
            print("-------------------------------------------------------------------------")
            for county in sorted(unmatched_counties):
                print(f"• {county}")
            print("\nValid Montana county names:")
            print("--------------------------------")
            for county in sorted(valid_counties):
                print(f"• {county}")
            print("-------------------------------------------------------------------------\n")
            
            messagebox.showwarning("County Name Mismatch",
                f"Some counties in your Excel file don't match the shapefile counties.\n\n"
                f"Number of unmatched counties: {len(unmatched_counties)}\n\n"
                f"Please check the console output for details and ensure county names match exactly."
            )
        
        # Create figure with calculated size
        fig = self.plt.figure(figsize=(12, 11))
        
        # Create main map axis with sufficient space for title and legend
        ax = fig.add_axes([0.1, 0.2, 0.8, 0.6])  # Adjusted to leave more space at top
        gdf_copy.boundary.plot(ax=ax, linewidth=1, edgecolor="black")
        gdf_copy.plot(ax=ax, color=gdf_copy["Color"], alpha=0.6)
        
        # Add title with dynamic font size
        title = f"{fam.title()} > {gen.title()} > {spec.lower()}"
        if isinstance(year, int):
            subtitle = f"\nYear: {year}"
            title += subtitle
        
        # Calculate dynamic font size based on figure width
        title_fontsize = min(15, max(8, fig.get_figwidth() * 1.5))
        ax.set_title(title, fontsize=title_fontsize, pad=25, wrap=True)  # Increased padding
        ax.axis("off")
        
        # Create a separate axis for the legend
        legend_ax = fig.add_axes([0.2, 0.02, 0.6, 0.12])
        legend_ax.axis('off')
        
        # Add a box around the legend
        legend_box = self.plt.Rectangle((0, 0), 1, 1, 
                                      facecolor='white', edgecolor='black',
                                      transform=legend_ax.transAxes)
        legend_ax.add_patch(legend_box)
        
        if isinstance(year, int):
            # Add horizontal bars with their descriptions for year analysis
            bar_length = 0.15
            bar_height = 0.1
            
            # Pre-year color
            legend_ax.add_patch(
                self.plt.Rectangle((0.2, 0.5), bar_length, bar_height, 
                                 facecolor=self.pre_color.get(), 
                                 alpha=0.6,
                                 edgecolor='black')
            )
            legend_ax.text(0.4, 0.55, f"Records ≤ {year}", fontsize=10, va='center')
            
            # Post-year color
            legend_ax.add_patch(
                self.plt.Rectangle((0.2, 0.2), bar_length, bar_height, 
                                 facecolor=self.post_color.get(), 
                                 alpha=0.6,
                                 edgecolor='black')
            )
            legend_ax.text(0.4, 0.25, f"Records > {year}", fontsize=10, va='center')
        else:
            # Add single bar for all records
            bar_length = 0.15
            bar_height = 0.1
            legend_ax.add_patch(
                self.plt.Rectangle((0.2, 0.35), bar_length, bar_height, 
                                 facecolor=self.all_color.get(), 
                                 alpha=0.6,
                                 edgecolor='black')
            )
            legend_ax.text(0.4, 0.4, "All Records", fontsize=10, va='center')
        
        # Set the legend axis limits
        legend_ax.set_xlim(0, 1)
        legend_ax.set_ylim(0, 1)
        
        # Adjust main plot to make room for legend
        fig.subplots_adjust(bottom=0.2, top=0.9)
        
        self.map_canvas = self.FigureCanvasTkAgg(fig, master=self.right_panel)
        self.map_canvas.draw()
        self.map_canvas.get_tk_widget().pack(fill='both', expand=True)
        
        self.current_fig = fig
        self.download_button.config(state="normal")
        print("✅ Map generated successfully!")
    
    def download_map(self):
        if self.current_fig:
            try:
                # Get Downloads folder path
                downloads_path = str(Path.home() / "Downloads")
                
                # Get taxonomic info and year
                fam = self.selected_family.get().strip().title()
                gen = self.selected_genus.get().strip().title()
                spec = self.selected_species.get().strip().lower()
                year = self.year_var.get().strip()
                
                # Create timestamp with current date (YYYYMMDD_HHMM)
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
                
                # Create filename with taxonomic info
                year_info = f"_{year}" if year else ""
                filename = f"{fam}-{gen}-{spec}{year_info}_{timestamp}.tiff"
                
                # Full path for the file
                file_path = os.path.join(downloads_path, filename)
                
                # Save the figure
                self.current_fig.savefig(file_path, format="tiff", dpi=300)
                
                # Show toast notification
                self.toast.show_toast(f"Map saved as {filename}")
                
                print(f"✅ TIFF map saved as '{file_path}'")
                
            except Exception as e:
                messagebox.showerror("Error", 
                    f"Error saving file:\n{str(e)}\n\n"
                    "Please try again."
                )
    
    def on_window_resize(self, event=None):
        try:
            # Get current window dimensions
            current_width = self.root.winfo_width()
            current_height = self.root.winfo_height()
            
            # Only process if we have valid dimensions
            if current_width <= 1 or current_height <= 1:
                return
                
            # Store current dimensions for restoring from maximize
            if not hasattr(self, 'last_valid_size'):
                self.last_valid_size = {
                    'width': current_width,
                    'height': current_height,
                    'x': self.root.winfo_x(),
                    'y': self.root.winfo_y()
                }
            
            # Update panels if needed
            self.update_panel_sizes()
                
        except Exception as e:
            print(f"Warning: Resize handling error: {str(e)}")  # For debugging
            pass  # Silently handle any errors during resize

    def update_panel_sizes(self):
        try:
            current_width = self.root.winfo_width()
            
            # Calculate new left panel width (20% of window width, max 300px)
            new_left_width = min(300, max(200, int(current_width * 0.2)))
            
            # Update left panel container if it exists
            if hasattr(self, 'left_panel_container'):
                self.left_panel_container.configure(width=new_left_width)
            
            # Update scroll canvas if it exists and has content
            if hasattr(self, 'scroll_canvas'):
                try:
                    # Only update if canvas exists and has content
                    if self.scroll_canvas.find_all():
                        self.scroll_canvas.configure(width=new_left_width-20)
                        self.scroll_canvas.itemconfig(
                            self.scroll_canvas.find_all()[0],
                            width=new_left_width-20
                        )
                except tk.TclError:
                    pass  # Handle case where canvas is being destroyed
            
            # Redraw map if it exists
            if hasattr(self, 'map_canvas') and self.map_canvas:
                try:
                    self.map_canvas.draw()
                except Exception:
                    pass  # Handle any drawing errors silently
                    
            # Force geometry update
            self.root.update_idletasks()
                
        except Exception as e:
            print(f"Warning: Panel update error: {str(e)}")  # For debugging
            pass  # Silently handle any errors

    def initialize_gui(self):
        # Configure style
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TLabel', font=('Helvetica', 10))
        style.configure('TButton', font=('Helvetica', 10))
        style.configure('Header.TLabel', font=('Helvetica', 12, 'bold'))
        style.configure('Title.TLabel', font=('Helvetica', 24, 'bold'))  # Increased from 16 to 24
        style.configure('Back.TButton', 
                       font=('Helvetica', 10, 'bold'),
                       padding=8)
        
        # Create main container with padding
        main_container = ttk.Frame(self.root, padding="20")
        main_container.pack(fill='both', expand=True)
        
        # Title
        title_frame = ttk.Frame(main_container)
        title_frame.pack(fill='x', pady=(0, 20))
        title_label = ttk.Label(
            title_frame, 
            text="Montana Bee County Map Generator - Single Year Analysis", 
            style='Title.TLabel',
            foreground='dark green'
        )
        title_label.pack()
        
        # Calculate left and right panel widths based on screen size
        screen_width = self.root.winfo_screenwidth()
        left_panel_width = min(300, int(screen_width * 0.2))  # 20% of screen width or 300px, whichever is smaller
        
        # Create left panel for controls with scrollbar
        self.left_panel_container = ttk.Frame(main_container, width=left_panel_width)  # Start with default width
        self.left_panel_container.pack(side='left', fill='y', padx=(0, 20))
        self.left_panel_container.pack_propagate(False)  # Prevent the frame from shrinking
        
        # Back button at the top of left panel container (outside scroll area)
        back_button = ttk.Button(
            self.left_panel_container,
            text="← Back to Selection",
            command=self.go_back,
            style='Back.TButton'
        )
        back_button.pack(fill='x', pady=(0, 10))
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(self.left_panel_container)
        scrollbar.pack(side='right', fill='y')
        
        # Create canvas for scrollable content with matching background
        bg_color = style.lookup('TFrame', 'background')  # Get ttk frame background color
        if not bg_color:  # Fallback if style lookup fails
            bg_color = self.left_panel_container.cget('background')
        self.scroll_canvas = Canvas(
            self.left_panel_container,
            yscrollcommand=scrollbar.set,
            width=left_panel_width-20,
            bg=bg_color,
            highlightthickness=0  # Remove border
        )
        self.scroll_canvas.pack(side='left', fill='y')
        
        scrollbar.config(command=self.scroll_canvas.yview)
        
        # Create frame for controls inside canvas
        left_panel = ttk.Frame(self.scroll_canvas)
        self.scroll_canvas.create_window((0, 0), window=left_panel, anchor='nw', width=left_panel_width-20)
        
        # Excel Load Section
        excel_frame = ttk.LabelFrame(left_panel, text="Data Input", padding="10")
        excel_frame.pack(fill='x', pady=(0, 20))
        
        # Create a frame for the button and file info
        file_info_frame = ttk.Frame(excel_frame)
        file_info_frame.pack(fill='x', expand=True)
        
        # Load button
        load_button = ttk.Button(
            file_info_frame, 
            text="Load Excel File", 
            command=self.load_excel, 
            style='TButton'
        )
        load_button.pack(fill='x', pady=(0, 5))
        
        # File info label
        file_label = ttk.Label(
            file_info_frame, 
            textvariable=self.selected_file_var,
            style='FileInfo.TLabel',
            wraplength=250
        )
        file_label.pack(fill='x')
        
        # Configure style for file info
        style.configure('FileInfo.TLabel', 
                       font=('Helvetica', 9),
                       foreground='dark green')
        
        # Color Selection Section
        color_frame = ttk.LabelFrame(left_panel, text="Color Settings", padding="10")
        color_frame.pack(fill='x', pady=(0, 20))
        
        # Common color options
        color_options = ["red", "blue", "grey", "black", "yellow", "purple", "orange", "pink", "brown"]
        
        # Pre Color
        ttk.Label(color_frame, text="Pre-Year Color:", style='TLabel').pack(fill='x')
        pre_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.pre_color, 
            values=color_options,
            state="normal"
        )
        pre_color_combo.pack(fill='x', pady=(0, 10))
        pre_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        pre_color_combo.bind('<Return>', self.on_color_change)
        pre_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Post Color
        ttk.Label(color_frame, text="Post-Year Color:", style='TLabel').pack(fill='x')
        post_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.post_color, 
            values=color_options,
            state="normal"
        )
        post_color_combo.pack(fill='x', pady=(0, 10))
        post_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        post_color_combo.bind('<Return>', self.on_color_change)
        post_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # All Color
        ttk.Label(color_frame, text="All Records Color:", style='TLabel').pack(fill='x')
        all_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.all_color, 
            values=color_options,
            state="normal"
        )
        all_color_combo.pack(fill='x', pady=(0, 10))
        all_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        all_color_combo.bind('<Return>', self.on_color_change)
        all_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Add helper text for custom colors
        helper_label = ttk.Label(
            color_frame, 
            text="Tip: You can type any valid color name or hex code (e.g., #FF5733)",
            style='TLabel',
            wraplength=250
        )
        helper_label.pack(fill='x', pady=(5, 0))
        
        # Year Input Section
        year_frame = ttk.LabelFrame(left_panel, text="Year Filter", padding="10")
        year_frame.pack(fill='x', pady=(0, 20))
        ttk.Label(year_frame, text="Year (Optional):", style='TLabel').pack(fill='x')
        ttk.Entry(year_frame, textvariable=self.year_var).pack(fill='x')
        
        # Species Selection Section
        species_frame = ttk.LabelFrame(left_panel, text="Species Selection", padding="10")
        species_frame.pack(fill='x', pady=(0, 20))
        
        # Family
        ttk.Label(species_frame, text="Family:", style='TLabel').pack(fill='x')
        self.family_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_family, state="readonly")
        self.family_dropdown.pack(fill='x', pady=(0, 10))
        
        # Genus
        ttk.Label(species_frame, text="Genus:", style='TLabel').pack(fill='x')
        self.genus_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_genus, state="readonly")
        self.genus_dropdown.pack(fill='x', pady=(0, 10))
        
        # Species
        ttk.Label(species_frame, text="Species:", style='TLabel').pack(fill='x')
        self.species_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_species, state="readonly")
        self.species_dropdown.pack(fill='x', pady=(0, 10))
        
        # Button Section
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill='x', pady=(0, 20))
        
        # Style configuration for buttons
        style.configure('Action.TButton', 
                       font=('Helvetica', 10, 'bold'),
                       padding=10)
        
        generate_button = ttk.Button(
            button_frame, 
            text="Generate Map", 
            command=self.generate_map, 
            style='Action.TButton'
        )
        generate_button.pack(fill='x', pady=(0, 5))
        
        self.download_button = ttk.Button(
            button_frame, 
            text="Download Map", 
            command=self.download_map, 
            state="disabled", 
            style='Action.TButton'
        )
        self.download_button.pack(fill='x', pady=(0, 5))
        
        # Create right panel for map display with dynamic width
        self.right_panel = ttk.Frame(main_container)
        self.right_panel.pack(side='left', fill='both', expand=True)
        
        # Bind resize event to the main update function
        self.root.bind('<Configure>', self.on_window_resize)
        
        # Store initial window size
        self.last_valid_size = {
            'width': self.root.winfo_width(),
            'height': self.root.winfo_height(),
            'x': self.root.winfo_x(),
            'y': self.root.winfo_y()
        }
        
        # Remove the local update_panel_sizes function and use the class method instead
        
        # Configure scrolling
        def on_configure(event):
            self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox('all'))
        
        left_panel.bind('<Configure>', on_configure)
        
        # Enable mousewheel scrolling
        def on_mousewheel(event):
            self.scroll_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        
        self.scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        # Bind dropdowns
        self.family_dropdown.bind("<<ComboboxSelected>>", self.update_genus_dropdown)
        self.genus_dropdown.bind("<<ComboboxSelected>>", self.update_species_dropdown)
    
    def go_back(self):
        """Return to the selection screen"""
        # Store current position and state
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        current_state = self.root.state()
        self.root.destroy()
        selection = SelectionScreen(self.root.master, self.main_app, from_analysis=True)
        # Set position first
        selection.root.geometry(f"{selection.root.winfo_screenwidth()}x{selection.root.winfo_screenheight()}+{x}+{y}")
        selection.root.update()  # Force geometry update
        # Restore the previous window state
        if current_state == 'zoomed':
            selection.root.state('zoomed')

class DualYearAnalysis:
    def __init__(self, parent, main_app):
        self.root = tk.Toplevel(parent)
        self.root.title("Dual Year Analysis - Montana Bee County Map Generator")
        
        # Allow the window to be moved to any screen
        self.root.attributes('-alpha', 1.0)
        self.root.attributes('-topmost', False)
        
        # Store main_app reference
        self.main_app = main_app
        
        # Get dependencies from main_app
        self.pd = main_app.pd
        self.plt = main_app.plt
        self.FigureCanvasTkAgg = main_app.FigureCanvasTkAgg
        
        # Set window icon
        if getattr(sys, 'frozen', False):
            base_dir = sys._MEIPASS
        else:
            base_dir = os.path.dirname(os.path.abspath(__file__))
        icon_path = os.path.join(base_dir, "app_icon.ico")
        if os.path.exists(icon_path):
            self.root.iconbitmap(icon_path)
        
        # Get screen dimensions
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        
        # Set initial size but don't maximize yet - that will be done by the calling method
        self.root.geometry(f"{screen_width}x{screen_height}")
        
        # Ensure window can be moved and resized
        self.root.resizable(True, True)
        self.root.minsize(800, 600)
        
        # Initialize variables
        self.map_canvas = None
        self.current_fig = None
        
        # Initialize StringVar variables
        self.first_color = StringVar(self.root)
        self.second_color = StringVar(self.root)
        self.third_color = StringVar(self.root)
        self.first_year_var = StringVar(self.root)
        self.second_year_var = StringVar(self.root)
        self.selected_family = StringVar(self.root)
        self.selected_genus = StringVar(self.root)
        self.selected_species = StringVar(self.root)
        self.selected_file_var = StringVar(self.root)
        
        # Set default values
        self.first_color.set("grey")  # For records ≤ first year
        self.second_color.set("red")  # For records between years
        self.third_color.set("yellow")  # For records > second year
        self.selected_file_var.set("No file selected")
        
        # Create toast notification instance
        self.toast = ToastNotification(self.root)
        
        # Initialize pandas DataFrame
        self.df = self.pd.DataFrame()
        
        # Get the shapefile data from parent
        self.gdf = main_app.gdf.copy()
        
        # Initialize GUI
        self.initialize_gui()
        
        # Bind window state change
        self.root.bind("<Configure>", self.on_window_resize)

    def standardize_county_names(self, county_series):
        """
        Standardize county names by:
        1. Converting '&' to 'and'
        2. Stripping whitespace
        3. Converting to lowercase
        """
        return county_series.str.strip().str.lower().str.replace('&', 'and')

    def load_excel(self):
        path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
        if not path:
            return
        
        # Show loading dialog
        loading_window = tk.Tk()
        loading_window.title("Loading")
        screen_width = loading_window.winfo_screenwidth()
        screen_height = loading_window.winfo_screenheight()
        window_width = 300
        window_height = 100
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        loading_window.geometry(f"{window_width}x{window_height}+{x}+{y}")
        loading_window.overrideredirect(True)
        
        frame = ttk.Frame(loading_window, padding="20", relief="raised")
        frame.pack(fill='both', expand=True)
        
        loading_label = ttk.Label(frame, text="Loading file...\nPlease wait", font=('Helvetica', 10))
        loading_label.pack(pady=10)
        
        progress = ttk.Progressbar(frame, mode='indeterminate')
        progress.pack(fill='x', pady=5)
        progress.start(10)
        
        loading_window.update()
        
        try:
            # Load the Excel file
            self.df = self.pd.read_excel(path, sheet_name=0)
            
            # # Print all column names immediately after loading
            # print("\nAvailable columns in Excel file:")
            # print("--------------------------------")
            # for i, col in enumerate(self.df.columns, 1):
            #     print(f"{i}. {col}")
            # print("--------------------------------\n")
            
            # Now continue with column processing
            self.df.columns = self.df.columns.str.strip()
            
            # Get just the filename from the path
            filename = path.split('/')[-1]
            
            # Validate required columns
            required_columns = ["county", "family", "genus", "species", "year"]
            missing_columns = [col for col in required_columns if col not in self.df.columns]
            
            if missing_columns:
                progress.stop()
                loading_window.destroy()
                self.selected_file_var.set("No file selected")
                messagebox.showerror("Error", 
                    f"Missing required columns: {', '.join(missing_columns)}\n\n"
                    "The following columns are required:\n"
                    "- county: for mapping locations\n"
                    "- family: for taxonomic classification\n"
                    "- genus: for taxonomic classification\n"
                    "- species: for taxonomic classification\n"
                    "- year: for temporal analysis\n\n"
                    "Please check your Excel file and try again."
                )
                return
            
            # Process the data
            # First standardize county names
            self.df["county"] = self.standardize_county_names(self.df["county"])
            
            # Process other columns
            for col in ["family", "genus", "species"]:
                self.df[col] = self.df[col].astype(str).str.strip().str.lower()
            self.df["year"] = self.df["year"].astype(str).str.extract(r'(\d{4})').astype(float)
            
            # Get valid Montana counties from shapefile
            valid_counties = set(self.standardize_county_names(self.gdf["County"]))
            
            # Filter DataFrame to only include valid Montana counties
            montana_records = self.df[self.df["county"].isin(valid_counties)]
            
            if len(montana_records) == 0:
                progress.stop()
                loading_window.destroy()
                self.selected_file_var.set("No file selected")
                messagebox.showerror("Error", 
                    "No valid Montana county records found in the Excel file.\n\n"
                    "Please check that your data contains Montana county records."
                )
                return
            
            # Replace the main DataFrame with only Montana records
            self.df = montana_records
            
            # Calculate statistics using Montana records
            num_records = len(montana_records)
            num_families = len(montana_records["family"].unique())
            num_genera = len(montana_records["genus"].unique())
            num_species = len(montana_records["species"].unique())
            year_range = f"{int(montana_records['year'].min())} - {int(montana_records['year'].max())}"
            num_counties = len(montana_records["county"].unique())
            
            # Get valid families (non-empty/non-null values)
            valid_families = sorted(montana_records["family"].dropna().unique())
            valid_families = [f for f in valid_families if str(f).strip() and str(f).lower() != 'nan']  # Remove empty strings and 'nan'
            
            # Capitalize family names
            family_values = ["All"] + [f.title() for f in valid_families]
            
            # Update Family dropdown
            self.family_dropdown["values"] = family_values
            self.family_dropdown.set("All")
            
            # Trigger genus dropdown update
            self.update_genus_dropdown()
            
            # Set genus to "All" and trigger species update
            self.genus_dropdown.set("All")
            self.update_species_dropdown()
            
            # Set species to "all"
            self.species_dropdown.set("all")
            
            # Reset year fields
            self.first_year_var.set("")
            self.second_year_var.set("")
            
            # Update file info display
            self.selected_file_var.set(f"✓ {filename}\n{num_records:,} Montana records loaded")
            
            # Stop progress bar and close loading window
            progress.stop()
            loading_window.destroy()
            
            # Show success message with detailed statistics
            messagebox.showinfo("Success", 
                f"File loaded successfully!\n\n"
                f"Montana Dataset Summary:\n"
                f"• Total Records: {num_records:,}\n"
                f"• Unique Families: {num_families}\n"
                f"• Unique Genera: {num_genera}\n"
                f"• Unique Species: {num_species}\n"
                f"• Counties Covered: {num_counties}\n"
                f"• Year Range: {year_range}\n\n"
                "Please select a Family to continue."
            )
            
            print("✅ Excel file loaded successfully!")
            
        except Exception as e:
            progress.stop()
            loading_window.destroy()
            self.selected_file_var.set("No file selected")
            error_message = str(e)
            if "No sheet named" in error_message:
                error_message = "Invalid Excel file format. Please ensure your data is in the first sheet."
            elif "Invalid file" in error_message:
                error_message = "Invalid file format. Please ensure you're uploading a valid Excel (.xlsx) file."
            
            messagebox.showerror("Error", 
                f"Error loading file:\n{error_message}\n\n"
                "Please check your Excel file format and try again."
            )
    
    def update_genus_dropdown(self, event=None):
        family = self.selected_family.get().strip()
        
        if family == "Select Family":
            self.genus_dropdown["values"] = []
            self.genus_dropdown.set("Select Genus")
            return
        
        # Filter based on family selection
        if family == "All":
            # Get all non-empty genus values
            filtered = self.df[self.df["genus"].notna() & (self.df["genus"].str.strip() != "")]
        else:
            # Get genus for specific family (case-insensitive)
            filtered = self.df[self.df["family"].str.lower() == family.lower()]
        
        # Get valid genera (non-empty/non-null values)
        valid_genera = sorted(filtered["genus"].dropna().unique())
        valid_genera = [g for g in valid_genera if str(g).strip() and str(g).lower() != 'nan']  # Remove empty strings and 'nan'
        
        # Create genus list with special options
        genus_values = ["All"] + [g.title() for g in valid_genera]
        
        # Update Genus dropdown
        self.genus_dropdown["values"] = genus_values
        self.genus_dropdown.set("Select Genus")
        
        # Reset species dropdown
        self.species_dropdown.set("Select Species")
        self.species_dropdown["values"] = []
    
    def update_species_dropdown(self, event=None):
        family = self.selected_family.get().strip()
        genus = self.selected_genus.get().strip()
        
        if family == "Select Family" or genus == "Select Genus":
            self.species_dropdown["values"] = []
            self.species_dropdown.set("Select Species")
            return
        
        # Start with base DataFrame
        filtered = self.df
        
        # Apply family filter
        if family == "All":
            filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
        else:
            filtered = filtered[filtered["family"].str.lower() == family.lower()]
        
        # Apply genus filter
        if genus == "All":
            filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
        else:
            filtered = filtered[filtered["genus"].str.lower() == genus.lower()]
        
        # Get valid species (non-empty/non-null values)
        valid_species = sorted(filtered["species"].dropna().unique())
        valid_species = [s for s in valid_species if str(s).strip() and str(s).lower() != 'nan']  # Remove empty strings and 'nan'
        
        # Create species list with special options - note lowercase for species
        species_values = ["all"] + valid_species
        
        # Update Species dropdown
        self.species_dropdown["values"] = species_values
        self.species_dropdown.set("Select Species")
    
    def is_valid_color(self, color):
        """Validate if a color string is a valid matplotlib color"""
        try:
            # Convert color to RGB
            self.plt.matplotlib.colors.to_rgb(color)
            return True
        except ValueError:
            return False

    def validate_colors(self):
        """Validate all selected colors"""
        colors = {
            'First Period': self.first_color.get(),
            'Second Period': self.second_color.get(),
            'Third Period': self.third_color.get()
        }
        
        invalid_colors = []
        for name, color in colors.items():
            if not self.is_valid_color(color):
                invalid_colors.append(f"{name} color '{color}'")
        
        if invalid_colors:
            error_msg = "Invalid colors detected:\n" + "\n".join(invalid_colors)
            self.toast.show_toast(error_msg, duration=5000, error=True)
            return False
        return True

    def validate_years(self):
        """Validate that both years are provided and first year is less than second year"""
        try:
            first_year = int(self.first_year_var.get().strip())
            second_year = int(self.second_year_var.get().strip())
            
            if first_year >= second_year:
                messagebox.showerror("Invalid Years", 
                    "First year must be less than second year.\n\n"
                    f"You entered:\n"
                    f"First Year: {first_year}\n"
                    f"Second Year: {second_year}"
                )
                return False
            return True
        except ValueError:
            messagebox.showerror("Invalid Years",
                "Both years are required and must be valid numbers.\n\n"
                "Please enter valid years (e.g., 1990, 2000)"
            )
            return False

    def on_color_change(self, event=None):
        """Validate colors when they change"""
        self.validate_colors()

    def generate_map(self):
        if self.map_canvas:
            self.map_canvas.get_tk_widget().destroy()
        
        # Validate colors first
        if not self.validate_colors():
            self.download_button.config(state="disabled")
            return
        
        # Validate years
        if not self.validate_years():
            self.download_button.config(state="disabled")
            return
        
        gdf_copy = self.gdf.copy()
        # Initialize all counties as white
        gdf_copy["Color"] = "white"
        
        fam = self.selected_family.get().strip()
        gen = self.selected_genus.get().strip()
        spec = self.selected_species.get().strip()
        first_year = int(self.first_year_var.get().strip())
        second_year = int(self.second_year_var.get().strip())
        
        if not fam or fam == "Select Family" or not gen or gen == "Select Genus" or not spec or spec == "Select Species":
            messagebox.showerror("Missing Input", "Please select Family, Genus, and Species.")
            self.download_button.config(state="disabled")
            return
        
        # Start with base DataFrame
        filtered = self.df
        
        # Apply family filter
        if fam == "All":
            filtered = filtered[filtered["family"].notna() & (filtered["family"].str.strip() != "")]
        elif fam == "Not Specified":
            filtered = filtered[filtered["family"].isna() | (filtered["family"].str.strip() == "")]
        else:
            filtered = filtered[filtered["family"].str.lower() == fam.lower()]
        
        # Apply genus filter
        if gen == "All":
            filtered = filtered[filtered["genus"].notna() & (filtered["genus"].str.strip() != "")]
        elif gen == "Not Specified":
            filtered = filtered[filtered["genus"].isna() | (filtered["genus"].str.strip() == "")]
        else:
            filtered = filtered[filtered["genus"].str.lower() == gen.lower()]
        
        # Apply species filter
        if spec == "all":
            filtered = filtered[filtered["species"].notna() & (filtered["species"].str.strip() != "")]
        elif spec == "not specified":
            filtered = filtered[filtered["species"].isna() | (filtered["species"].str.strip() == "")]
        else:
            filtered = filtered[filtered["species"].str.lower() == spec.lower()]
        
        # Create a set of valid county names from the shapefile for quick lookup
        valid_counties = set(self.standardize_county_names(gdf_copy["County"]))
        
        # Track unmatched counties to report to user
        unmatched_counties = set()
        
        # First pass: Mark counties with records > second_year with third_color (lowest priority)
        for county in filtered[filtered["year"] > second_year]["county"].unique():
            county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
            if county_lower in valid_counties:
                mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                if gdf_copy.loc[mask, "Color"].iloc[0] == "white":
                    gdf_copy.loc[mask, "Color"] = self.third_color.get()
            else:
                unmatched_counties.add(county)
        
        # Second pass: Mark counties with records between years with second_color (medium priority)
        for county in filtered[(filtered["year"] > first_year) & (filtered["year"] <= second_year)]["county"].unique():
            county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
            if county_lower in valid_counties:
                mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                gdf_copy.loc[mask, "Color"] = self.second_color.get()
            else:
                unmatched_counties.add(county)
        
        # Third pass: Mark counties with records ≤ first_year with first_color (highest priority)
        for county in filtered[filtered["year"] <= first_year]["county"].unique():
            county_lower = self.standardize_county_names(pd.Series([county])).iloc[0]
            if county_lower in valid_counties:
                mask = self.standardize_county_names(gdf_copy["County"]) == county_lower
                gdf_copy.loc[mask, "Color"] = self.first_color.get()
            else:
                unmatched_counties.add(county)
        
        # Report any unmatched counties
        if unmatched_counties:
            print("\nWarning: The following counties in your Excel file don't match the shapefile counties:")
            print("-------------------------------------------------------------------------")
            for county in sorted(unmatched_counties):
                print(f"• {county}")
            print("\nValid Montana county names:")
            print("--------------------------------")
            for county in sorted(valid_counties):
                print(f"• {county}")
            print("-------------------------------------------------------------------------\n")
            
            messagebox.showwarning("County Name Mismatch",
                f"Some counties in your Excel file don't match the shapefile counties.\n\n"
                f"Number of unmatched counties: {len(unmatched_counties)}\n\n"
                f"Please check the console output for details and ensure county names match exactly."
            )
        
        # Create figure with calculated size
        fig = self.plt.figure(figsize=(12, 11))
        
        # Create main map axis with sufficient space for title and legend
        ax = fig.add_axes([0.1, 0.2, 0.8, 0.6])  # Adjusted to leave more space at top
        gdf_copy.boundary.plot(ax=ax, linewidth=1, edgecolor="black")
        gdf_copy.plot(ax=ax, color=gdf_copy["Color"], alpha=0.6)
        
        # Add title with dynamic font size
        title = f"{fam.title()} > {gen.title()} > {spec.lower()}\nYears: {first_year} - {second_year}"
        ax.set_title(title, fontsize=15, pad=20)
        ax.axis("off")
        
        # Create a separate axis for the legend
        legend_ax = fig.add_axes([0.2, 0.02, 0.6, 0.12])
        legend_ax.axis('off')
        
        # Add a box around the legend
        legend_box = self.plt.Rectangle((0, 0), 1, 1, 
                                      facecolor='white', edgecolor='black',
                                      transform=legend_ax.transAxes)
        legend_ax.add_patch(legend_box)
        
        # Add horizontal bars with their descriptions
        bar_length = 0.15
        bar_height = 0.1
        
        # First period color (highest priority)
        legend_ax.add_patch(
            self.plt.Rectangle((0.2, 0.7), bar_length, bar_height, 
                             facecolor=self.first_color.get(), 
                             alpha=0.6,
                             edgecolor='black')
        )
        legend_ax.text(0.4, 0.75, f"Records ≤ {first_year}", fontsize=10, va='center')
        
        # Second period color (medium priority)
        legend_ax.add_patch(
            self.plt.Rectangle((0.2, 0.4), bar_length, bar_height, 
                             facecolor=self.second_color.get(), 
                             alpha=0.6,
                             edgecolor='black')
        )
        legend_ax.text(0.4, 0.45, f"Records {first_year+1} - {second_year}", fontsize=10, va='center')
        
        # Third period color (lowest priority)
        legend_ax.add_patch(
            self.plt.Rectangle((0.2, 0.1), bar_length, bar_height, 
                             facecolor=self.third_color.get(), 
                             alpha=0.6,
                             edgecolor='black')
        )
        legend_ax.text(0.4, 0.15, f"Records > {second_year}", fontsize=10, va='center')
        
        # Set the legend axis limits
        legend_ax.set_xlim(0, 1)
        legend_ax.set_ylim(0, 1)
        
        # Adjust main plot to make room for legend
        fig.subplots_adjust(bottom=0.2, top=0.9)
        
        self.map_canvas = self.FigureCanvasTkAgg(fig, master=self.right_panel)
        self.map_canvas.draw()
        self.map_canvas.get_tk_widget().pack(fill='both', expand=True)
        
        self.current_fig = fig
        self.download_button.config(state="normal")
        print("✅ Map generated successfully!")
    
    def download_map(self):
        if self.current_fig:
            try:
                # Get Downloads folder path
                downloads_path = str(Path.home() / "Downloads")
                
                # Get taxonomic info and years
                fam = self.selected_family.get().strip().title()
                gen = self.selected_genus.get().strip().title()
                spec = self.selected_species.get().strip().lower()
                first_year = self.first_year_var.get().strip()
                second_year = self.second_year_var.get().strip()
                
                # Create timestamp with current date (YYYYMMDD_HHMM)
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
                
                # Create filename with taxonomic info
                filename = f"{fam}-{gen}-{spec}_{first_year}-{second_year}_{timestamp}.tiff"
                
                # Full path for the file
                file_path = os.path.join(downloads_path, filename)
                
                # Save the figure
                self.current_fig.savefig(file_path, format="tiff", dpi=300)
                
                # Show toast notification
                self.toast.show_toast(f"Map saved as {filename}")
                
                print(f"✅ TIFF map saved as '{file_path}'")
                
            except Exception as e:
                messagebox.showerror("Error", 
                    f"Error saving file:\n{str(e)}\n\n"
                    "Please try again."
                )
    
    def on_window_resize(self, event=None):
        try:
            # Get current window dimensions
            current_width = self.root.winfo_width()
            current_height = self.root.winfo_height()
            
            # Only process if we have valid dimensions
            if current_width <= 1 or current_height <= 1:
                return
                
            # Store current dimensions for restoring from maximize
            if not hasattr(self, 'last_valid_size'):
                self.last_valid_size = {
                    'width': current_width,
                    'height': current_height,
                    'x': self.root.winfo_x(),
                    'y': self.root.winfo_y()
                }
            
            # Update panels if needed
            self.update_panel_sizes()
                
        except Exception as e:
            print(f"Warning: Resize handling error: {str(e)}")  # For debugging
            pass  # Silently handle any errors during resize

    def update_panel_sizes(self):
        try:
            current_width = self.root.winfo_width()
            
            # Calculate new left panel width (20% of window width, max 300px)
            new_left_width = min(300, max(200, int(current_width * 0.2)))
            
            # Update left panel container if it exists
            if hasattr(self, 'left_panel_container'):
                self.left_panel_container.configure(width=new_left_width)
            
            # Update scroll canvas if it exists and has content
            if hasattr(self, 'scroll_canvas'):
                try:
                    # Only update if canvas exists and has content
                    if self.scroll_canvas.find_all():
                        self.scroll_canvas.configure(width=new_left_width-20)
                        self.scroll_canvas.itemconfig(
                            self.scroll_canvas.find_all()[0],
                            width=new_left_width-20
                        )
                except tk.TclError:
                    pass  # Handle case where canvas is being destroyed
            
            # Redraw map if it exists
            if hasattr(self, 'map_canvas') and self.map_canvas:
                try:
                    self.map_canvas.draw()
                except Exception:
                    pass  # Handle any drawing errors silently
                    
            # Force geometry update
            self.root.update_idletasks()
                
        except Exception as e:
            print(f"Warning: Panel update error: {str(e)}")  # For debugging
            pass  # Silently handle any errors

    def initialize_gui(self):
        # Configure style
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TLabel', font=('Helvetica', 10))
        style.configure('TButton', font=('Helvetica', 10))
        style.configure('Header.TLabel', font=('Helvetica', 12, 'bold'))
        style.configure('Title.TLabel', font=('Helvetica', 24, 'bold'))  # Increased from 16 to 24
        style.configure('Back.TButton', 
                       font=('Helvetica', 10, 'bold'),
                       padding=8)
        
        # Create main container with padding
        main_container = ttk.Frame(self.root, padding="20")
        main_container.pack(fill='both', expand=True)
        
        # Title
        title_frame = ttk.Frame(main_container)
        title_frame.pack(fill='x', pady=(0, 20))
        title_label = ttk.Label(
            title_frame, 
            text="Montana Bee County Map Generator - Dual Year Analysis", 
            style='Title.TLabel',
            foreground='dark green'
        )
        title_label.pack()
        
        # Calculate left and right panel widths based on screen size
        screen_width = self.root.winfo_screenwidth()
        left_panel_width = min(300, int(screen_width * 0.2))  # 20% of screen width or 300px, whichever is smaller
        
        # Create left panel for controls with scrollbar
        self.left_panel_container = ttk.Frame(main_container, width=left_panel_width)  # Start with default width
        self.left_panel_container.pack(side='left', fill='y', padx=(0, 20))
        self.left_panel_container.pack_propagate(False)  # Prevent the frame from shrinking
        
        # Back button at the top of left panel container (outside scroll area)
        back_button = ttk.Button(
            self.left_panel_container,
            text="← Back to Selection",
            command=self.go_back,
            style='Back.TButton'
        )
        back_button.pack(fill='x', pady=(0, 10))
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(self.left_panel_container)
        scrollbar.pack(side='right', fill='y')
        
        # Create canvas for scrollable content with matching background
        bg_color = style.lookup('TFrame', 'background')  # Get ttk frame background color
        if not bg_color:  # Fallback if style lookup fails
            bg_color = self.left_panel_container.cget('background')
        self.scroll_canvas = Canvas(
            self.left_panel_container,
            yscrollcommand=scrollbar.set,
            width=left_panel_width-20,
            bg=bg_color,
            highlightthickness=0  # Remove border
        )
        self.scroll_canvas.pack(side='left', fill='y')
        
        scrollbar.config(command=self.scroll_canvas.yview)
        
        # Create frame for controls inside canvas
        left_panel = ttk.Frame(self.scroll_canvas)
        self.scroll_canvas.create_window((0, 0), window=left_panel, anchor='nw', width=left_panel_width-20)
        
        # Excel Load Section
        excel_frame = ttk.LabelFrame(left_panel, text="Data Input", padding="10")
        excel_frame.pack(fill='x', pady=(0, 20))
        
        # Create a frame for the button and file info
        file_info_frame = ttk.Frame(excel_frame)
        file_info_frame.pack(fill='x', expand=True)
        
        # Load button
        load_button = ttk.Button(
            file_info_frame, 
            text="Load Excel File", 
            command=self.load_excel, 
            style='TButton'
        )
        load_button.pack(fill='x', pady=(0, 5))
        
        # File info label
        file_label = ttk.Label(
            file_info_frame, 
            textvariable=self.selected_file_var,
            style='FileInfo.TLabel',
            wraplength=250
        )
        file_label.pack(fill='x')
        
        # Configure style for file info
        style.configure('FileInfo.TLabel', 
                       font=('Helvetica', 9),
                       foreground='dark green')
        
        # Color Selection Section
        color_frame = ttk.LabelFrame(left_panel, text="Color Settings", padding="10")
        color_frame.pack(fill='x', pady=(0, 20))
        
        # Common color options
        color_options = ["red", "blue", "grey", "black", "yellow", "purple", "orange", "pink", "brown"]
        
        # First Color (Highest Priority)
        ttk.Label(color_frame, text="First Period Color (Highest Priority):", style='TLabel').pack(fill='x')
        first_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.first_color, 
            values=color_options,
            state="normal"
        )
        first_color_combo.pack(fill='x', pady=(0, 10))
        first_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        first_color_combo.bind('<Return>', self.on_color_change)
        first_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Second Color (Medium Priority)
        ttk.Label(color_frame, text="Second Period Color (Medium Priority):", style='TLabel').pack(fill='x')
        second_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.second_color, 
            values=color_options,
            state="normal"
        )
        second_color_combo.pack(fill='x', pady=(0, 10))
        second_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        second_color_combo.bind('<Return>', self.on_color_change)
        second_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Third Color (Lowest Priority)
        ttk.Label(color_frame, text="Third Period Color (Lowest Priority):", style='TLabel').pack(fill='x')
        third_color_combo = ttk.Combobox(
            color_frame, 
            textvariable=self.third_color, 
            values=color_options,
            state="normal"
        )
        third_color_combo.pack(fill='x', pady=(0, 10))
        third_color_combo.bind('<<ComboboxSelected>>', self.on_color_change)
        third_color_combo.bind('<Return>', self.on_color_change)
        third_color_combo.bind('<FocusOut>', self.on_color_change)
        
        # Add helper text for custom colors
        helper_label = ttk.Label(
            color_frame, 
            text="Tip: You can type any valid color name or hex code (e.g., #FF5733)",
            style='TLabel',
            wraplength=250
        )
        helper_label.pack(fill='x', pady=(5, 0))
        
        # Year Input Section
        year_frame = ttk.LabelFrame(left_panel, text="Year Filters", padding="10")
        year_frame.pack(fill='x', pady=(0, 20))
        
        # First Year (Required)
        ttk.Label(year_frame, text="First Year (Required):", style='TLabel').pack(fill='x')
        ttk.Entry(year_frame, textvariable=self.first_year_var).pack(fill='x', pady=(0, 10))
        
        # Second Year (Required)
        ttk.Label(year_frame, text="Second Year (Required):", style='TLabel').pack(fill='x')
        ttk.Entry(year_frame, textvariable=self.second_year_var).pack(fill='x')
        
        # Add helper text for years
        year_helper = ttk.Label(
            year_frame,
            text="Note: Records will be colored based on:\n" +
                 "• First color: Records ≤ First Year\n" +
                 "• Second color: Records between years\n" +
                 "• Third color: Records > Second Year",
            style='TLabel',
            wraplength=250
        )
        year_helper.pack(fill='x', pady=(10, 0))
        
        # Species Selection Section
        species_frame = ttk.LabelFrame(left_panel, text="Species Selection", padding="10")
        species_frame.pack(fill='x', pady=(0, 20))
        
        # Family
        ttk.Label(species_frame, text="Family:", style='TLabel').pack(fill='x')
        self.family_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_family, state="readonly")
        self.family_dropdown.pack(fill='x', pady=(0, 10))
        
        # Genus
        ttk.Label(species_frame, text="Genus:", style='TLabel').pack(fill='x')
        self.genus_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_genus, state="readonly")
        self.genus_dropdown.pack(fill='x', pady=(0, 10))
        
        # Species
        ttk.Label(species_frame, text="Species:", style='TLabel').pack(fill='x')
        self.species_dropdown = ttk.Combobox(species_frame, textvariable=self.selected_species, state="readonly")
        self.species_dropdown.pack(fill='x', pady=(0, 10))
        
        # Button Section
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill='x', pady=(0, 20))
        
        # Style configuration for buttons
        style.configure('Action.TButton', 
                       font=('Helvetica', 10, 'bold'),
                       padding=10)
        
        generate_button = ttk.Button(
            button_frame, 
            text="Generate Map", 
            command=self.generate_map, 
            style='Action.TButton'
        )
        generate_button.pack(fill='x', pady=(0, 5))
        
        self.download_button = ttk.Button(
            button_frame, 
            text="Download Map", 
            command=self.download_map, 
            state="disabled", 
            style='Action.TButton'
        )
        self.download_button.pack(fill='x', pady=(0, 5))
        
        # Create right panel for map display with dynamic width
        self.right_panel = ttk.Frame(main_container)
        self.right_panel.pack(side='left', fill='both', expand=True)
        
        # Bind resize event to the main update function
        self.root.bind('<Configure>', self.on_window_resize)
        
        # Store initial window size
        self.last_valid_size = {
            'width': self.root.winfo_width(),
            'height': self.root.winfo_height(),
            'x': self.root.winfo_x(),
            'y': self.root.winfo_y()
        }
        
        # Remove the local update_panel_sizes function and use the class method instead
        
        # Configure scrolling
        def on_configure(event):
            self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox('all'))
        
        left_panel.bind('<Configure>', on_configure)
        
        # Enable mousewheel scrolling
        def on_mousewheel(event):
            self.scroll_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        
        self.scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        # Bind dropdowns
        self.family_dropdown.bind("<<ComboboxSelected>>", self.update_genus_dropdown)
        self.genus_dropdown.bind("<<ComboboxSelected>>", self.update_species_dropdown)
    
    def go_back(self):
        """Return to the selection screen"""
        # Store current position and state
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        current_state = self.root.state()
        self.root.destroy()
        selection = SelectionScreen(self.root.master, self.main_app, from_analysis=True)
        # Set position first
        selection.root.geometry(f"{selection.root.winfo_screenwidth()}x{selection.root.winfo_screenheight()}+{x}+{y}")
        selection.root.update()  # Force geometry update
        # Restore the previous window state
        if current_state == 'zoomed':
            selection.root.state('zoomed')

if __name__ == "__main__":
    app = MainApplication()