import xarray as xr
import rioxarray
import geopandas as gpd
import rasterio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import netCDF4
import os
import json
from matplotlib import colors as mcolors
from matplotlib import patches as mpatches

# Load and process chlorophyll data
chl_dataset = xr.open_dataset("C:\\Users\\claud\\OneDrive\\Documents\\SALENTO\\Semester 2\\Modeling Ecology\\Project\\PugliaFishingEffort\\Raw\\CopernicusData\\ChlA\\cmems_mod_med_bgc-pft_anfc_4.2km_P1M-m_1761296079827.nc")

# Extract total chlorophyll for the surface layer
chl_surface = chl_dataset['chl'].isel(depth=0)
print("\nChlorophyll data information:")
print("Time range:", chl_surface.time.values[0], "to", chl_surface.time.values[-1])
print("Latitude range:", chl_surface.latitude.values.min(), "to", chl_surface.latitude.values.max())
print("Longitude range:", chl_surface.longitude.values.min(), "to", chl_surface.longitude.values.max())
print("Shape of surface chlorophyll data:", chl_surface.shape)

# Create a multi-panel plot showing chlorophyll components
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
fig.suptitle('Chlorophyll Components - Surface Layer', fontsize=14)

# Plot total chlorophyll
chl_surface.isel(time=0).plot(
    ax=axes[0,0],
    cmap='viridis',
    robust=True,
    cbar_kwargs={'label': 'Total Chl-a (mg/m³)'}
)
axes[0,0].set_title('Total Chlorophyll-a')

# Plot diatoms
chl_dataset['diatoChla'].isel(time=0, depth=0).plot(
    ax=axes[0,1],
    cmap='viridis',
    robust=True,
    cbar_kwargs={'label': 'Diatom Chl-a (mg/m³)'}
)
axes[0,1].set_title('Diatom Chlorophyll-a')

# Plot nano-phytoplankton
chl_dataset['nanoChla'].isel(time=0, depth=0).plot(
    ax=axes[1,0],
    cmap='viridis',
    robust=True,
    cbar_kwargs={'label': 'Nano Chl-a (mg/m³)'}
)
axes[1,0].set_title('Nano-phytoplankton Chlorophyll-a')

# Plot pico-phytoplankton
chl_dataset['picoChla'].isel(time=0, depth=0).plot(
    ax=axes[1,1],
    cmap='viridis',
    robust=True,
    cbar_kwargs={'label': 'Pico Chl-a (mg/m³)'}
)
axes[1,1].set_title('Pico-phytoplankton Chlorophyll-a')

plt.tight_layout()
plt.savefig('chlorophyll_components.png', dpi=300, bbox_inches='tight')
print("\nPlot saved as chlorophyll_components.png")
# Load habitat data (user-provided path)
habitat_path = r"C:\Users\claud\OneDrive\Documents\SALENTO\Semester 2\Modeling Ecology\Project\PugliaFishingEffort\Processed\ClippedHabitatMap_QGIS\Final\HabitatMap_EUNIS2007_Box_FINAL.shp"
if os.path.exists(habitat_path):
    try:
        habitat_gdf = gpd.read_file(habitat_path)
        print("\nLoaded habitat shapefile from:", habitat_path)
        print("Habitat CRS:", habitat_gdf.crs)
        print("Number of features:", len(habitat_gdf))
        print("Geometry types:", habitat_gdf.geom_type.unique())
        try:
            bounds = habitat_gdf.total_bounds
            print("Habitat bounds (minx, miny, maxx, maxy):", bounds)
        except Exception:
            pass

        # Ensure CRS is geographic (lon/lat) to match datasets (EPSG:4326)
        try:
            if habitat_gdf.crs is not None and habitat_gdf.crs.to_string() not in ("EPSG:4326", "CRS84"):
                habitat_gdf = habitat_gdf.to_crs(epsg=4326)
                print("Reprojected habitat to EPSG:4326 for overlay")
        except Exception as e:
            print("Warning: could not check/reproject CRS:", e)

        # Build or load a persistent color mapping for categories in 'allcomb'
        color_map_file = 'habitat_colors.json'
        color_map = {}
        if os.path.exists(color_map_file):
            try:
                with open(color_map_file, 'r', encoding='utf-8') as f:
                    color_map = json.load(f)
            except Exception:
                color_map = {}

        # Detect the category field case-insensitively (e.g., 'allcomb', 'ALLCOMB', 'AllComb')
        target_name = 'allcomb'
        cat_field = None
        lower_cols = {c.lower(): c for c in habitat_gdf.columns}
        if target_name in lower_cols:
            cat_field = lower_cols[target_name]

        # Fallback: choose first likely categorical/text column if allcomb not found
        if cat_field is None:
            possible = [c for c in habitat_gdf.columns if c.lower() not in ['geometry'] and habitat_gdf[c].dtype == 'object']
            if possible:
                cat_field = possible[0]
                print(f"Warning: column 'allcomb' not found. Using '{cat_field}' instead for coloring.")
            else:
                print("Warning: no suitable category/text column found; using uniform color.")

        # Legend title used for maps
        legend_title = "Marine Biome"

        if cat_field is not None and cat_field in habitat_gdf.columns:
            categories = sorted([str(v) for v in habitat_gdf[cat_field].dropna().unique().tolist()])
            print(f"Using habitat category field: {cat_field} with {len(categories)} categories")
            # assign new colors for any categories not already in map
            palette = plt.get_cmap('tab20')
            next_color_idx = 0
            for cat in categories:
                if cat not in color_map:
                    color_map[cat] = mcolors.to_hex(palette(next_color_idx % palette.N))
                    next_color_idx += 1
            # persist updated mapping
            try:
                with open(color_map_file, 'w', encoding='utf-8') as f:
                    json.dump(color_map, f, indent=2)
            except Exception as e:
                print("Warning: could not save color map:", e)

            # Attach color per feature
            habitat_gdf['_color'] = habitat_gdf[cat_field].astype(str).map(color_map).fillna('#cccccc')
        else:
            print("Warning: category field not available; using uniform color")
            habitat_gdf['_color'] = '#66c2a5'
            categories = []
            color_map = {}

        # Save a colored habitat-only map with legend
        fig, ax = plt.subplots(figsize=(10, 8))
        habitat_gdf.plot(ax=ax, color=habitat_gdf['_color'], edgecolor='black', linewidth=0.2, alpha=0.7)
        ax.set_title('Habitat Map by allcomb')
        # Legend patches (place BELOW the chart in two columns)
        if color_map:
            patches = [mpatches.Patch(color=color_map[c], label=c) for c in categories]
            fig.legend(
                handles=patches,
                title=legend_title,
                loc='lower center',
                bbox_to_anchor=(0.5, -0.02),
                ncol=2,
                fontsize='small',
                title_fontsize='small',
                frameon=True,
                columnspacing=0.8,
                labelspacing=0.4,
                borderaxespad=0.2,
            )
            # Make space at the bottom for the legend row
            fig.subplots_adjust(bottom=0.28)
        plt.xlabel('Longitude'); plt.ylabel('Latitude')
        plt.tight_layout()
        plt.savefig('habitat_map.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Habitat map saved as habitat_map.png")

        # Also save an overlay with total chlorophyll at t0
        try:
            fig2, ax2 = plt.subplots(figsize=(10, 8))
            chl_surface.isel(time=0).plot(ax=ax2, cmap='viridis', robust=True, cbar_kwargs={'label': 'Total Chl-a (mg/m³)'})
            habitat_gdf.plot(ax=ax2, color=habitat_gdf['_color'], edgecolor='black', linewidth=0.2, alpha=0.3)
            ax2.set_title('Total Chlorophyll (t0) with Habitat Overlay')
            # Legend BELOW in two columns
            if color_map:
                patches = [mpatches.Patch(color=color_map[c], label=c) for c in categories]
                fig2.legend(
                    handles=patches,
                    title=legend_title,
                    loc='lower center',
                    bbox_to_anchor=(0.5, -0.02),
                    ncol=2,
                    fontsize='small',
                    title_fontsize='small',
                    frameon=True,
                    columnspacing=0.8,
                    labelspacing=0.4,
                    borderaxespad=0.2,
                )
                fig2.subplots_adjust(bottom=0.28)
            plt.xlabel('Longitude'); plt.ylabel('Latitude')
            plt.tight_layout()
            plt.savefig('chlorophyll_habitat_overlay.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("Overlay saved as chlorophyll_habitat_overlay.png")
        except Exception as e:
            print("Warning: failed to create chlorophyll-habitat overlay:", e)
    except Exception as e:
        print("Warning: failed to load or plot habitat shapefile:", e)
        habitat_gdf = None
else:
    print("\nHabitat file not found at:", habitat_path)
    habitat_gdf = None

# Load fishing effort data (user-provided path)
fishing_path = r"C:\Users\claud\OneDrive\Documents\SALENTO\Semester 2\Modeling Ecology\Project\PugliaFishingEffort\Raw\FishingEffort\layer-activity-data-0\public-global-fishing-effort-v3.0.tif"
fishing_da = None
if os.path.exists(fishing_path):
    try:
        fishing_da = rioxarray.open_rasterio(fishing_path, masked=True)
        print("\nLoaded fishing effort GeoTIFF from:", fishing_path)
        print("Fishing effort dims:", fishing_da.dims)
        print("Fishing effort shape:", fishing_da.shape)
        try:
            print("Fishing CRS:", fishing_da.rio.crs)
            bounds = fishing_da.rio.bounds()
            print("Fishing bounds (minx, miny, maxx, maxy):", bounds)
        except Exception:
            pass

        # Ensure CRS is EPSG:4326 for overlay
        try:
            if fishing_da.rio.crs is not None and str(fishing_da.rio.crs) not in ("EPSG:4326", "CRS84"):
                fishing_da = fishing_da.rio.reproject("EPSG:4326")
                print("Reprojected fishing effort to EPSG:4326")
        except Exception as e:
            print("Warning: could not reproject fishing effort:", e)

        # Save a standalone fishing effort map
        try:
            fig_fish, ax_fish = plt.subplots(figsize=(10, 8))
            # If multi-band, select first band
            if 'band' in fishing_da.dims and fishing_da.sizes['band'] > 1:
                fishing_plot = fishing_da.isel(band=0).squeeze()
            else:
                fishing_plot = fishing_da.squeeze()

            fishing_plot.plot(ax=ax_fish, cmap='YlOrRd', robust=True, cbar_kwargs={'label': 'Fishing Hours'})
            ax_fish.set_title('Fishing Effort (Hours)')
            ax_fish.set_xlabel('Longitude')
            ax_fish.set_ylabel('Latitude')

            # Overlay habitat if available
            if 'habitat_gdf' in globals() and habitat_gdf is not None:
                habitat_gdf.boundary.plot(ax=ax_fish, edgecolor='black', linewidth=0.3)

            plt.tight_layout()
            plt.savefig('fishing_effort.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("Fishing effort map saved as fishing_effort.png")
        except Exception as e:
            print("Warning: failed to plot fishing effort:", e)
    except Exception as e:
        print("Warning: failed to load fishing effort GeoTIFF:", e)
        fishing_da = None
else:
    print("\nFishing effort file not found at:", fishing_path)
    fishing_da = None

# Ask for the correct paths (update these if needed)
print("\nPlease verify the following file paths:")
print("1. Currents data file")
print("2. Fishing effort data file")
print("3. Habitat shapefile")

# Set currents path (user-provided)
currents_path = r"C:\Users\claud\OneDrive\Documents\SALENTO\Semester 2\Modeling Ecology\Project\PugliaFishingEffort\Raw\CurrentData\cmems_mod_med_phy-cur_anfc_4.2km_P1M-m_1761297806777.nc"
currents = None
if os.path.exists(currents_path):
    try:
        currents = xr.open_dataset(currents_path)
        print("\nLoaded currents dataset from:", currents_path)
    except Exception as e:
        print("Warning: failed to open currents dataset:", e)
        currents = None
else:
    print("\nCurrents file not found at: ", currents_path)

# If currents loaded, process and plot; otherwise skip gracefully
if currents is None:
    print("\nSkipping currents processing because 'currents' is not available.")
else:
    # Debug information about the currents dataset
    print("\nDataset information:")
    print(currents)

    # Extract and print time coordinates
    if 'time' in currents.coords:
        print("\nTime coordinates in currents dataset:")
        print(currents.time.values)

    # Get u and v components (check availability)
    if 'uo' in currents and 'vo' in currents:
        u = currents['uo']  # Zonal current component
        v = currents['vo']  # Meridional current component

        print("\nU variable information:")
        print(u)
        print("\nV variable information:")
        print(v)

        #analysis start
        # Print debug information about the current data
        if 'time' in u.coords and 'time' in v.coords:
            print("U time coordinates:", u.time.values)
            print("V time coordinates:", v.time.values)

        # Align/select compatible time slices
        # Strategy: pick the intersection of time coordinates if they differ
        try:
            u_times = u.coords['time'].values if 'time' in u.coords else None
            v_times = v.coords['time'].values if 'time' in v.coords else None
            if u_times is not None and v_times is not None and len(u_times) != len(v_times):
                # find common times
                common_times = np.intersect1d(u_times, v_times)
                if common_times.size == 0:
                    print("No common time steps between u and v. Selecting first timestep of each.")
                    u_single = u.isel(time=0)
                    v_single = v.isel(time=0)
                else:
                    # select the first common time
                    t0 = common_times[0]
                    u_single = u.sel(time=t0)
                    v_single = v.sel(time=t0)
            else:
                # same length or no time dimension: select first timestep if present
                if 'time' in u.coords:
                    u_single = u.isel(time=0)
                else:
                    u_single = u
                if 'time' in v.coords:
                    v_single = v.isel(time=0)
                else:
                    v_single = v
        except Exception as e:
            print("Warning while aligning times:", e)
            u_single = u.isel(time=0) if 'time' in u.coords else u
            v_single = v.isel(time=0) if 'time' in v.coords else v

        print("Single timestep shapes - u:", u_single.shape, "v:", v_single.shape)

        # Now align the single timestep data
        u_aligned = u_single.broadcast_like(v_single)
        v_aligned = v_single.broadcast_like(u_single)

        print("Aligned shapes:", u_aligned.shape, v_aligned.shape)

        # Calculate speed as an xarray.DataArray (keeps coords/dims)
        speed_da = np.sqrt(u_aligned**2 + v_aligned**2)
        print("Speed DataArray dims:", speed_da.dims, "shape:", speed_da.shape)

        # If there's a depth dimension, select the surface (shallowest) layer
        if 'depth' in speed_da.dims:
            # choose index 0 as the surface (modify if another depth index is preferred)
            speed2d = speed_da.isel(depth=0)
            print("Selected depth index 0 for surface; new shape:", speed2d.shape)
        else:
            speed2d = speed_da

        # If speed2d still has unexpected dims (e.g., extra singleton dims), squeeze
        speed2d = speed2d.squeeze()

        # Create figure with specific size
        plt.figure(figsize=(12, 8))

        # Use xarray's plotting which understands coords (preferred)
        try:
            speed2d.plot(cmap='viridis', robust=True)
        except Exception:
            # Fallback to pcolormesh: ensure we pass 2D arrays for X, Y
            try:
                X = u_aligned.longitude.values
                Y = u_aligned.latitude.values
                C = speed2d.values
                plt.pcolormesh(X, Y, C, cmap='viridis', shading='auto')
            except Exception as e:
                print("Failed to plot with fallback pcolormesh:", e)
                raise

        plt.title('Surface Current Speed - Selected Timestep')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        # Overlay habitat polygons if available
        try:
            if 'habitat_gdf' in globals() and habitat_gdf is not None:
                habitat_gdf.plot(ax=plt.gca(), color=(habitat_gdf['_color'] if '_color' in habitat_gdf.columns else '#66c2a5'),
                                 edgecolor='black', linewidth=0.2, alpha=0.3)
                if 'color_map' in globals() and color_map:
                    patches = [mpatches.Patch(color=color_map[c], label=c) for c in categories]
                    fig_cs = plt.gcf()
                    fig_cs.legend(
                        handles=patches,
                        title=(legend_title if 'legend_title' in globals() else 'Marine Biome'),
                        loc='lower center',
                        bbox_to_anchor=(0.5, -0.02),
                        ncol=2,
                        fontsize='small',
                        title_fontsize='small',
                        frameon=True,
                        columnspacing=0.8,
                        labelspacing=0.4,
                        borderaxespad=0.2,
                    )
                    fig_cs.subplots_adjust(bottom=0.28)
        except Exception as e:
            print("Warning: failed to overlay habitat on speed plot:", e)

        # Save the plot (legend outside requires bbox_inches='tight')
        plt.savefig('current_speed.png', dpi=300, bbox_inches='tight')
        print("Plot saved as current_speed.png")

        # Show and then close
        plt.show()
        plt.close()
    else:
        print("Currents dataset does not contain 'uo' and 'vo' variables.")

# Close datasets to free resources
try:
    chl_dataset.close()
except Exception:
    pass
try:
    if 'currents' in locals() and currents is not None:
        currents.close()
except Exception:
    pass