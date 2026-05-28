import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures'

# Coordinates (latitude, longitude, name)
coordinates = [
    (46.48424, -119.43821, 'EX'),   # WA, USA
    (46.43398, -119.44972, 'EY'),   # WA, USA
    (46.45517, -119.40764, 'IY'),   # WA, USA
    (46.46170, -119.41768, 'VA'),   # WA, USA
]

# Create a new figure and axes object
fig = plt.figure(figsize=(4.5, 4.5))
ax = plt.axes(projection=ccrs.PlateCarree())

# Add satellite imagery background using Cartopy's OSM tiles
img = cimgt.QuadtreeTiles()  # This will add OSM satellite tiles
ax.add_image(img, 12)

# Add borders, coastlines, and other features
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.COASTLINE)

# Add gridlines
gridlines = ax.gridlines(draw_labels=True, linewidth=0.5)
gridlines.top_labels = False
gridlines.right_labels = False

# Plot the coordinates as markers
for lat, lon, name in coordinates:
    ax.plot(lon, lat, marker='^', color='blue', markersize=14, markeredgecolor='black', markeredgewidth=0.2, transform=ccrs.PlateCarree())
    if 'E' not in name:
        ax.plot(lon, lat, marker='^', color='white', markersize=14, markeredgecolor='blue', markeredgewidth=1.5, transform=ccrs.PlateCarree())
    ax.text(lon+0.001, lat, f'{name}', color='white', fontsize=11, ha='left', transform=ccrs.PlateCarree())

# Set map extent to focus on the region of interest (WA, USA)
ax.set_extent([-119.465, -119.39, 46.4225, 46.4925], crs=ccrs.PlateCarree())



# Display the plot with title
plt.title('(c) Hanford site', loc='left', fontsize=11)
fig.savefig(f"{root_savefig}/MAP_LIV_BRS_4Paper.png", bbox_inches='tight', dpi=300)
plt.show()


