import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from fontTools.ttLib.woff2 import bboxFormat

root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures'

# Coordinates (latitude, longitude, name)
coordinates = [
    (30.55197, -90.81396, 'EX'),   # WA, USA
    (30.52860, -90.76153, 'EY'),   # WA, USA
    (30.56281, -90.77427, 'IY'),   # WA, USA
    (30.56290, -90.77451, 'IX'),   # WA, USA
    #(46.46170, -119.41768, 'VA'),   # WA, USA
]

# Create a new figure and axes object
fig = plt.figure(figsize=(4, 4))
ax = plt.axes(projection=ccrs.PlateCarree())

# Add satellite imagery background using Cartopy's OSM tiles
img = cimgt.QuadtreeTiles()  # This will add OSM satellite tiles
ax.add_image(img, 12)

# Add borders, coastlines, and other features
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.COASTLINE)

# Plot the coordinates as markers
for lat, lon, name in coordinates:
    ax.plot(lon, lat, marker='^', color='blue', markersize=14, markeredgecolor='white', markeredgewidth=0.2, transform=ccrs.PlateCarree())
    if name== 'IX':
        ax.text(lon-0.004, lat, f'{name}', color='white', fontsize=15, ha='left', transform=ccrs.PlateCarree())
    else:
        ax.text(lon+0.001, lat, f'{name}', color='white', fontsize=15, ha='left', transform=ccrs.PlateCarree())

# Set map extent to focus on the region of interest (WA, USA)
ax.set_extent([-90.8275, -90.7525, 30.51, 30.58], crs=ccrs.PlateCarree())

# Add gridlines
gridlines = ax.gridlines(draw_labels=True, linewidth=0.5)
gridlines.top_labels = False
gridlines.right_labels = False

# Display the plot with title
plt.title('(d) Livingston site', loc='left')
fig.savefig(f"{root_savefig}/MAP_HAN_BRS_4Paper.png", bbox_inches='tight', dpi=300)
plt.show()


