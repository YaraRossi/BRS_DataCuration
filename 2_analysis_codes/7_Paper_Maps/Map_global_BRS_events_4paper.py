import pygmt
import pandas as pd

root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures'

# Define earthquake(s) locations
earthquakes = [
    (-74.454, -15.828, "Mw 7.2", 'CM'),  # Atiquipa
    (-130.612, 51.609, "Mw 6.5", 'RM'),  # Port McNeill
    (-119.383, 46.269, "Ml 1.8", 'LM')  # Explosion
]

path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

events = f'{path}/event_times_LIV.txt'  # Replace with your actual file path
data = pd.read_csv(events)


# Define seismic stations (global examples)
stations = [
    (-119.40785307788653, 46.45511825164643, 'HAN', 'black'),   # WA, USA
    #(-118.12602285654852, 34.136726696087834, 'CAL'),  # CA, USA
    #(-116.492470, 33.582325, 'BSPF', 'blue'),                   # CA, USA
    (-90.77469251757329, 30.562775410180993, 'LIV', 'black'),   # LA, USA
    #(11.275483826919913, 48.16291747145096, 'ROMY', 'purple4')   # Europe
    #(12.8780,49.1444, 'G-RING') # Germany
]

# Create PyGMT figure
fig = pygmt.Figure()

# Global map using Robinson projection
#fig.basemap(region="g", projection="R-0/10i", frame=["af"])
with pygmt.config(FONT_ANNOT_PRIMARY="20p", FONT_LABEL="20p",MAP_FRAME_WIDTH="2c"):
    fig.basemap(region="g", projection='N10i', frame=["af"]) # projection='N10i'
    #fig.basemap(region="g", projection='H12c', frame=["af","+tb)"]) # projection='N10i'


# add title
fig.text(text="(b)",position="TL",
    offset="-0.0c/-0.1c",justify="TL",font="20p,Helvetica")

# Add coastlines, water, and borders
fig.coast(land="gray50", water="white", shorelines="1/0.5p,black", borders=[1])

# plot all events as black stars
for _, data_event in data.iterrows():
    if 'tornado' in data_event[' type']:
        continue
    if 'eruption' in data_event[' type']:
        continue
    source_latitude = float(data_event[' latitude'])
    source_longitude = float(data_event[' longitude'])

    fig.plot(x=source_longitude, y=source_latitude, style="a0.8c", fill="yellow", pen="black")


# Plot all earthquakes as red stars
'''for lon, lat, label, justify in earthquakes:
    fig.plot(x=lon, y=lat, style="a0.8c", fill="orange", pen="black")
    fig.text(x=lon+5, y=lat + 10, text=label, font="20p,Helvetica-Bold,black", justify=justify)
    fig.text(x=lon+5, y=lat + 10, text=label, font="20p,Helvetica-Bold, orange", justify=justify)'''

# Plot all stations as blue triangles
for lon, lat, label, color in stations:
    if label == 'BSPF':
        fig.text(x=lon-5, y=lat, text=label, font=f"20p,Helvetica-Bold,{color}", justify="RM")
        fig.plot(x=lon, y=lat, style="t0.8c", fill=color, pen="black")
    else:
        fig.text(x=lon+5, y=lat, text=label, font=f"20p,Helvetica-Bold,{color}", justify="LM")
        fig.plot(x=lon, y=lat, style="t0.8c", fill=color, pen="black")

'''# Use an inset to place the title outside the main map.
# The inset creates a small panel that we can control independently.
with fig.inset(position="jTL+w1i", margin="0i", box="+gwhite+p1p,white"):
    # Inset coordinates are relative to the inset region.
    # Here we put the text at (x=0, y=0) which will be at the left‐most, bottom of the inset box.
    fig.text(x=0, y=0, text="(a)", font="20p,Helvetica,black", justify="LM")'''


# Save and optionally show the figure
fig.savefig(f"{root_savefig}/MAP_Global_BRS_4paper.png")

fig.show()