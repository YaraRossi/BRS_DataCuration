from pathlib import Path
import re
from datetime import datetime

# -------- CONFIG --------
base_dir = Path("/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/zenodo_upload/miniseed_LIGO")  # CHANGE THIS
dry_run = False  # set to False to actually rename
# ------------------------

month_map = {
    "jan": 1, "january": 1,
    "feb": 2, "february": 2,
    "mar": 3, "march": 3,
    "apr": 4, "april": 4,
    "may": 5,
    "jun": 6, "june": 6,
    "jul": 7, "july": 7,
    "aug": 8, "august": 8,
    "sep": 9, "sept": 9, "september": 9,
    "oct": 10, "october": 10,
    "nov": 11, "november": 11,
    "dec": 12, "december": 12,
}

# Pattern: MonthDay.Year_optionalAddon
pattern = re.compile(
    r"^(?P<month>[A-Za-z]+)(?P<day>\d{1,2})\.(?P<year>\d{4})(?P<addon>[_\-].+)?$"
)

for folder in base_dir.iterdir():
    if not folder.is_dir():
        continue

    match = pattern.match(folder.name)
    if not match:
        print(f"Skipping (no match): {folder.name}")
        continue

    month_str = match.group("month").lower()
    day = int(match.group("day"))
    year = int(match.group("year"))
    addon = match.group("addon") or ""  # keep the underscore/dash if present

    if month_str not in month_map:
        print(f"Skipping (unknown month): {folder.name}")
        continue

    month = month_map[month_str]

    try:
        date = datetime(year, month, day)
    except ValueError as e:
        print(f"Skipping (invalid date): {folder.name} → {e}")
        continue

    new_name = date.strftime("%Y-%m-%d") + addon
    new_path = folder.parent / new_name

    print(f"{folder.name}  →  {new_name}")

    if not dry_run:
        folder.rename(new_path)
