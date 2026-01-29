from pathlib import Path
import re

# -------- CONFIG --------
project_dir = Path("/BRS_DataCuration")
extensions = [".py", ".txt"]  # files to process
dry_run = False  # True = just print changes, False = actually overwrite
# ------------------------

# Match old folder pattern: MonthDay.Year[_addon]
pattern = re.compile(r"\b([A-Za-z]+)(\d{1,2})\.(\d{4})([_\-]\w+)?\b")

month_map = {
    "jan": "01", "january": "01",
    "feb": "02", "february": "02",
    "mar": "03", "march": "03",
    "apr": "04", "april": "04",
    "may": "05",
    "jun": "06", "june": "06",
    "jul": "07", "july": "07",
    "aug": "08", "august": "08",
    "sep": "09", "sept": "09", "september": "09",
    "oct": "10", "october": "10",
    "nov": "11", "november": "11",
    "dec": "12", "december": "12",
}

def replace_folder_name(match):
    month_str = match.group(1).lower()
    day = int(match.group(2))
    year = match.group(3)
    addon = match.group(4) or ""

    month = month_map.get(month_str)
    if not month:
        return match.group(0)  # leave unchanged if unknown month

    return f"{year}-{month}-{day:02d}{addon}"

# Iterate over all files with given extensions
for file_path in project_dir.rglob("*"):
    if not file_path.is_file() or file_path.suffix not in extensions:
        continue

    text = file_path.read_text()
    new_text = pattern.sub(replace_folder_name, text)

    if new_text != text:
        print(f"Changes in {file_path}")
        if not dry_run:
            file_path.write_text(new_text)
