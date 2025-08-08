# Data Acquisition Instructions

This project depends on three publicly available astronomical catalogs. Below are the exact sources, the minimal columns the analysis requires, and two supported retrieval methods: (1) using the **VizieR web interface** and (2) programmatically via **Astroquery** in Python.  
The file `example_data/all_hmxb.tsv` included in this repository already contains the HMXB candidate list, so no external download is needed for those.

## Catalogs & VizieR Identifiers

1. **Muno X-ray Catalog** (Chandra Galactic Center)  
   VizieR ID: `J/ApJS/181/110`  
   Needed columns: `RAJ2000`, `DEJ2000`, `ePos` (positional uncertainty in arcseconds)

2. **GALACTICNUCLEUS Infrared Catalog**  
   VizieR ID: `J/A+A/653/A133`  
   Needed columns: `RAJ2000`, `DEJ2000`, `Ksmag`

3. **VVV Infrared Catalog**  
   VizieR ID: `II/348`  
   Needed columns: `RAJ2000`, `DEJ2000`, and the appropriate Ks-band magnitude (e.g., `Ksmag3`).

## Option A: Manual via VizieR Web

1. Go to https://vizier.u-strasbg.fr  
2. Enter the catalog identifier (e.g., `J/ApJS/181/110`) in the "By catalog identifier" box.  
3. Click the catalog, then:
   - Use the interface to select the required columns.
   - Export/save the table in **Tab-separated (TSV)** or **ASCII** format.
   - Rename or ensure the columns match expectations:

## Option B: Programmatic (Python + Astroquery)

Install dependencies if needed:

```bash
pip install astroquery pandas numpy
```

Minimal example to fetch and clean:

```python
from astroquery.vizier import Vizier
import pandas as pd

Vizier.ROW_LIMIT = -1  # fetch full catalog

def fetch_catalog(catalog_id):
    viz = Vizier(columns=['*'])
    table = viz.get_catalogs(catalog_id)[0]
    return table.to_pandas()

# 1. Muno
muno = fetch_catalog('J/ApJS/181/110')
muno['RAJ2000'] = pd.to_numeric(muno['RAJ2000'], errors='coerce')
muno['DEJ2000'] = pd.to_numeric(muno['DEJ2000'], errors='coerce')
muno['ePos'] = pd.to_numeric(muno['ePos'], errors='coerce')
muno.to_csv('muno_catalog_clean.tsv', sep='\t', index=False)

# 2. GALACTICNUCLEUS
gn = fetch_catalog('J/A+A/653/A133')
gn['RAJ2000'] = pd.to_numeric(gn['RAJ2000'], errors='coerce')
gn['DEJ2000'] = pd.to_numeric(gn['DEJ2000'], errors='coerce')
gn['Ksmag'] = pd.to_numeric(gn['Ksmag'], errors='coerce')
gn.to_csv('gn_catalog_clean.tsv', sep='\t', index=False)

# 3. VVV
vvv = fetch_catalog('II/376')
vvv['RAJ2000'] = pd.to_numeric(vvv['RAJ2000'], errors='coerce')
vvv['DEJ2000'] = pd.to_numeric(vvv['DEJ2000'], errors='coerce')
# Replace 'Ks1ap3' below with the actual Ks column you choose
vvv['Ksmag3'] = pd.to_numeric(vvv['Ks1ap3'], errors='coerce')
vvv.to_csv('vvv_catalog_clean.tsv', sep='\t', index=False)
```

## Column Expectations for the Notebook

- **Muno (X-ray)**: `RAJ2000`, `DEJ2000`, `ePos`
- **GALACTICNUCLEUS (IR)**: `RAJ2000`, `DEJ2000`, `Ksmag`
- **VVV (IR)**: `RAJ2000`, `DEJ2000`, `Ksmag3`
- **HMXB candidates**: Provided in `all_hmxb.tsv` with:
  - `name`, `ra`, `dec`, `err_ellipse_r0`, `err_ellipse_r1`, `err_ellipse_ang`

## Attribution

> Catalog data retrieved from VizieR: Muno X-ray (`J/ApJS/181/110`), GALACTICNUCLEUS (`J/A+A/653/A133`), and VVV (`II/348`).

