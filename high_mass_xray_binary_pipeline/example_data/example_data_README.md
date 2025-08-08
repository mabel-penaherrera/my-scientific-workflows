# Example Data for HMXB Pipeline

These are mostly synthetic stub versions of the real input catalogs used in the HMXB cross-matching analysis. They are minimal, small, and intended to allow the notebook to execute end-to-end for demonstration purposes—they do **NOT** produce science-grade results by themselves.

## Included files

- `muno3.csv`: simplified X-ray source catalog (stub).  
- `1.tsv` through `7.tsv`: simplified GALACTICNUCLEUS region files (semicolon-delimited, stubs).  
- `VVV_DR2_CDS.csv`: simplified near-infrared catalog (stub).  
- `all_hmxb.tsv`: **authentic** HMXB candidate list used in the analysis. This file is real (not synthetic) and contains the five candidates referenced throughout the notebook; no external download is required for it.

## Usage

Run `hmxb_crossmatching.ipynb` to exercise the pipeline. The stubs let you verify the code flow and plotting. To perform the actual scientific analysis, replace the stub catalogs with the full real catalogs (see `GET_DATA.md` for retrieval instructions), ensuring the notebook’s paths/column names are adjusted if needed.

## Notes

- The `all_hmxb.tsv` file supplies the genuine HMXB candidates; all other files are placeholders for demonstration.  
- Real catalogs (Muno, GALACTICNUCLEUS, VVV) are needed to reproduce published results. Refer to `GET_DATA.md` for how to obtain them.  
- If you replace any stub with a full catalog, clean numeric columns (e.g., magnitudes) with `pd.to_numeric(..., errors="coerce")` and verify coordinate units (expected in degrees).
