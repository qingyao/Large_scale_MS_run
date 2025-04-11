0. `download_metadata.py` on all pxd IDs which pass the LLM model
1. `parse_metadata_for_filegroup.py` parse filenames into groups. PXD with readable file names will be grouped into `pxd_group_fnames.yaml` 
2. `parse_pxd_taxID.py` parse taxIDs. Check if exists in STRING and map taxIDs for the missing taxIDs to STRING. Write `pxdID_taxID.yaml`
3. `download_files.py` download raw files
4. `write_manifest_workflow.py` write manifest and workflow (depends on file availability, use `pxd_group_fnames.yaml`, manual editing when multiple taxIDs present)
5. `./fetch_fasta.bash $year_month_folder $string_ver $paxdb_ver` download fasta and `~/fragpipe_decoy_generation/run_decoy.bash` run DecoyDatabase to `fasta+decoy` folder
6. `run_fragpipe_wrapper.sh` runs fragpipe, logs the CPU/RAM usage with `monitor_pid.py`
7. `rewrite_scope_label.py` check if data scope is ok, finalize the label name, tissue.
8. `compute_abu.py` 
9. `compute_score.py`
10. `finalize_metadata.py` put all output files in converted/ folder into a table.
11. clean up the space, raw data, keep only necessary files.