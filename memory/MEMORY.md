# Memory Index

- [LV Scar Segmentation Project](project_lv_scar.md) - Current project context for `lv-scar-segmentation`: repo layout, active notebooks, results folders, data utilities, and working conventions
- [Polar Map Bull's Eye Pipeline](polar_map_pipeline.md) - 7-stage pipeline: 3D LV → 2D bull's eye (64×128); long axis, PCA septal ref, cylindrical projection, AHA-17 overlay; anatomical alignment (`_compute_transform`), pre-flattening diagnostic (sec 4a), sex-stratified + ring-density outputs
- [ISOMAP Shape Analysis Pipeline](isomap_pipeline.md) - Per-fragment intrinsic shape analysis: connected components → noise filter → ISOMAP embedding → normalised 2D shape; contrasts with polar map (where vs. what shape)
