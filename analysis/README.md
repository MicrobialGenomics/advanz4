### README

Main directives:
- avoid exposing data openly in public using github.
- Use own branch for code development
- Do not modify master branch output
- Put data sources (including rds, csv metadata and related files in /data directory). data/* should be in .gitignore

General Recommendations

Instructions for code Organization

000-100 - Initial Data Setup
  002 - Data Transfer to s3
  004 - Metadata Curation
  006 - Data Curation, using sequence analysis results and curated metadata - produces mre object
  008 - MRE analysis population
  020 - Create Report using metaExplorer
  030 - Analysis Plan

100-200 - Standard non-project oriented Analysis


200-999 - Project oriented analysis


Use of workflowr package to render and manage code.
