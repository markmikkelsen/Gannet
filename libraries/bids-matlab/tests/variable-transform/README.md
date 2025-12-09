# variable-transform

Variable transformation for the BIDS stats model:
https://bids-standard.github.io/stats-models/

This repository contains validation data for several of the variable
transformation used in the BIDS stats model by bids matlab and pybids.

pybids transformation documentation:
https://docs.google.com/document/d/1uxN6vPWbC7ciAx2XWtT5Y-lBrdckZKpPdNUNpwRxHoU/edit#heading=h.kuzdziksbkpm

bids-matlab transformation documentation:
https://bids-matlab.readthedocs.io/en/latest/transformers.html

The transformations are split in 2 categories:

- [compute](./spec/compute.md)
- [munge](.spec/munge.md)

Each test within those categories is contained within one folder containting:

- input.tsv: defines the input data to be transformed
- output.tsv: shows the expected output of the transformation
- transform.json: specifies the transformation to be applied

## dev

See transformations.tsv and the online version:
[table of transformations](https://docs.google.com/spreadsheets/d/1_BykMB9Uybe9javHDwGtQVBC7jqi36eTGpibJsGC8t4/edit#gid=0)
