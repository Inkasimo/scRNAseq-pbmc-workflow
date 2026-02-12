# Results layout

```text
results/
├── qc/
│   ├── fastqc/
│   │   └── raw/
│   │       ├── donor1/
│   │       ├── donor2/
│   │       ├── donor3/
│   │       └── donor4/
│   └── multiqc/
│       └── raw/
│           ├── multiqc_report.html
│           └── multiqc_report_data/
├── alignment/
│   └── starsolo/
│       └── raw/
│           ├── donor1/
│           ├── donor2/
│           ├── donor3/
│           └── donor4/
├── downstream/
│   ├── seurat/
│   │   └── untrimmed/
│   │       ├── donor1/
│   │       ├── donor2/
│   │       ├── donor3/
│   │       └── donor4/
│   ├── deg_and_tost/
│   │   └── untrimmed/
│   │       └── deg_and_tost/
│   └── networks/
│       └── untrimmed/
│           ├── consensus/
│           ├── per_donor/
│           ├── plots/
│           └── tables/
└── logs/
    ├── qc/
    ├── alignment/
    ├── download/
    └── ref/


```
