## Run the pipeline from the workflow

```bash
nohup snakemake \
  --configfile ../config/config-glicid-scratch.yml \
  --profile profiles/default/ \
  >> snake.out &
```

## Generate the models summary report

```bash
nohup snakemake \
  --configfile ../config/config-glicid-scratch.yml \
  --profile profiles/default/ \
  --report /scratch/nautilus/users/sanglier-l@univ-nantes.fr/Midasbolism/pipeline/results/report.zip  \
  >> snake.out &
```
