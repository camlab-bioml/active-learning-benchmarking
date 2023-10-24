# Active learning benchmarking

To reproduce the analysis in this repository and the associated paper follow these steps:

1. Install `pipenv` using `pip install pipenv`
2. Clone this repository `git clone git@github.com:camlab-bioml/active-learning-benchmarking.git`
3. Move into the folder that was created as part of the previous step: `cd active-learning-benchmarking`
4. Install the environment required for this analysis to work using `pipenv install`
5. Create a snakemake profile by following these instructions: https://github.com/Snakemake-Profiles/slurm
6. Inside the root directory of this repository run the following commands
```
module load singularity
snakemake --profile [your_snakemake_profile] -j [n_jobs] --rerun-incomplete
```

where
- [your_snakemake_profile] is the name of the snakemake profile you generated in step 5
- [n_jobs] is the total number of jobs you would like to be submitted at a time
