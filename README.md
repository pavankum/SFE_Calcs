# SFE_Calcs

create_top_gro.ipynb - This notebook creates inputs for the SFE run. You can change the force field and charge model as you want, highlighted the code block for MH's purpose. Code contributions from Melissa Plata, Meghan Osato, Pavan Behara.

sfe_calcs.yml - you can install the conda environment to run this book via `mamba env create -f sfe_calcs.yml`, add any extra packages needed

For the equlibrium calculations the 0001.top and 0001.gro in `run_eq.sh` has to be replaced with the appropriate directory, can use
```
for i in $(seq -f "%04g" 1 115);do cd $i;sed -i "s/0001/$i/g" run_eq.sh;sbatch run_eq.sh;cd ../;done
```

Post-processing scripts from Meghan will be added soon.
